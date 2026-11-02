/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2026-     Equinor ASA
//
//  ResInsight is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ResInsight is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.
//
//  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
//  for more details.
//
/////////////////////////////////////////////////////////////////////////////////

#include "RicApproximateWellPathTools.h"

#include "RiaColorTables.h"
#include "RiaLogging.h"

#include "CompletionExportCommands/RicTransmissibilityCalculator.h"

#include "RifRestartWellConnectionReader.h"

#include "RigCell.h"
#include "RigEclipseCaseData.h"
#include "RigHexIntersectionTools.h"
#include "RigMainGrid.h"
#include "Well/RigSimWellData.h"
#include "Well/RigWellLogExtractor.h"
#include "Well/RigWellPath.h"
#include "Well/RigWellPathGeometryTools.h"
#include "Well/RigWellPathIntersectionTools.h"
#include "Well/RigWellResultFrame.h"
#include "Well/RigWellResultPoint.h"

#include "RimEclipseCase.h"
#include "RimEclipseView.h"
#include "RimPerforationCollection.h"
#include "RimPerforationInterval.h"
#include "RimPointBasedWellPath.h"
#include "RimProject.h"
#include "RimSegmentCollection.h"
#include "RimSimWellInView.h"
#include "RimTools.h"
#include "RimValveTemplate.h"
#include "RimValveTemplateCollection.h"
#include "RimWellPathCollection.h"
#include "RimWellPathCompletionSettings.h"
#include "RimWellPathValve.h"

#include <QFileInfo>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <tuple>

namespace
{
using Connection = RifRestartWellConnectionReader::Connection;
using WellData   = RifRestartWellConnectionReader::Well;
using FaceType   = cvf::StructGridInterface::FaceType;

//--------------------------------------------------------------------------------------------------
/// The connections ResInsight imported, for results without a readable restart file. Diameter and skin are
/// not known.
//--------------------------------------------------------------------------------------------------
WellData wellFromResultFrame( const RigMainGrid& mainGrid, const RigSimWellData& simWell, size_t timeStep )
{
    WellData well;
    well.name = simWell.m_wellName;

    if ( !simWell.hasWellResult( timeStep ) ) return well;
    const RigWellResultFrame* frame = simWell.wellResultFrame( timeStep );
    if ( !frame ) return well;

    if ( auto head = frame->wellHeadOrStartCell(); head.isCell() && head.gridIndex() == 0 )
    {
        if ( auto ijk = mainGrid.ijkFromCellIndex( head.cellIndex() ) )
        {
            well.headI = static_cast<int>( ijk->i() );
            well.headJ = static_cast<int>( ijk->j() );
        }
    }

    std::set<size_t> cells;
    for ( const auto& point : frame->allResultPoints() )
    {
        if ( !point.isCell() ) continue;
        if ( point.gridIndex() != 0 )
        {
            well.isInLgr = true;
            continue;
        }
        if ( !cells.insert( point.cellIndex() ).second ) continue;

        auto ijk = mainGrid.ijkFromCellIndex( point.cellIndex() );
        if ( !ijk ) continue;

        Connection connection;
        connection.i                = static_cast<int>( ijk->i() );
        connection.j                = static_cast<int>( ijk->j() );
        connection.k                = static_cast<int>( ijk->k() );
        connection.isOpen           = point.isOpen();
        connection.connectionFactor = point.connectionFactor();
        well.connections.push_back( connection );
    }
    return well;
}

//--------------------------------------------------------------------------------------------------
/// Parameters along p1 -> p2 (0 at p1, 1 at p2) where the segment crosses the faces of the cell
//--------------------------------------------------------------------------------------------------
std::vector<double> faceCrossings( const cvf::Vec3d& p1, const cvf::Vec3d& p2, const std::array<cvf::Vec3d, 8>& hexCorners )
{
    std::vector<HexIntersectionInfo> intersections;
    RigHexIntersectionTools::lineHexCellIntersection( p1, p2, hexCorners, 0, &intersections );

    const cvf::Vec3d    segment       = p2 - p1;
    const double        lengthSquared = segment.lengthSquared();
    std::vector<double> parameters;
    if ( lengthSquared <= 0.0 ) return parameters;

    for ( const auto& intersection : intersections )
    {
        parameters.push_back( std::clamp( ( intersection.m_intersectionPoint - p1 ) * segment / lengthSquared, 0.0, 1.0 ) );
    }
    return parameters;
}

//--------------------------------------------------------------------------------------------------
/// A short and smooth path//--------------------------------------------------------------------------------------------------
/// A short and smooth path through the connected cells, in connection order: straight inside each cell, from
/// the face where it enters the cell to the face where it leaves it. The path minimises its length plus a penalty on
/// its bends. Consecutive cells that touch share the
/// crossing point, which lies on their common face (or the part of it that the cells share across a fault).
/// Consecutive cells that do not touch get an exit point in the first and an entry point in the second, on the faces
/// that point most towards the other cell, joined by an unperforated straight line. The first cell is entered through
/// the face opposite to the next cell, the last cell left through the face opposite to the previous one.
//--------------------------------------------------------------------------------------------------
struct FaceRoute
{
    std::vector<cvf::Vec3d>                nodes;
    std::vector<std::pair<size_t, size_t>> cellNodes; // entry and exit node of each cell
};

FaceRoute faceRoute( const RigMainGrid& mainGrid, const std::vector<size_t>& cells )
{
    constexpr double edgeMargin = 0.05; // crossing points are kept this far (in face parameters) from the face edges
    constexpr int    samples    = 9; // per face direction and refinement pass
    constexpr int    passes     = 6;
    constexpr double centring   = 0.02; // weight of the distance to the face centre, which settles otherwise free points
    constexpr double bending    = 1000.0; // length [m] that a 90 degree bend costs: (1 - cos) of the angle times this

    const std::array<FaceType, 6> allFaces = { FaceType::NEG_I, FaceType::POS_I, FaceType::NEG_J, FaceType::POS_J, FaceType::NEG_K, FaceType::POS_K };

    auto bilinear = []( const std::array<cvf::Vec3d, 4>& c, double u, double v )
    { return c[0] * ( ( 1 - u ) * ( 1 - v ) ) + c[1] * ( u * ( 1 - v ) ) + c[2] * ( u * v ) + c[3] * ( ( 1 - u ) * v ); };

    // The faces of cell 'from' that point towards 'towards'
    auto facesTowards = [&]( size_t from, const cvf::Vec3d& towards )
    {
        const RigCell&        cell   = mainGrid.cell( from );
        const cvf::Vec3d      centre = cell.center();
        std::vector<FaceType> faces;
        FaceType              best      = FaceType::POS_K;
        double                bestAlign = -std::numeric_limits<double>::max();
        cvf::Vec3d            direction = towards - centre;
        direction.normalize();
        for ( FaceType face : allFaces )
        {
            cvf::Vec3d toFace = cell.faceCenter( face ) - centre;
            if ( !toFace.normalize() ) continue;
            const double align = toFace * direction;
            if ( align > 0.0 ) faces.push_back( face );
            if ( align > bestAlign )
            {
                bestAlign = align;
                best      = face;
            }
        }
        return std::make_pair( faces, best );
    };

    struct Candidate
    {
        cvf::Vec3d point;
        FaceType   face;
        double     u, v;
    };
    struct Node
    {
        size_t                 cell; // the cell whose faces carry the node
        std::vector<FaceType>  faces;
        size_t                 sharedWith = std::numeric_limits<size_t>::max(); // the other cell at a shared face
        std::vector<Candidate> candidates;
    };

    auto isShared = [&]( const Node& node, const cvf::Vec3d& point )
    {
        if ( node.sharedWith == std::numeric_limits<size_t>::max() ) return true;
        const cvf::Vec3d inside = point + ( mainGrid.cell( node.sharedWith ).center() - point ) * 0.02;
        return RigHexIntersectionTools::isPointInCell( inside, mainGrid.cellCornerVertices( node.sharedWith ) );
    };

    auto sample = [&]( Node& node, FaceType face, double u0, double u1, double v0, double v1 )
    {
        const auto corners = mainGrid.cell( node.cell ).faceCorners( face );
        for ( int a = 0; a < samples; ++a )
        {
            for ( int b = 0; b < samples; ++b )
            {
                const double u     = u0 + ( u1 - u0 ) * a / ( samples - 1 );
                const double v     = v0 + ( v1 - v0 ) * b / ( samples - 1 );
                const auto   point = bilinear( corners, u, v );
                if ( isShared( node, point ) ) node.candidates.push_back( { point, face, u, v } );
            }
        }
    };

    FaceRoute route;
    if ( cells.empty() ) return route;

    std::vector<Node> nodes;
    auto              upperK = [&]( size_t cell, bool upper )
    {
        const RigCell& c   = mainGrid.cell( cell );
        const bool     neg = c.faceCenter( FaceType::NEG_K ).z() > c.faceCenter( FaceType::POS_K ).z();
        return ( neg == upper ) ? FaceType::NEG_K : FaceType::POS_K;
    };

    // Entry of the first cell
    {
        const FaceType face =
            cells.size() > 1 ? facesTowards( cells[0], mainGrid.cell( cells[0] ).center() * 2.0 - mainGrid.cell( cells[1] ).center() ).second
                             : upperK( cells[0], true );
        nodes.push_back( { cells[0], { face } } );
    }
    route.cellNodes.assign( cells.size(), { 0, 0 } );
    for ( size_t c = 0; c + 1 < cells.size(); ++c )
    {
        Node shared{ cells[c], facesTowards( cells[c], mainGrid.cell( cells[c + 1] ).center() ).first, cells[c + 1] };
        for ( FaceType face : shared.faces )
        {
            sample( shared, face, edgeMargin, 1.0 - edgeMargin, edgeMargin, 1.0 - edgeMargin );
        }
        route.cellNodes[c].second = nodes.size();
        if ( !shared.candidates.empty() )
        {
            route.cellNodes[c + 1].first = nodes.size();
            nodes.push_back( shared );
        }
        else
        {
            nodes.push_back( { cells[c], { facesTowards( cells[c], mainGrid.cell( cells[c + 1] ).center() ).second } } );
            route.cellNodes[c + 1].first = nodes.size();
            nodes.push_back( { cells[c + 1], { facesTowards( cells[c + 1], mainGrid.cell( cells[c] ).center() ).second } } );
        }
    }
    // Exit of the last cell
    {
        const size_t   last = cells.back();
        const FaceType face =
            cells.size() > 1
                ? facesTowards( last, mainGrid.cell( last ).center() * 2.0 - mainGrid.cell( cells[cells.size() - 2] ).center() ).second
                : upperK( last, false );
        route.cellNodes.back().second = nodes.size();
        nodes.push_back( { last, { face } } );
    }
    for ( auto& node : nodes )
    {
        if ( node.candidates.empty() )
        {
            for ( FaceType face : node.faces )
            {
                sample( node, face, edgeMargin, 1.0 - edgeMargin, edgeMargin, 1.0 - edgeMargin );
            }
        }
    }

    // The cell that the segment ending at each node runs through, if any (not between cells that do not touch)
    std::vector<std::optional<std::array<cvf::Vec3d, 8>>> segmentCell( nodes.size() );
    for ( size_t c = 0; c < cells.size(); ++c )
    {
        const auto [entry, exit] = route.cellNodes[c];
        if ( exit == entry + 1 ) segmentCell[exit] = mainGrid.cellCornerVertices( cells[c] );
    }
    // A segment between points on the faces of a cell must not leave it, which it can through the curved faces of thin
    // cells at grazing angles
    auto staysInCell = [&]( size_t node, const cvf::Vec3d& from, const cvf::Vec3d& to )
    {
        if ( !segmentCell[node] ) return true;
        for ( double t : faceCrossings( from, to, *segmentCell[node] ) )
        {
            if ( t > 0.002 && t < 0.998 ) return false;
        }
        return true;
    };

    // The best path over the candidates by dynamic programming over pairs of consecutive points, refined around the
    // chosen points in later passes. Besides the length, bends are penalised, so that the path keeps its direction
    // where the cells leave it free, e.g. between connections that do not touch.
    constexpr double maxCost = std::numeric_limits<double>::max();
    auto             bend    = []( const cvf::Vec3d& a, const cvf::Vec3d& b, const cvf::Vec3d& c )
    {
        cvf::Vec3d in  = b - a;
        cvf::Vec3d out = c - b;
        if ( !in.normalize() || !out.normalize() ) return 0.0;
        return 1.0 - in * out;
    };
    std::vector<size_t> chosen( nodes.size(), 0 );
    double              halfWidth = ( 1.0 - 2.0 * edgeMargin ) / ( samples - 1 );
    for ( int pass = 0; pass < passes; ++pass )
    {
        auto centreCost = [&]( size_t n, size_t a )
        {
            const auto& candidate = nodes[n].candidates[a];
            return centring * ( candidate.point - mainGrid.cell( nodes[n].cell ).faceCenter( candidate.face ) ).length();
        };

        // cost[n][a][b]: the best path ending with candidate b of node n-1 and candidate a of node n
        std::vector<std::vector<std::vector<double>>> cost( nodes.size() );
        std::vector<std::vector<std::vector<size_t>>> from( nodes.size() );
        for ( size_t n = 1; n < nodes.size(); ++n )
        {
            const auto& current  = nodes[n].candidates;
            const auto& previous = nodes[n - 1].candidates;
            cost[n].assign( current.size(), std::vector<double>( previous.size(), maxCost ) );
            from[n].assign( current.size(), std::vector<size_t>( previous.size(), 0 ) );

            std::vector<std::vector<char>> allowed( current.size(), std::vector<char>( previous.size(), 1 ) );
            bool                           anyAllowed = false;
            for ( size_t a = 0; a < current.size(); ++a )
            {
                for ( size_t b = 0; b < previous.size(); ++b )
                {
                    allowed[a][b] = staysInCell( n, previous[b].point, current[a].point );
                    anyAllowed    = anyAllowed || allowed[a][b];
                }
            }
            // No segment stays in the cell: take the best ones anyway
            if ( !anyAllowed )
            {
                for ( auto& row : allowed )
                {
                    std::fill( row.begin(), row.end(), 1 );
                }
            }

            for ( size_t a = 0; a < current.size(); ++a )
            {
                const double centre = centreCost( n, a );
                for ( size_t b = 0; b < previous.size(); ++b )
                {
                    if ( !allowed[a][b] ) continue;
                    const double segment = ( current[a].point - previous[b].point ).length() + centre;
                    if ( n == 1 )
                    {
                        cost[n][a][b] = centreCost( 0, b ) + segment;
                        continue;
                    }
                    for ( size_t c = 0; c < nodes[n - 2].candidates.size(); ++c )
                    {
                        if ( cost[n - 1][b][c] == maxCost ) continue;
                        const double total = cost[n - 1][b][c] + segment +
                                             bending * bend( nodes[n - 2].candidates[c].point, previous[b].point, current[a].point );
                        if ( total < cost[n][a][b] )
                        {
                            cost[n][a][b] = total;
                            from[n][a][b] = c;
                        }
                    }
                }
            }
        }
        if ( nodes.size() > 1 )
        {
            const size_t last = nodes.size() - 1;
            double       best = maxCost;
            for ( size_t a = 0; a < cost[last].size(); ++a )
            {
                for ( size_t b = 0; b < cost[last][a].size(); ++b )
                {
                    if ( cost[last][a][b] < best )
                    {
                        best             = cost[last][a][b];
                        chosen[last]     = a;
                        chosen[last - 1] = b;
                    }
                }
            }
            for ( size_t n = last; n > 1; --n )
            {
                chosen[n - 2] = from[n][chosen[n]][chosen[n - 1]];
            }
        }
        if ( pass + 1 == passes ) break;

        for ( size_t n = 0; n < nodes.size(); ++n )
        {
            const Candidate best = nodes[n].candidates[chosen[n]];
            nodes[n].candidates.clear();
            sample( nodes[n],
                    best.face,
                    std::max( edgeMargin, best.u - halfWidth ),
                    std::min( 1.0 - edgeMargin, best.u + halfWidth ),
                    std::max( edgeMargin, best.v - halfWidth ),
                    std::min( 1.0 - edgeMargin, best.v + halfWidth ) );
            nodes[n].candidates.push_back( best );
        }
        halfWidth *= 2.0 / ( samples - 1 );
    }

    for ( size_t n = 0; n < nodes.size(); ++n )
    {
        route.nodes.push_back( nodes[n].candidates[chosen[n]].point );
    }
    return route;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimWellPathCompletionSettings::WellType wellType( const WellData& well )
{
    using WT = RimWellPathCompletionSettings::WellType;
    switch ( well.preferredPhase )
    {
        case 1:
            return WT::OIL;
        case 2:
            return WT::WATER;
        case 3:
            return WT::GAS;
        case 4:
            return WT::LIQUID;
        default:
            break;
    }
    if ( well.wellType == 3 ) return WT::WATER;
    if ( well.wellType == 4 ) return WT::GAS;
    return WT::OIL;
}

//--------------------------------------------------------------------------------------------------
/// The segment with the given number
//--------------------------------------------------------------------------------------------------
const RifRestartWellConnectionReader::Segment* findSegment( const WellData& well, int number )
{
    if ( number >= 1 && static_cast<size_t>( number ) <= well.segments.size() && well.segments[number - 1].number == number )
    {
        return &well.segments[number - 1];
    }
    for ( const auto& segment : well.segments )
    {
        if ( segment.number == number ) return &segment;
    }
    return nullptr;
}

//--------------------------------------------------------------------------------------------------
/// The tubing segment that a segment drains to: the segment itself, or the first regular segment downstream of a device
//--------------------------------------------------------------------------------------------------
const RifRestartWellConnectionReader::Segment* tubingSegment( const WellData& well, int number )
{
    using Type    = RifRestartWellConnectionReader::Segment::Type;
    auto* segment = findSegment( well, number );
    for ( size_t guard = 0; segment && segment->type != Type::REGULAR && guard < well.segments.size(); ++guard )
    {
        segment = findSegment( well, segment->outlet );
    }
    return segment;
}

//--------------------------------------------------------------------------------------------------
/// The main branch of a multi-segment well: the tubing length and depth of its segment nodes
//--------------------------------------------------------------------------------------------------
struct MainBore
{
    std::vector<double> lengths;
    std::vector<double> depths;
    double              diameter  = 0.0;
    double              roughness = 0.0;

    bool empty() const { return lengths.empty(); }

    double depthAt( double length ) const
    {
        if ( lengths.size() == 1 ) return depths.front();
        size_t i          = std::upper_bound( lengths.begin(), lengths.end(), length ) - lengths.begin();
        i                 = std::clamp<size_t>( i, 1, lengths.size() - 1 );
        const double span = lengths[i] - lengths[i - 1];
        const double t    = span > 0.0 ? ( length - lengths[i - 1] ) / span : 0.0;
        return depths[i - 1] + t * ( depths[i] - depths[i - 1] );
    }
};

MainBore mainBore( const WellData& well )
{
    using Type = RifRestartWellConnectionReader::Segment::Type;

    std::vector<const RifRestartWellConnectionReader::Segment*> segments;
    for ( const auto& segment : well.segments )
    {
        if ( segment.branch == 1 && segment.type == Type::REGULAR ) segments.push_back( &segment );
    }
    std::sort( segments.begin(), segments.end(), []( auto* a, auto* b ) { return a->length < b->length; } );

    MainBore bore;
    for ( auto* segment : segments )
    {
        bore.lengths.push_back( segment->length );
        bore.depths.push_back( segment->depth );
        // The top segment has no tubing
        if ( segment->number != 1 && bore.diameter == 0.0 && segment->diameter > 0.0 )
        {
            bore.diameter  = segment->diameter;
            bore.roughness = segment->roughness;
        }
    }
    return bore;
}

//--------------------------------------------------------------------------------------------------
/// An ICD valve template with the flow coefficient of a valve segment (WSEGVALV) and the given area, reused when it exists
//--------------------------------------------------------------------------------------------------
RimValveTemplate* valveTemplateFor( const RifRestartWellConnectionReader::Segment& segment,
                                    double                                         area,
                                    RiaDefines::EclipseUnitSystem                  unitSystem,
                                    const QString&                                 wellName )
{
    auto* collection = RimTools::valveTemplateCollection();
    if ( !collection ) return nullptr;

    // Orifice diameter in mm (METRIC) or inches (FIELD), from the area in m2 or ft2
    const double  unitsPerLength = unitSystem == RiaDefines::EclipseUnitSystem::UNITS_FIELD ? 12.0 : 1000.0;
    const double  diameter       = std::sqrt( 4.0 * area / M_PI ) * unitsPerLength;
    const QString label =
        QString( "SIM_%1 Cv %2 d %3" ).arg( wellName ).arg( segment.valveFlowCoefficient, 0, 'g', 4 ).arg( diameter, 0, 'g', 5 );

    for ( auto* existing : collection->items() )
    {
        if ( existing->type() == RiaDefines::WellPathComponentType::ICD && existing->name().endsWith( label ) ) return existing;
    }

    auto* valveTemplate = new RimValveTemplate();
    valveTemplate->setType( RiaDefines::WellPathComponentType::ICD );
    valveTemplate->setUnitSystem( unitSystem );
    valveTemplate->setUserLabel( label );
    valveTemplate->setName( valveTemplate->fullLabel() );
    valveTemplate->setDefaultValuesFromUnits();
    valveTemplate->setOrificeDiameter( diameter );
    valveTemplate->setFlowCoefficient( segment.valveFlowCoefficient );
    collection->addItem( valveTemplate );
    collection->updateConnectedEditors();
    return valveTemplate;
}

//--------------------------------------------------------------------------------------------------
/// Perforations are kept this far from where the path enters and leaves a cell, so that a completion export does not see
/// the neighbour cells, also on other grids with nearly the same faces (e.g. refined grids)
//--------------------------------------------------------------------------------------------------
double perforationMargin( double length )
{
    return std::min( 0.05, 0.1 * length );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int lastTimeStepWithConnections( const RigSimWellData* simWell, size_t timeStepCount )
{
    if ( !simWell ) return static_cast<int>( timeStepCount ) - 1;

    for ( size_t t = timeStepCount; t-- > 0; )
    {
        if ( !simWell->hasWellResult( t ) ) continue;
        const auto* frame = simWell->wellResultFrame( t );
        if ( !frame ) continue;
        const auto points = frame->allResultPoints();
        if ( std::any_of( points.begin(), points.end(), []( const RigWellResultPoint& p ) { return p.isCell(); } ) )
        {
            return static_cast<int>( t );
        }
    }
    return -1;
}
} // namespace

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RicApproximateWellPathTools::wellPathName( const QString& simulationWellName )
{
    return "SIM_" + simulationWellName;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<RimWellPath*, QString> RicApproximateWellPathTools::createApproximateWellPath( RimEclipseCase* eclipseCase,
                                                                                             const QString&  simulationWellName,
                                                                                             int             timeStep,
                                                                                             bool            useMultiSegmentData )
{
    if ( !eclipseCase || !eclipseCase->eclipseCaseData() || !eclipseCase->eclipseCaseData()->mainGrid() )
    {
        return std::unexpected( "The case has no grid." );
    }
    RigEclipseCaseData* caseData = eclipseCase->eclipseCaseData();
    const RigMainGrid&  mainGrid = *caseData->mainGrid();

    RimWellPathCollection* wellPathCollection = RimWellPathCollection::instance();
    if ( !wellPathCollection ) return std::unexpected( "The project has no well path collection." );

    const RigSimWellData*        simWell = caseData->findSimWellData( simulationWellName );
    const std::vector<QDateTime> dates   = eclipseCase->timeStepDates();
    if ( dates.empty() ) return std::unexpected( "The case has no time steps." );

    if ( timeStep < 0 ) timeStep = lastTimeStepWithConnections( simWell, dates.size() );
    if ( timeStep < 0 ) return std::unexpected( QString( "Well '%1' has no connections at any time step." ).arg( simulationWellName ) );
    if ( static_cast<size_t>( timeStep ) >= dates.size() )
    {
        return std::unexpected( QString( "Time step %1 is out of range, the case has %2." ).arg( timeStep ).arg( dates.size() ) );
    }
    const QDate date = dates[timeStep].date();

    // The connections of the restart file, else those ResInsight imported
    std::optional<WellData> well;
    QStringList             warnings;
    {
        auto wells =
            RifRestartWellConnectionReader::readWells( RifRestartWellConnectionReader::restartFiles( eclipseCase->gridFileName() ), date );
        if ( wells )
        {
            for ( const auto& w : *wells )
            {
                if ( w.name == simulationWellName ) well = w;
            }
        }
        if ( !well && simWell )
        {
            warnings.push_back( QString( "No restart data for well '%1' at %2 (%3); the connection diameters and skin factors are "
                                         "defaulted." )
                                    .arg( simulationWellName )
                                    .arg( date.toString( Qt::ISODate ) )
                                    .arg( wells ? "the well is not in the restart file" : wells.error() ) );
            well = wellFromResultFrame( mainGrid, *simWell, static_cast<size_t>( timeStep ) );
        }
    }
    if ( !well ) return std::unexpected( QString( "The case has no well named '%1'." ).arg( simulationWellName ) );
    if ( well->isInLgr )
    {
        return std::unexpected( QString( "Well '%1' is completed in an LGR, which is not supported." ).arg( simulationWellName ) );
    }
    const MainBore bore = well->isMultiSegment && useMultiSegmentData ? mainBore( *well ) : MainBore();
    if ( well->isMultiSegment && useMultiSegmentData && bore.empty() )
    {
        warnings.push_back( QString( "Well '%1' is a multi-segment well without segment data in the restart file; only its connections "
                                     "are used." )
                                .arg( simulationWellName ) );
    }

    // Connections in cells of the main grid
    struct ConnectionCell
    {
        Connection connection;
        size_t     cellIndex;
    };
    std::vector<ConnectionCell> connectionCells;
    for ( const auto& connection : well->connections )
    {
        if ( connection.i < 0 || connection.j < 0 || connection.k < 0 || static_cast<size_t>( connection.i ) >= mainGrid.cellCountI() ||
             static_cast<size_t>( connection.j ) >= mainGrid.cellCountJ() || static_cast<size_t>( connection.k ) >= mainGrid.cellCountK() )
        {
            warnings.push_back( QString( "Connection (%1, %2, %3) is outside the grid and was skipped." )
                                    .arg( connection.i + 1 )
                                    .arg( connection.j + 1 )
                                    .arg( connection.k + 1 ) );
            continue;
        }
        const size_t cellIndex = mainGrid.cellIndexFromIJK( connection.i, connection.j, connection.k );
        if ( mainGrid.cell( cellIndex ).isInvalid() )
        {
            warnings.push_back( QString( "Connection (%1, %2, %3) is in a cell without geometry and was skipped." )
                                    .arg( connection.i + 1 )
                                    .arg( connection.j + 1 )
                                    .arg( connection.k + 1 ) );
            continue;
        }
        if ( !bore.empty() && connection.segment > 0 )
        {
            const auto* tubing = tubingSegment( *well, connection.segment );
            if ( tubing && tubing->branch != 1 )
            {
                warnings.push_back(
                    QString( "Connection (%1, %2, %3) is on lateral branch %4, which is not reconstructed, and was skipped." )
                        .arg( connection.i + 1 )
                        .arg( connection.j + 1 )
                        .arg( connection.k + 1 )
                        .arg( tubing->branch ) );
                continue;
            }
        }
        connectionCells.push_back( { connection, cellIndex } );
    }
    if ( connectionCells.empty() )
    {
        return std::unexpected(
            QString( "Well '%1' has no connections in the grid at %2." ).arg( simulationWellName ).arg( date.toString( Qt::ISODate ) ) );
    }

    // The well path, new or replacing a previous approximation
    const QString name     = wellPathName( simulationWellName );
    auto*         existing = wellPathCollection->wellPathByName( name );
    auto*         wellPath = dynamic_cast<RimPointBasedWellPath*>( existing );
    if ( existing && !wellPath )
    {
        return std::unexpected( QString( "The project has a well path named '%1' that is not an approximate well path." ).arg( name ) );
    }
    if ( !wellPath )
    {
        wellPath = new RimPointBasedWellPath();
        wellPath->setName( name );
        wellPath->setWellPathColor(
            RiaColorTables::editableWellPathsPaletteColors().cycledColor3f( wellPathCollection->allWellPaths().size() ) );
        wellPath->setUnitSystem( caseData->unitsType() );
        wellPathCollection->addWellPath( wellPath, false );
    }
    const double defaultDiameter = caseData->unitsType() == RiaDefines::EclipseUnitSystem::UNITS_FIELD ? 0.709 : 0.216;

    // The path: from sea level above the well head (or the first connected cell), then through the connected cells
    size_t headCell = connectionCells.front().cellIndex;
    if ( well->headI >= 0 && well->headJ >= 0 && static_cast<size_t>( well->headI ) < mainGrid.cellCountI() &&
         static_cast<size_t>( well->headJ ) < mainGrid.cellCountJ() )
    {
        for ( size_t k = 0; k < mainGrid.cellCountK(); ++k )
        {
            const size_t cellIndex = mainGrid.cellIndexFromIJK( well->headI, well->headJ, k );
            if ( !mainGrid.cell( cellIndex ).isInvalid() )
            {
                headCell = cellIndex;
                break;
            }
        }
    }
    const RigCell&   head    = mainGrid.cell( headCell );
    const cvf::Vec3d headTop = head.faceCenter( FaceType::NEG_K ).z() > head.faceCenter( FaceType::POS_K ).z()
                                   ? head.faceCenter( FaceType::NEG_K )
                                   : head.faceCenter( FaceType::POS_K );

    std::vector<cvf::Vec3d> points;
    points.emplace_back( headTop.x(), headTop.y(), std::max( 0.0, headTop.z() + 1.0 ) );
    points.push_back( headTop );

    std::vector<std::pair<size_t, size_t>> cellPointIndices; // entry and exit point of each cell
    // Connection intervals along the tubing closer than this follow each other (exported values are rounded)
    constexpr double                                      contiguityTolerance = 0.5;
    bool                                                  intervalsUsed       = false;
    std::vector<std::optional<std::pair<double, double>>> tubingInterval( connectionCells.size() ); // used connection intervals
    std::optional<size_t>                                 topPointIndex;
    bool                                                  measuredDepthCalibrated   = false;
    double                                                measuredDepthAtFirstPoint = 0.0;
    {
        std::vector<size_t> cells;
        for ( const auto& connectionCell : connectionCells )
        {
            cells.push_back( connectionCell.cellIndex );
        }
        // Multi-segment wells: where the connections are along the tubing (COMPSEGS) gives the measured depths of the path
        // when the intervals follow each other, as when they were made from the path in each cell. Otherwise, e.g. when
        // the intervals are where devices are, the measured depths are fitted to the depths of the tubing. (Fitting the path
        // itself to the tubing lengths and depths gave connection factors on refined grids further from those of the
        // true well paths in tests, as the lateral position of the path is still free.)
        if ( !bore.empty() && connectionCells.size() > 1 )
        {
            // An interval is used when it is after the previous used one along the tubing (not e.g. inside it)
            std::vector<bool> valid( connectionCells.size(), false );
            double            lastEnd    = -std::numeric_limits<double>::max();
            size_t            contiguous = 0;
            for ( size_t c = 0; c < connectionCells.size(); ++c )
            {
                const auto& connection = connectionCells[c].connection;
                valid[c] = connection.segmentEnd > connection.segmentStart && connection.segmentStart >= lastEnd - contiguityTolerance;
                if ( valid[c] ) lastEnd = connection.segmentEnd;
                if ( c > 0 && std::abs( connection.segmentStart - connectionCells[c - 1].connection.segmentEnd ) < contiguityTolerance )
                    contiguous++;
            }
            intervalsUsed = 2 * contiguous >= connectionCells.size() - 1;
            if ( !intervalsUsed )
            {
                warnings.push_back(
                    QString( "The connection intervals along the tubing of well '%1' do not follow each other (%2 of %3 do); "
                             "they are not used for the path." )
                        .arg( simulationWellName )
                        .arg( contiguous )
                        .arg( connectionCells.size() - 1 ) );
            }
            else if ( std::count( valid.begin(), valid.end(), false ) > 0 )
            {
                warnings.push_back(
                    QString( "%1 connection interval(s) along the tubing of well '%2' are not after the previous ones and are "
                             "not used for the path." )
                        .arg( std::count( valid.begin(), valid.end(), false ) )
                        .arg( simulationWellName ) );
            }

            for ( size_t c = 0; c < connectionCells.size() && intervalsUsed; ++c )
            {
                if ( valid[c] )
                    tubingInterval[c] = std::make_pair( connectionCells[c].connection.segmentStart, connectionCells[c].connection.segmentEnd );
            }
            intervalsUsed = intervalsUsed && valid.front();
        }

        const FaceRoute route = faceRoute( mainGrid, cells );

        // Straight down to the first cell when the well head is above it
        const auto headIjk  = mainGrid.ijkFromCellIndex( headCell );
        const auto firstIjk = mainGrid.ijkFromCellIndex( cells.front() );
        if ( headIjk && firstIjk && headIjk->i() == firstIjk->i() && headIjk->j() == firstIjk->j() )
        {
            const cvf::Vec3d& entry = route.nodes.front();
            points                  = { cvf::Vec3d( entry.x(), entry.y(), std::max( 0.0, headTop.z() + 1.0 ) ) };
        }

        // Multi-segment wells: the measured depths of the path follow the tubing. The measured depth where the path enters
        // the first cell is the start of its connection interval along the tubing, when the intervals are used. Otherwise it
        // is the one that best fits the depths of the path to the depths of the tubing nodes along the tubing.
        std::optional<double> entryLength;
        if ( intervalsUsed && route.nodes.size() > 1 )
        {
            const auto&  first         = connectionCells.front().connection;
            const double lengthInFirst = ( route.nodes[route.cellNodes.front().second] - route.nodes.front() ).length();
            const bool   firstIsFull   = connectionCells.size() > 1 &&
                                     std::abs( connectionCells[1].connection.segmentStart - first.segmentEnd ) < contiguityTolerance;
            entryLength = firstIsFull ? first.segmentEnd - lengthInFirst : first.segmentStart;
        }
        else if ( !bore.empty() && bore.lengths.size() > 2 && route.nodes.size() > 2 )
        {
            std::vector<double> along( 1, 0.0 );
            for ( size_t n = 1; n < route.nodes.size(); ++n )
            {
                along.push_back( along.back() + ( route.nodes[n] - route.nodes[n - 1] ).length() );
            }
            double bestMisfit = std::numeric_limits<double>::max();
            for ( double start = bore.lengths.front(); start <= bore.lengths.back(); start += 0.25 )
            {
                double misfit = 0.0;
                size_t count  = 0;
                for ( size_t n = 0; n < route.nodes.size(); ++n )
                {
                    if ( start + along[n] > bore.lengths.back() ) break;
                    const double difference = bore.depthAt( start + along[n] ) + route.nodes[n].z();
                    misfit += difference * difference;
                    count++;
                }
                if ( count >= 3 && misfit / count < bestMisfit )
                {
                    bestMisfit  = misfit / count;
                    entryLength = start;
                }
            }
        }

        // The path passes the top segment at its length and depth, and is vertical above it, with the measured depth at sea
        // level making up the difference
        if ( entryLength && route.nodes.size() > 1 )
        {
            const cvf::Vec3d& entry     = route.nodes.front();
            const double      topLength = bore.lengths.front();
            const double      topDepth  = bore.depths.front();
            const double      alongTop  = *entryLength - topLength;
            const double      belowTop  = -entry.z() - topDepth;

            cvf::Vec3d backwards = entry - route.nodes[1];
            backwards.z()        = 0.0;
            if ( !backwards.normalize() ) backwards = cvf::Vec3d( 1, 0, 0 );

            if ( alongTop > 0.0 && alongTop >= std::abs( belowTop ) )
            {
                const double     horizontal = std::sqrt( alongTop * alongTop - belowTop * belowTop );
                const cvf::Vec3d top        = entry + backwards * horizontal + cvf::Vec3d( 0, 0, belowTop );
                points                      = { cvf::Vec3d( top.x(), top.y(), 0.0 ), top };
                measuredDepthAtFirstPoint   = topLength - topDepth;
                topPointIndex               = 1;
                measuredDepthCalibrated     = true;
            }
            else
            {
                points                    = { cvf::Vec3d( entry.x(), entry.y(), 0.0 ) };
                measuredDepthAtFirstPoint = *entryLength + entry.z();
                warnings.push_back( QString( "The top segment of well '%1' is not within reach of the first connection along the tubing; "
                                             "the path above the first connection is vertical." )
                                        .arg( simulationWellName ) );
            }
            if ( measuredDepthAtFirstPoint < 0.0 ) measuredDepthAtFirstPoint = 0.0;
        }
        const size_t offset = points.size();
        points.insert( points.end(), route.nodes.begin(), route.nodes.end() );
        for ( const auto& [entry, exit] : route.cellNodes )
        {
            cellPointIndices.emplace_back( offset + entry, offset + exit );
        }
        // A little beyond the exit of the last cell
        cvf::Vec3d direction = points.back() - points[points.size() - 2];
        if ( direction.normalize() ) points.push_back( points.back() + direction );
    }

    // Measured depths along the path, or along the tubing where the connection intervals are known: the path then has the
    // tubing length at the top segment and where it enters and leaves the connected cells, and in between the lengths
    // along the path scaled to fit
    std::vector<double> measuredDepths = RigWellPathGeometryTools::calculateMeasuredDepth( points, measuredDepthAtFirstPoint );
    if ( intervalsUsed )
    {
        std::vector<std::pair<size_t, double>> anchors;
        if ( topPointIndex ) anchors.emplace_back( *topPointIndex, bore.lengths.front() );
        for ( size_t c = 0; c < connectionCells.size(); ++c )
        {
            if ( !tubingInterval[c] ) continue;
            anchors.emplace_back( cellPointIndices[c].first, tubingInterval[c]->first );
            anchors.emplace_back( cellPointIndices[c].second, tubingInterval[c]->second );
        }
        std::stable_sort( anchors.begin(), anchors.end(), []( const auto& a, const auto& b ) { return a.first < b.first; } );
        std::vector<std::pair<size_t, double>> increasing;
        for ( const auto& anchor : anchors )
        {
            if ( increasing.empty() || ( anchor.first > increasing.back().first && anchor.second > increasing.back().second ) )
            {
                increasing.push_back( anchor );
            }
        }

        const std::vector<double> alongPath = measuredDepths;
        for ( size_t i = 0; i < measuredDepths.size() && !increasing.empty(); ++i )
        {
            auto next =
                std::lower_bound( increasing.begin(), increasing.end(), i, []( const auto& a, size_t index ) { return a.first < index; } );
            if ( next == increasing.begin() )
            {
                measuredDepths[i] = alongPath[i] + increasing.front().second - alongPath[increasing.front().first];
            }
            else if ( next == increasing.end() )
            {
                measuredDepths[i] = alongPath[i] + increasing.back().second - alongPath[increasing.back().first];
            }
            else
            {
                const auto&  previous = *( next - 1 );
                const double span     = alongPath[next->first] - alongPath[previous.first];
                const double t        = span > 0.0 ? ( alongPath[i] - alongPath[previous.first] ) / span : 0.0;
                measuredDepths[i]     = previous.second + t * ( next->second - previous.second );
            }
        }
    }

    wellPath->setUnitSystem( caseData->unitsType() );
    wellPath->setTrajectoryPoints( points, measuredDepths );

    const RigWellPath* geometry = wellPath->wellPathGeometry();
    if ( !geometry ) return std::unexpected( "Failed to create the well path geometry." );
    const auto& mds = geometry->measuredDepths();

    auto* completionSettings = wellPath->completionSettings();
    completionSettings->setWellNameForExport( simulationWellName );
    if ( !well->groupName.isEmpty() ) completionSettings->setGroupName( well->groupName );
    completionSettings->setReferenceDepth( well->referenceDepth );
    completionSettings->setDrainageRadius( well->drainageRadius );
    completionSettings->setWellType( wellType( *well ) );
    completionSettings->setAllowWellCrossFlow( well->allowCrossFlow );

    // One perforation per connection, over the path in the cell

    std::vector<RimPerforationInterval*>              perforations;
    std::vector<std::vector<RimPerforationInterval*>> perforationsOfConnection( connectionCells.size() );
    double                                            minFactor = std::numeric_limits<double>::max();
    double                                            maxFactor = 0.0;
    int                                               shutCount = 0;
    // The path in each connected cell: between its entry and exit points
    std::vector<std::pair<double, double>> pathInCell;
    for ( const auto& [entry, exit] : cellPointIndices )
    {
        pathInCell.emplace_back( mds[entry], mds[exit] );
    }

    // Cells that overlap, e.g. across faults over pinched out layers, share a part of the path. That part is left out
    // of both perforations, as a completion export assumes that the path is in one cell at a time.
    for ( size_t c = 1; c < pathInCell.size(); ++c )
    {
        auto& previous = pathInCell[c - 1];
        auto& current  = pathInCell[c];
        if ( current.first < previous.second )
        {
            std::swap( current.first, previous.second );
            previous.second = std::max( previous.second, previous.first );
            current.first   = std::min( current.first, current.second );
        }
    }

    const bool hasLgrs = mainGrid.gridCount() > 1;
    for ( size_t c = 0; c < connectionCells.size(); ++c )
    {
        const auto&  connection      = connectionCells[c].connection;
        const size_t cellIndex       = connectionCells[c].cellIndex;
        const auto   corners         = mainGrid.cellCornerVertices( cellIndex );
        const auto [entryMd, exitMd] = pathInCell[c];

        // Kept off the cell faces, so that a completion export does not see the neighbour cells
        const double margin  = perforationMargin( exitMd - entryMd );
        const double startMd = entryMd + margin;
        const double endMd   = exitMd - margin;
        if ( !( endMd > startMd ) )
        {
            warnings.push_back( QString( "Connection (%1, %2, %3) got no perforation; the path does not pass through the cell." )
                                    .arg( connection.i + 1 )
                                    .arg( connection.j + 1 )
                                    .arg( connection.k + 1 ) );
            continue;
        }

        const double diameter = connection.diameter > 0.0 ? connection.diameter : defaultDiameter;

        auto* perforation = new RimPerforationInterval();
        perforation->setStartAndEndMD( startMd, endMd );
        perforation->setDiameter( diameter );
        perforation->setSkinFactor( connection.skinFactor );
        perforation->setCheckState( connection.isOpen );
        if ( !connection.isOpen ) shutCount++;

        // The connection factor of the export is proportional to the perforated length in the cell. The cell intersections
        // are those of a completion export, unless LGRs hide the cell from it.
        auto transmissibility = [&]( const cvf::Vec3d& lengths )
        {
            return RicTransmissibilityCalculator::calculateTransmissibilityData( eclipseCase,
                                                                                 wellPath,
                                                                                 lengths,
                                                                                 connection.skinFactor,
                                                                                 diameter / 2.0,
                                                                                 cellIndex,
                                                                                 false )
                .connectionFactor();
        };

        double                                 unscaled = 0.0;
        std::vector<std::pair<double, double>> piecesInCell; // the perforated parts, through cell faces
        if ( !hasLgrs )
        {
            const auto [clippedPoints, clippedMds] = geometry->clippedPointSubset( startMd, endMd );
            std::set<size_t> otherCells;
            for ( const auto& intersection :
                  RigWellPathIntersectionTools::findCellIntersectionInfosAlongPath( caseData, name, clippedPoints, clippedMds ) )
            {
                if ( intersection.globCellIndex == cellIndex )
                {
                    unscaled += transmissibility( intersection.intersectionLengthsInCellCS );
                    if ( !piecesInCell.empty() && intersection.startMD - piecesInCell.back().second < 1.0e-6 )
                    {
                        piecesInCell.back().second = intersection.endMD;
                    }
                    else
                    {
                        piecesInCell.emplace_back( intersection.startMD, intersection.endMD );
                    }
                }
                else
                {
                    otherCells.insert( intersection.globCellIndex );
                }
            }
            // Parts of the path outside the cell, through its curved faces, are left unperforated
            if ( !otherCells.empty() && !piecesInCell.empty() )
            {
                otherCells.clear();
                unscaled = 0.0;
                for ( auto& [pieceStart, pieceEnd] : piecesInCell )
                {
                    const double pieceMargin = perforationMargin( pieceEnd - pieceStart );
                    pieceStart += pieceMargin;
                    pieceEnd -= pieceMargin;
                    const auto [piecePoints, pieceMds] = geometry->clippedPointSubset( pieceStart, pieceEnd );
                    for ( const auto& intersection :
                          RigWellPathIntersectionTools::findCellIntersectionInfosAlongPath( caseData, name, piecePoints, pieceMds ) )
                    {
                        if ( intersection.globCellIndex == cellIndex )
                        {
                            unscaled += transmissibility( intersection.intersectionLengthsInCellCS );
                        }
                    }
                }
            }
            else
            {
                piecesInCell.clear();
            }
            // A completion export does not see the cell along the path, e.g. where it overlaps other cells across a fault.
            // The connection is left out rather than moved to the other cells.
            if ( !otherCells.empty() )
            {
                warnings.push_back( QString( "Connection (%1, %2, %3) got no perforation; a completion export does not see the path in the "
                                             "cell, but in %4 other cell(s)." )
                                        .arg( connection.i + 1 )
                                        .arg( connection.j + 1 )
                                        .arg( connection.k + 1 )
                                        .arg( otherCells.size() ) );
                delete perforation;
                continue;
            }
        }
        if ( !( unscaled > 0.0 ) )
        {
            unscaled =
                transmissibility( RigWellPathIntersectionTools::calculateLengthInCell( corners,
                                                                                       geometry->interpolatedPointAlongWellPath( startMd ),
                                                                                       geometry->interpolatedPointAlongWellPath( endMd ) ) );
        }
        double factor = 1.0;
        if ( connection.connectionFactor > 0.0 && unscaled > 0.0 )
        {
            factor = connection.connectionFactor / unscaled;
        }
        else if ( connection.connectionFactor > 0.0 )
        {
            warnings.push_back( QString( "Connection (%1, %2, %3): no connection factor could be computed in the cell, the perforated "
                                         "length is not adjusted." )
                                    .arg( connection.i + 1 )
                                    .arg( connection.j + 1 )
                                    .arg( connection.k + 1 ) );
        }
        perforation->setPerforatedLengthFactor( factor );
        minFactor = std::min( minFactor, factor );
        maxFactor = std::max( maxFactor, factor );

        if ( piecesInCell.empty() )
        {
            perforations.push_back( perforation );
            perforationsOfConnection[c].push_back( perforation );
            continue;
        }
        for ( const auto& [pieceStart, pieceEnd] : piecesInCell )
        {
            auto* piece = new RimPerforationInterval();
            piece->setStartAndEndMD( pieceStart, pieceEnd );
            piece->setDiameter( diameter );
            piece->setSkinFactor( connection.skinFactor );
            piece->setCheckState( connection.isOpen );
            piece->setPerforatedLengthFactor( factor );
            perforations.push_back( piece );
            perforationsOfConnection[c].push_back( piece );
        }
        delete perforation;
    }

    wellPath->perforationIntervalCollection()->replacePerforations( perforations );

    // Multi-segment wells: the tubing and the valves
    int valveCount = 0;
    if ( !bore.empty() )
    {
        using PressureDrop = RimCompletionEnums::PressureDropType;
        completionSettings->setPressureDrop( well->pressureDropModel == 2   ? PressureDrop::HYDROSTATIC
                                             : well->pressureDropModel == 1 ? PressureDrop::HYDROSTATIC_FRICTION
                                                                            : PressureDrop::HYDROSTATIC_FRICTION_ACCELERATION );
        completionSettings->setLengthAndDepth( RimCompletionEnums::LengthAndDepthType::ABS );

        auto* segments = wellPath->segmentCollection();
        if ( bore.diameter > 0.0 ) segments->setLinerDiameter( bore.diameter );
        if ( bore.roughness > 0.0 ) segments->setRoughnessFactor( bore.roughness );
        if ( measuredDepthCalibrated )
        {
            segments->setReferenceMDType( RimCompletionEnums::ReferenceMDType::MANUAL_REFERENCE_MD );
            segments->setManualReferenceMD( bore.lengths.front() );
        }
        else
        {
            segments->setReferenceMDType( RimCompletionEnums::ReferenceMDType::AUTO_REFERENCE_MD );
        }

        // Connections that share a valve segment share its area
        std::map<int, int> connectionsOfSegment;
        for ( const auto& connectionCell : connectionCells )
        {
            connectionsOfSegment[connectionCell.connection.segment]++;
        }

        std::set<QString> unsupportedDevices;
        for ( size_t c = 0; c < connectionCells.size(); ++c )
        {
            const auto* segment = findSegment( *well, connectionCells[c].connection.segment );
            if ( !segment || segment->type == RifRestartWellConnectionReader::Segment::Type::REGULAR ) continue;
            if ( segment->type != RifRestartWellConnectionReader::Segment::Type::VALVE )
            {
                unsupportedDevices.insert( segment->type == RifRestartWellConnectionReader::Segment::Type::AUTO_ICD     ? "AICD"
                                           : segment->type == RifRestartWellConnectionReader::Segment::Type::SPIRAL_ICD ? "SICD"
                                                                                                                        : "other" );
                continue;
            }
            if ( perforationsOfConnection[c].empty() || !( segment->valveArea > 0.0 ) ) continue;

            // One valve over each perforation of the connection, sharing the area by length. A completion export shares the area
            // of a valve between the cells along it, e.g. the cells of a refined grid.
            double totalLength = 0.0;
            for ( auto* perforation : perforationsOfConnection[c] )
            {
                totalLength += perforation->endMD() - perforation->startMD();
            }
            if ( !( totalLength > 0.0 ) ) continue;

            const double area = segment->valveArea / std::max( 1, connectionsOfSegment[segment->number] );
            for ( auto* perforation : perforationsOfConnection[c] )
            {
                const double length = perforation->endMD() - perforation->startMD();
                auto* valveTemplate = valveTemplateFor( *segment, area * length / totalLength, caseData->unitsType(), simulationWellName );
                if ( !valveTemplate ) continue;

                auto* valve = new RimWellPathValve;
                valve->setValveTemplate( valveTemplate );
                perforation->addValve( valve );
                valve->setMeasuredDepthAndCount( perforation->startMD(), length, 1 );
                valve->applyValveLabelAndIcon();
            }
            valveCount++;
        }
        for ( const auto& device : unsupportedDevices )
        {
            warnings.push_back( QString( "Well '%1' has %2 devices, which are not reconstructed." ).arg( simulationWellName ).arg( device ) );
        }
    }

    wellPathCollection->updateConnectedEditors();
    wellPathCollection->scheduleRedrawAffectedViews();
    if ( auto* project = RimProject::current() ) project->scheduleCreateDisplayModelAndRedrawAllViews();

    RiaLogging::info( QString( "Created the approximate well path '%1' for well '%2' at %3: %4 perforations (%5 shut), perforated "
                               "length factors %6 - %7%8." )
                          .arg( name )
                          .arg( simulationWellName )
                          .arg( date.toString( Qt::ISODate ) )
                          .arg( perforations.size() )
                          .arg( shutCount )
                          .arg( perforations.empty() ? 0.0 : minFactor, 0, 'g', 3 )
                          .arg( perforations.empty() ? 0.0 : maxFactor, 0, 'g', 3 )
                          .arg( bore.empty() ? QString()
                                             : QString( ", multi-segment data from %1 segments (%2), %3 valves" )
                                                   .arg( well->segments.size() )
                                                   .arg( intervalsUsed ? "with the connection intervals along the tubing" : "the tubing only" )
                                                   .arg( valveCount ) )
                          .toStdString() );
    for ( const auto& warning : warnings )
    {
        RiaLogging::warning( warning.toStdString() );
    }

    return wellPath;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int RicApproximateWellPathTools::viewTimeStepWithConnections( const RimSimWellInView& simWell )
{
    auto* view = simWell.firstAncestorOrThisOfType<RimEclipseView>();
    if ( !view || !simWell.simWellData() ) return -1;

    const int timeStep = view->currentTimeStep();
    if ( timeStep < 0 || !simWell.simWellData()->hasWellResult( timeStep ) ) return -1;

    const auto* frame = simWell.simWellData()->wellResultFrame( timeStep );
    if ( !frame ) return -1;
    const auto points = frame->allResultPoints();
    return std::any_of( points.begin(), points.end(), []( const RigWellResultPoint& p ) { return p.isCell(); } ) ? timeStep : -1;
}
