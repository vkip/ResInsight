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
#include "Well/RigWellPathIntersectionTools.h"
#include "Well/RigWellResultFrame.h"
#include "Well/RigWellResultPoint.h"

#include "RimEclipseCase.h"
#include "RimPerforationCollection.h"
#include "RimPerforationInterval.h"
#include "RimPointBasedWellPath.h"
#include "RimProject.h"
#include "RimWellPathCollection.h"
#include "RimWellPathCompletionSettings.h"

#include <QFileInfo>

#include <algorithm>
#include <limits>
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
/// The centre of the face of the cell that the direction points at most
//--------------------------------------------------------------------------------------------------
cvf::Vec3d exitFaceCenter( const RigCell& cell, const cvf::Vec3d& direction )
{
    const cvf::Vec3d center = cell.center();
    cvf::Vec3d       unitDirection( direction );
    if ( !unitDirection.normalize() ) unitDirection = cvf::Vec3d( 0, 0, -1 );

    cvf::Vec3d best      = center;
    double     bestAlign = -std::numeric_limits<double>::max();
    for ( FaceType face : { FaceType::NEG_I, FaceType::POS_I, FaceType::NEG_J, FaceType::POS_J, FaceType::NEG_K, FaceType::POS_K } )
    {
        cvf::Vec3d toFace = cell.faceCenter( face ) - center;
        if ( !toFace.normalize() ) continue;
        const double align = toFace * unitDirection;
        if ( align > bestAlign )
        {
            bestAlign = align;
            best      = cell.faceCenter( face );
        }
    }
    return best;
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
std::expected<RimWellPath*, QString>
    RicApproximateWellPathTools::createApproximateWellPath( RimEclipseCase* eclipseCase, const QString& simulationWellName, int timeStep )
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
    if ( well->isMultiSegment )
    {
        warnings.push_back( QString( "Well '%1' is a multi-segment well; only its connections are used." ).arg( simulationWellName ) );
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
        connectionCells.push_back( { connection, cellIndex } );
    }
    if ( connectionCells.empty() )
    {
        return std::unexpected(
            QString( "Well '%1' has no connections in the grid at %2." ).arg( simulationWellName ).arg( date.toString( Qt::ISODate ) ) );
    }

    // The path: from sea level above the well head to the top of its column, then through the cell centres
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

    std::vector<size_t> centerPointIndices;
    for ( const auto& connectionCell : connectionCells )
    {
        centerPointIndices.push_back( points.size() );
        points.push_back( mainGrid.cell( connectionCell.cellIndex ).center() );
    }
    {
        const RigCell& last = mainGrid.cell( connectionCells.back().cellIndex );
        points.push_back( exitFaceCenter( last, points.back() - points[points.size() - 2] ) );
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
    wellPath->setUnitSystem( caseData->unitsType() );
    wellPath->setTrajectoryPoints( points );

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
    const double defaultDiameter = caseData->unitsType() == RiaDefines::EclipseUnitSystem::UNITS_FIELD ? 0.709 : 0.216;

    std::vector<RimPerforationInterval*> perforations;
    double                               minFactor = std::numeric_limits<double>::max();
    double                               maxFactor = 0.0;
    int                                  shutCount = 0;
    // The path in each connected cell: from the last face crossing before the centre to the first one after it
    std::vector<std::pair<double, double>> pathInCell;
    for ( size_t c = 0; c < connectionCells.size(); ++c )
    {
        const size_t centre  = centerPointIndices[c];
        const auto   corners = mainGrid.cellCornerVertices( connectionCells[c].cellIndex );

        double entryMd = mds[centre - 1];
        if ( auto crossings = faceCrossings( points[centre - 1], points[centre], corners ); !crossings.empty() )
        {
            entryMd = mds[centre - 1] + *std::max_element( crossings.begin(), crossings.end() ) * ( mds[centre] - mds[centre - 1] );
        }
        double exitMd = mds[centre + 1];
        if ( auto crossings = faceCrossings( points[centre], points[centre + 1], corners ); !crossings.empty() )
        {
            exitMd = mds[centre] + *std::min_element( crossings.begin(), crossings.end() ) * ( mds[centre + 1] - mds[centre] );
        }
        pathInCell.emplace_back( entryMd, exitMd );
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
        const double margin  = 1.0e-3 * ( exitMd - entryMd );
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

        double unscaled = 0.0;
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
                }
                else
                {
                    otherCells.insert( intersection.globCellIndex );
                }
            }
            if ( !otherCells.empty() )
            {
                warnings.push_back( QString( "The perforation of connection (%1, %2, %3) also crosses %4 other cell(s), which a completion "
                                             "export includes." )
                                        .arg( connection.i + 1 )
                                        .arg( connection.j + 1 )
                                        .arg( connection.k + 1 )
                                        .arg( otherCells.size() ) );
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

        perforations.push_back( perforation );
    }

    wellPath->perforationIntervalCollection()->replacePerforations( perforations );

    wellPathCollection->updateConnectedEditors();
    wellPathCollection->scheduleRedrawAffectedViews();
    if ( auto* project = RimProject::current() ) project->scheduleCreateDisplayModelAndRedrawAllViews();

    RiaLogging::info( QString( "Created the approximate well path '%1' for well '%2' at %3: %4 perforations (%5 shut), perforated "
                               "length factors %6 - %7." )
                          .arg( name )
                          .arg( simulationWellName )
                          .arg( date.toString( Qt::ISODate ) )
                          .arg( perforations.size() )
                          .arg( shutCount )
                          .arg( perforations.empty() ? 0.0 : minFactor, 0, 'g', 3 )
                          .arg( perforations.empty() ? 0.0 : maxFactor, 0, 'g', 3 )
                          .toStdString() );
    for ( const auto& warning : warnings )
    {
        RiaLogging::warning( warning.toStdString() );
    }

    return wellPath;
}
