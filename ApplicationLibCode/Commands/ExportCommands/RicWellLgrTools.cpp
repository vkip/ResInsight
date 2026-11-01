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

#include "RicWellLgrTools.h"

#include "RiaLogging.h"

#include "RicCreateTemporaryLgrFeature.h"
#include "RicDeleteTemporaryLgrsFeature.h"
#include "RicExportLgrFeature.h"

#include "RigCell.h"
#include "RigEclipseCaseData.h"
#include "RigLocalGrid.h"
#include "RigMainGrid.h"
#include "RigNNCData.h"
#include "RigNonUniformRefinement.h"
#include "RigReservoirGridTools.h"
#include "RigUniformRefinement.h"
#include "Well/RigSimWellData.h"
#include "Well/RigWellResultFrame.h"
#include "Well/RigWellResultPoint.h"

#include "RimEclipseCase.h"

#include <QRegularExpression>
#include <QTextStream>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace caf
{
template <>
void RicWellLgrTools::LateralRefinementEnum::setUp()
{
    addItem( RicWellLgrTools::LateralRefinement::UNIFORM, "UNIFORM", "Uniform" );
    addItem( RicWellLgrTools::LateralRefinement::LOGARITHMIC, "LOGARITHMIC", "Logarithmic" );
    setDefault( RicWellLgrTools::LateralRefinement::LOGARITHMIC );
}
} // namespace caf

namespace
{
//--------------------------------------------------------------------------------------------------
/// Main grid cell hosting a cell of any grid
//--------------------------------------------------------------------------------------------------
size_t hostMainGridCellIndex( const RigMainGrid& mainGrid, size_t reservoirCellIndex )
{
    size_t cellIndex = reservoirCellIndex;
    while ( cellIndex != cvf::UNDEFINED_SIZE_T && !mainGrid.cell( cellIndex ).hostGrid()->isMainGrid() )
    {
        cellIndex = mainGrid.cell( cellIndex ).parentCellIndex();
    }
    return cellIndex;
}

//--------------------------------------------------------------------------------------------------
/// Average width of the main grid cells in a slice of the box, one value per parent cell along the direction
//--------------------------------------------------------------------------------------------------
std::vector<double>
    parentCellWidths( const RigMainGrid& mainGrid, const caf::VecIjk0& boxMin, const caf::VecIjk0& boxMax, RigRefinement::Dimension dim )
{
    using FaceType               = cvf::StructGridInterface::FaceType;
    const FaceType negFaces[3]   = { FaceType::NEG_I, FaceType::NEG_J, FaceType::NEG_K };
    const FaceType posFaces[3]   = { FaceType::POS_I, FaceType::POS_J, FaceType::POS_K };
    const size_t   count         = boxMax[dim] - boxMin[dim] + 1;
    double         defaultLength = mainGrid.characteristicIJCellSize();
    if ( !( defaultLength > 0.0 ) ) defaultLength = 1.0;

    std::vector<double> widths( count, 0.0 );
    for ( size_t p = 0; p < count; ++p )
    {
        double sum = 0.0;
        size_t n   = 0;
        for ( size_t k = boxMin.k(); k <= boxMax.k(); ++k )
        {
            for ( size_t j = boxMin.j(); j <= boxMax.j(); ++j )
            {
                for ( size_t i = boxMin.i(); i <= boxMax.i(); ++i )
                {
                    caf::VecIjk0 ijk( i, j, k );
                    if ( ijk[dim] != boxMin[dim] + p ) continue;

                    const RigCell& cell = mainGrid.cell( mainGrid.cellIndexFromIJK( i, j, k ) );
                    if ( cell.isInvalid() ) continue;

                    const double length = ( cell.faceCenter( posFaces[dim] ) - cell.faceCenter( negFaces[dim] ) ).length();
                    if ( length > 0.0 )
                    {
                        sum += length;
                        n++;
                    }
                }
            }
        }
        widths[p] = n > 0 ? sum / n : defaultLength;
    }
    return widths;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<double> equalFractions( size_t count )
{
    return RigRefinement::generateEqualFractions( std::max( count, size_t( 1 ) ) );
}

//--------------------------------------------------------------------------------------------------
/// Whether all parent cells of a direction are split into the same number of equally sized cells
//--------------------------------------------------------------------------------------------------
bool isUniform( const RigRefinement& refinement, RigRefinement::Dimension dim )
{
    const size_t count = refinement.sectorSize( dim );
    if ( count == 0 ) return true;

    const auto reference = equalFractions( refinement.subcellCount( dim, 0 ) );
    for ( size_t p = 0; p < count; ++p )
    {
        const auto& fractions = refinement.cumulativeFractions( dim, p );
        if ( fractions.size() != reference.size() ) return false;
        for ( size_t s = 0; s < fractions.size(); ++s )
        {
            if ( std::abs( fractions[s] - reference[s] ) > 1.0e-9 ) return false;
        }
    }
    return true;
}

//--------------------------------------------------------------------------------------------------
/// Values of a keyword, eight per line
//--------------------------------------------------------------------------------------------------
template <typename T>
void writeValues( QTextStream& stream, const QString& keyword, const std::vector<T>& values, int precision )
{
    stream << keyword << "\n";
    for ( size_t n = 0; n < values.size(); ++n )
    {
        stream << ( n % 8 == 0 ? "  " : " " );
        if constexpr ( std::is_floating_point_v<T> )
            stream << QString::number( values[n], 'g', precision );
        else
            stream << values[n];
        if ( n % 8 == 7 || n + 1 == values.size() ) stream << "\n";
    }
    stream << "/\n";
}
} // namespace

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RicWellLgrTools::CellSet RicWellLgrTools::simulationWellCells( const RigEclipseCaseData& caseData, const QString& wellName, int timeStep )
{
    CellSet cells;

    const RigSimWellData* wellData = caseData.findSimWellData( wellName );
    const RigMainGrid*    mainGrid = caseData.mainGrid();
    if ( !wellData || !mainGrid ) return cells;

    auto addFrame = [&]( const RigWellResultFrame& frame )
    {
        for ( const auto& point : frame.allResultPoints() )
        {
            if ( !point.isCell() ) continue;

            const RigGridBase* grid = mainGrid->gridByIndex( point.gridIndex() );
            if ( !grid ) continue;

            const size_t hostIndex = hostMainGridCellIndex( *mainGrid, grid->reservoirCellIndex( point.cellIndex() ) );
            if ( auto ijk = mainGrid->ijkFromCellIndex( hostIndex ) ) cells.insert( *ijk );
        }
    };

    if ( timeStep >= 0 )
    {
        if ( const auto* frame = wellData->wellResultFrame( static_cast<size_t>( timeStep ) ) ) addFrame( *frame );
    }
    else
    {
        for ( const auto& frame : wellData->m_wellCellsTimeSteps )
        {
            addFrame( frame );
        }
    }

    return cells;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> RicWellLgrTools::gradedCumulativeFractions( const std::vector<double>& parentCellWidths,
                                                                             const std::set<size_t>&    wellParentCells,
                                                                             double                     minCellSize,
                                                                             double                     growthFactor )
{
    const size_t                     parentCount = parentCellWidths.size();
    std::vector<std::vector<double>> fractions( parentCount, std::vector<double>{ 1.0 } );
    if ( parentCount == 0 || wellParentCells.empty() || !( minCellSize > 0.0 ) ) return fractions;

    const double growthRate = std::max( growthFactor - 1.0, 1.0e-3 );

    std::vector<double> parentStart( parentCount + 1, 0.0 );
    for ( size_t p = 0; p < parentCount; ++p )
    {
        parentStart[p + 1] = parentStart[p] + parentCellWidths[p];
    }

    std::vector<double> wellPositions;
    for ( size_t p : wellParentCells )
    {
        if ( p < parentCount ) wellPositions.push_back( 0.5 * ( parentStart[p] + parentStart[p + 1] ) );
    }
    if ( wellPositions.empty() ) return fractions;

    auto cellSize = [&]( double x )
    {
        double distance = std::numeric_limits<double>::max();
        for ( double position : wellPositions )
        {
            distance = std::min( distance, std::abs( x - position ) );
        }
        return std::max( minCellSize, growthRate * distance );
    };

    // The cell boundaries divide the integral of 1/h into equal parts within each parent cell
    const size_t samplesPerCell    = 2000;
    const size_t maxCellsPerParent = 500;
    for ( size_t p = 0; p < parentCount; ++p )
    {
        const double width = parentCellWidths[p];
        if ( !( width > 0.0 ) ) continue;

        std::vector<double> cumulative( samplesPerCell + 1, 0.0 );
        const double        dx = width / samplesPerCell;
        for ( size_t s = 0; s < samplesPerCell; ++s )
        {
            cumulative[s + 1] = cumulative[s] + dx / cellSize( parentStart[p] + ( s + 0.5 ) * dx );
        }

        const double total = cumulative.back();
        const size_t count = std::clamp( static_cast<size_t>( std::lround( total ) ), size_t( 1 ), maxCellsPerParent );

        std::vector<double> parentFractions;
        size_t              s = 0;
        for ( size_t c = 1; c < count; ++c )
        {
            const double target = total * c / count;
            while ( s < samplesPerCell && cumulative[s + 1] < target )
            {
                s++;
            }
            const double segment = cumulative[s + 1] - cumulative[s];
            const double t       = segment > 0.0 ? ( target - cumulative[s] ) / segment : 0.0;
            parentFractions.push_back( ( s + t ) / samplesPerCell );
        }
        parentFractions.push_back( 1.0 );
        fractions[p] = parentFractions;
    }

    return fractions;
}

//--------------------------------------------------------------------------------------------------
/// Eclipse LGR names hold at most eight characters
//--------------------------------------------------------------------------------------------------
QString RicWellLgrTools::defaultLgrName( const QString& wellName )
{
    QString name = "L" + wellName.toUpper();
    name.remove( QRegularExpression( "[^A-Z0-9_\\-]" ) );
    return name.left( 8 );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<LgrInfo, QString>
    RicWellLgrTools::buildLgrForSimulationWell( RimEclipseCase& eclipseCase, const QString& wellName, const Settings& settings )
{
    const RigEclipseCaseData* caseData = eclipseCase.eclipseCaseData();
    if ( !caseData || !caseData->mainGrid() ) return std::unexpected( "The case has no grid." );
    const RigMainGrid& mainGrid = *caseData->mainGrid();

    if ( settings.refinementK < 1 || settings.refinementI < 1 || settings.refinementJ < 1 )
    {
        return std::unexpected( "Refinement must be at least 1 in all directions." );
    }
    if ( settings.lateralRefinement == LateralRefinement::LOGARITHMIC && ( !( settings.minCellSize > 0.0 ) || !( settings.growthFactor > 1.0 ) ) )
    {
        return std::unexpected( "Logarithmic refinement needs a positive minimum cell size and a growth factor above 1." );
    }

    const auto wellCells = simulationWellCells( *caseData, wellName, settings.timeStep );
    if ( wellCells.empty() ) return std::unexpected( QString( "Well '%1' has no connections in the main grid." ).arg( wellName ) );

    caf::VecIjk0 wellMin = *wellCells.begin();
    caf::VecIjk0 wellMax = *wellCells.begin();
    for ( const auto& c : wellCells )
    {
        for ( size_t d = 0; d < 3; ++d )
        {
            wellMin[d] = std::min( wellMin[d], c[d] );
            wellMax[d] = std::max( wellMax[d], c[d] );
        }
    }

    const size_t gridCount[3] = { mainGrid.cellCountI(), mainGrid.cellCountJ(), mainGrid.cellCountK() };
    const int    padding[3]   = { settings.paddingIJ, settings.paddingIJ, settings.paddingK };
    caf::VecIjk0 boxMin       = wellMin;
    caf::VecIjk0 boxMax       = wellMax;
    for ( size_t d = 0; d < 3; ++d )
    {
        const size_t pad = static_cast<size_t>( std::max( padding[d], 0 ) );
        boxMin[d]        = wellMin[d] > pad ? wellMin[d] - pad : 0;
        boxMax[d]        = std::min( wellMax[d] + pad, gridCount[d] - 1 );
    }

    // The LGR must not overlap other LGRs
    for ( size_t k = boxMin.k(); k <= boxMax.k(); ++k )
    {
        for ( size_t j = boxMin.j(); j <= boxMax.j(); ++j )
        {
            for ( size_t i = boxMin.i(); i <= boxMax.i(); ++i )
            {
                const RigGridBase* subGrid = mainGrid.cell( mainGrid.cellIndexFromIJK( i, j, k ) ).subGrid();
                if ( subGrid && subGrid->associatedWellPathName() != wellName.toStdString() )
                {
                    return std::unexpected( QString( "The LGR around well '%1' would overlap the LGR '%2' at cell (%3, %4, %5)." )
                                                .arg( wellName )
                                                .arg( QString::fromStdString( subGrid->gridName() ) )
                                                .arg( i + 1 )
                                                .arg( j + 1 )
                                                .arg( k + 1 ) );
                }
            }
        }
    }

    // The LGR must not hold connections of other wells, which would refer to refined host cells
    for ( const auto& wellData : caseData->wellResults() )
    {
        if ( wellData->m_wellName == wellName ) continue;
        for ( const auto& c : simulationWellCells( *caseData, wellData->m_wellName ) )
        {
            if ( c.i() >= boxMin.i() && c.i() <= boxMax.i() && c.j() >= boxMin.j() && c.j() <= boxMax.j() && c.k() >= boxMin.k() &&
                 c.k() <= boxMax.k() )
            {
                return std::unexpected( QString( "The LGR around well '%1' would hold connections of well '%2' at cell (%3, %4, %5). "
                                                 "Use a smaller padding." )
                                            .arg( wellName )
                                            .arg( wellData->m_wellName )
                                            .arg( c.i() + 1 )
                                            .arg( c.j() + 1 )
                                            .arg( c.k() + 1 ) );
            }
        }
    }

    // OPM Flow's default LGR backend cannot refine cells with more than six faces, e.g. cells at faults with throw
    if ( auto nncData = const_cast<RigMainGrid&>( mainGrid ).nncData() )
    {
        size_t faultedHostCells = 0;
        for ( const auto& connection : nncData->availableConnections() )
        {
            for ( size_t cellIndex : { connection.c1GlobIdx(), connection.c2GlobIdx() } )
            {
                if ( cellIndex >= mainGrid.cellCount() ) continue;
                if ( auto ijk = mainGrid.ijkFromCellIndex( cellIndex ); ijk && ijk->i() >= boxMin.i() && ijk->i() <= boxMax.i() &&
                                                                        ijk->j() >= boxMin.j() && ijk->j() <= boxMax.j() &&
                                                                        ijk->k() >= boxMin.k() && ijk->k() <= boxMax.k() )
                {
                    faultedHostCells++;
                }
            }
        }
        if ( faultedHostCells > 0 )
        {
            RiaLogging::warning(
                QString( "The LGR around well '%1' has %2 non-neighbour connections to its host cells, e.g. across faults. "
                         "OPM Flow can currently only refine cells with six faces." )
                    .arg( wellName )
                    .arg( faultedHostCells )
                    .toStdString() );
        }
    }

    const QString     lgrName = settings.lgrName.isEmpty() ? defaultLgrName( wellName ) : settings.lgrName.left( 8 );
    const cvf::Vec3st boxSize( boxMax.i() - boxMin.i() + 1, boxMax.j() - boxMin.j() + 1, boxMax.k() - boxMin.k() + 1 );

    if ( settings.lateralRefinement == LateralRefinement::UNIFORM )
    {
        cvf::Vec3st refinement( settings.refinementI, settings.refinementJ, settings.refinementK );
        return LgrInfo( 0, lgrName, wellName, refinement, boxMin, boxMax );
    }

    auto graded = std::make_shared<RigNonUniformRefinement>( boxSize );

    // A direction the well runs along over several cells is not graded, the others are graded around the well cells
    std::set<size_t> wellParents[2];
    for ( const auto& c : wellCells )
    {
        wellParents[0].insert( c.i() - boxMin.i() );
        wellParents[1].insert( c.j() - boxMin.j() );
    }
    const size_t alongWellThreshold = 3;
    int          alongWellDim       = -1;
    if ( std::max( wellParents[0].size(), wellParents[1].size() ) >= alongWellThreshold && wellParents[0].size() != wellParents[1].size() )
    {
        alongWellDim = wellParents[0].size() > wellParents[1].size() ? 0 : 1;
    }

    for ( int d = 0; d < 2; ++d )
    {
        const auto dim = static_cast<RigRefinement::Dimension>( d );
        if ( d == alongWellDim ) continue;

        const auto widths    = parentCellWidths( mainGrid, boxMin, boxMax, dim );
        const auto fractions = gradedCumulativeFractions( widths, wellParents[d], settings.minCellSize, settings.growthFactor );
        for ( size_t p = 0; p < fractions.size(); ++p )
        {
            graded->setCumulativeFractions( dim, p, fractions[p] );
        }
    }
    for ( size_t k = 0; k < boxSize.z(); ++k )
    {
        graded->setCumulativeFractions( RigRefinement::DimK, k, equalFractions( settings.refinementK ) );
    }

    LgrInfo info( 0, lgrName, wellName, cvf::Vec3st( 1, 1, settings.refinementK ), boxMin, boxMax );
    info.gradedRefinement = graded;
    return info;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<LgrInfo, QString>
    RicWellLgrTools::createTemporaryLgrForSimulationWell( RimEclipseCase* eclipseCase, const QString& wellName, const Settings& settings )
{
    if ( !eclipseCase || !eclipseCase->eclipseCaseData() ) return std::unexpected( "No case is available." );

    auto newLgr = buildLgrForSimulationWell( *eclipseCase, wellName, settings );
    if ( !newLgr ) return newLgr;

    // Deleting temporary LGRs reloads the grid, so the LGRs of other wells are recreated afterwards
    std::vector<LgrInfo> lgrsToKeep;
    int                  maxRealGridId = 0;
    {
        const RigMainGrid& mainGrid = *eclipseCase->eclipseCaseData()->mainGrid();
        for ( size_t g = 1; g < mainGrid.gridCount(); ++g )
        {
            const RigGridBase* grid = mainGrid.gridByIndex( g );
            if ( !grid->isTempGrid() )
            {
                maxRealGridId = std::max( maxRealGridId, grid->gridId() );
                if ( QString::fromStdString( grid->gridName() ).compare( newLgr->name, Qt::CaseInsensitive ) == 0 )
                {
                    return std::unexpected( QString( "The case already has an LGR named '%1'." ).arg( newLgr->name ) );
                }
                continue;
            }

            const auto* localGrid = dynamic_cast<const RigLocalGrid*>( grid );
            if ( !localGrid || localGrid->associatedWellPathName() == wellName.toStdString() ||
                 QString::fromStdString( grid->gridName() ).compare( newLgr->name, Qt::CaseInsensitive ) == 0 )
            {
                continue;
            }
            lgrsToKeep.push_back( lgrInfoFromTemporaryLgr( mainGrid, *localGrid ) );
        }
    }

    RicDeleteTemporaryLgrsFeature::deleteAllTemporaryLgrs( eclipseCase );

    lgrsToKeep.push_back( *newLgr );
    int nextId = maxRealGridId + 1;
    for ( auto& lgr : lgrsToKeep )
    {
        lgr.id = nextId++;
        RicCreateTemporaryLgrFeature::createLgr( lgr, eclipseCase->eclipseCaseData() );
    }

    RigReservoirGridTools::refreshEclipseCaseDataAndViews( eclipseCase );

    const auto cellCounts = lgrsToKeep.back().lgrCellCounts();
    RiaLogging::info( QString( "Created LGR '%1' around well '%2': host cells I %3-%4, J %5-%6, K %7-%8, %9 x %10 x %11 cells" )
                          .arg( lgrsToKeep.back().name )
                          .arg( wellName )
                          .arg( lgrsToKeep.back().mainGridStartCell.i() + 1 )
                          .arg( lgrsToKeep.back().mainGridEndCell.i() + 1 )
                          .arg( lgrsToKeep.back().mainGridStartCell.j() + 1 )
                          .arg( lgrsToKeep.back().mainGridEndCell.j() + 1 )
                          .arg( lgrsToKeep.back().mainGridStartCell.k() + 1 )
                          .arg( lgrsToKeep.back().mainGridEndCell.k() + 1 )
                          .arg( cellCounts.i() )
                          .arg( cellCounts.j() )
                          .arg( cellCounts.k() )
                          .toStdString() );

    return lgrsToKeep.back();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
const RigLocalGrid* RicWellLgrTools::temporaryLgrForWell( const RigMainGrid& mainGrid, const QString& wellName )
{
    for ( size_t g = 1; g < mainGrid.gridCount(); ++g )
    {
        const auto* localGrid = dynamic_cast<const RigLocalGrid*>( mainGrid.gridByIndex( g ) );
        if ( localGrid && localGrid->isTempGrid() && localGrid->associatedWellPathName() == wellName.toStdString() )
        {
            return localGrid;
        }
    }
    return nullptr;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
LgrInfo RicWellLgrTools::lgrInfoFromTemporaryLgr( const RigMainGrid& mainGrid, const RigLocalGrid& localGrid )
{
    const auto lgrInfos = RicExportLgrFeature::createLgrInfoListForTemporaryLgrs( &mainGrid );
    for ( const auto& [wellName, infos] : lgrInfos )
    {
        for ( const auto& info : infos )
        {
            if ( info.name.toStdString() == localGrid.gridName() ) return info;
        }
    }
    return LgrInfo( localGrid.gridId(),
                    QString::fromStdString( localGrid.gridName() ),
                    QString::fromStdString( localGrid.associatedWellPathName() ),
                    cvf::Vec3st( 1, 1, 1 ),
                    caf::VecIjk0( 0, 0, 0 ),
                    caf::VecIjk0( 0, 0, 0 ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::unique_ptr<RigRefinement> RicWellLgrTools::lgrRefinement( const LgrInfo& lgrInfo )
{
    if ( lgrInfo.gradedRefinement ) return lgrInfo.gradedRefinement->clone();

    const cvf::Vec3st boxSize( lgrInfo.mainGridEndCell.i() - lgrInfo.mainGridStartCell.i() + 1,
                               lgrInfo.mainGridEndCell.j() - lgrInfo.mainGridStartCell.j() + 1,
                               lgrInfo.mainGridEndCell.k() - lgrInfo.mainGridStartCell.k() + 1 );
    return std::make_unique<RigUniformRefinement>( lgrInfo.refinement, boxSize );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RicWellLgrTools::gradingKeywords( const RigRefinement& refinement )
{
    QString     text;
    QTextStream stream( &text );

    const char* dimNames[3] = { "X", "Y", "Z" };
    for ( size_t d = 0; d < 3; ++d )
    {
        const auto dim = static_cast<RigRefinement::Dimension>( d );
        if ( isUniform( refinement, dim ) ) continue;

        std::vector<int>    counts;
        std::vector<double> widths;
        for ( size_t p = 0; p < refinement.sectorSize( dim ); ++p )
        {
            const auto& fractions = refinement.cumulativeFractions( dim, p );
            counts.push_back( static_cast<int>( fractions.size() ) );
            double previous = 0.0;
            for ( double fraction : fractions )
            {
                widths.push_back( fraction - previous );
                previous = fraction;
            }
        }

        stream << "\n";
        writeValues( stream, QString( "N%1FIN" ).arg( dimNames[d] ), counts, 0 );
        stream << "\n";
        writeValues( stream, QString( "H%1FIN" ).arg( dimNames[d] ), widths, 5 );
    }

    return text;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RicWellLgrTools::carfinBlock( const LgrInfo& lgrInfo, int maxWellCount )
{
    const caf::VecIjk1 start  = lgrInfo.mainGridStartCell.toOneBased();
    const caf::VecIjk1 end    = lgrInfo.mainGridEndCell.toOneBased();
    const cvf::Vec3st  counts = lgrInfo.lgrCellCounts();

    QString text = QString( "-- LGR %1 around well %2\n"
                            "CARFIN\n"
                            "-- NAME I1 I2 J1 J2 K1 K2 NX NY NZ NWMAX\n"
                            "  '%1' %3 %4 %5 %6 %7 %8 %9 %10 %11 %12 /\n" )
                       .arg( lgrInfo.name )
                       .arg( lgrInfo.associatedWellPathName )
                       .arg( start.i() )
                       .arg( end.i() )
                       .arg( start.j() )
                       .arg( end.j() )
                       .arg( start.k() )
                       .arg( end.k() )
                       .arg( counts.i() )
                       .arg( counts.j() )
                       .arg( counts.k() )
                       .arg( maxWellCount );

    if ( lgrInfo.gradedRefinement ) text += gradingKeywords( *lgrInfo.gradedRefinement );

    text += "\nENDFIN\n";
    return text;
}
