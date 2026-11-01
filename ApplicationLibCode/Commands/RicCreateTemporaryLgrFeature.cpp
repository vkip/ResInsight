/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2017 Statoil ASA
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

#include "RicCreateTemporaryLgrFeature.h"

#include "RiaCellDividingTools.h"
#include "RiaGuiApplication.h"
#include "RiaLogging.h"
#include "RiuMessageDialog.h"

#include "ExportCommands/RicExportLgrFeature.h"
#include "ExportCommands/RicExportLgrUi.h"
#include "RicDeleteTemporaryLgrsFeature.h"

#include "RigActiveCellInfo.h"
#include "RigCaseCellResultsData.h"
#include "RigCell.h"
#include "RigEclipseCaseData.h"
#include "RigLocalGrid.h"
#include "RigMainGrid.h"
#include "RigReservoirGridTools.h"

#include "RimEclipseCase.h"
#include "RimEclipseView.h"
#include "RimGridCollection.h"
#include "RimMainPlotCollection.h"
#include "RimWellLogPlotCollection.h"
#include "RimWellPath.h"
#include "RimWellPathCompletions.h"

#include <QAction>
#include <QDir>

#include <cafSelectionManagerTools.h>
#include <cafUtils.h>
#include <cafVecIjk.h>

#include <array>
#include <set>

CAF_CMD_SOURCE_INIT( RicCreateTemporaryLgrFeature, "RicCreateTemporaryLgrFeature" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateTemporaryLgrFeature::createLgrsForWellPaths( std::vector<RimWellPath*>                          wellPaths,
                                                           RimEclipseCase*                                    eclipseCase,
                                                           size_t                                             timeStep,
                                                           const cvf::Vec3st&                                 refinement,
                                                           Lgr::SplitType                                     splitType,
                                                           const std::set<RigCompletionData::CompletionType>& completionTypes,
                                                           QStringList*                                       wellsIntersectingOtherLgrs )
{
    auto lgrs =
        RicExportLgrFeature::buildLgrsForWellPaths( wellPaths, eclipseCase, timeStep, refinement, splitType, completionTypes, wellsIntersectingOtherLgrs );

    for ( const auto& lgr : lgrs )
    {
        createLgr( lgr, eclipseCase->eclipseCaseData() );
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RicCreateTemporaryLgrFeature::isCommandEnabled() const
{
    std::vector<RimWellPathCompletions*> completions = caf::selectedObjectsByTypeStrict<RimWellPathCompletions*>();
    std::vector<RimWellPath*>            wellPaths   = caf::selectedObjectsByTypeStrict<RimWellPath*>();

    return !completions.empty() || !wellPaths.empty();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateTemporaryLgrFeature::onActionTriggered( bool isChecked )
{
    std::vector<RimWellPath*> wellPaths = RicExportLgrFeature::selectedWellPaths();
    if ( wellPaths.empty() ) return;

    QString dialogTitle = "Create Temporary LGR";

    RimEclipseCase* defaultEclipseCase = nullptr;
    int             defaultTimeStep    = 0;
    auto            activeView         = dynamic_cast<RimEclipseView*>( RiaApplication::instance()->activeGridView() );
    if ( activeView )
    {
        defaultEclipseCase = activeView->eclipseCase();
        defaultTimeStep    = activeView->currentTimeStep();
    }

    auto dialogData = RicExportLgrFeature::openDialog( dialogTitle, defaultEclipseCase, defaultTimeStep, true );
    if ( dialogData )
    {
        auto       eclipseCase     = dialogData->caseToApply();
        auto       refinement      = dialogData->refinement();
        size_t     timeStep        = dialogData->timeStep();
        auto       splitType       = dialogData->splitType();
        const auto completionTypes = dialogData->completionTypes();

        RicDeleteTemporaryLgrsFeature::deleteAllTemporaryLgrs( eclipseCase );

        QStringList wellsIntersectingOtherLgrs;

        createLgrsForWellPaths( wellPaths, eclipseCase, timeStep, refinement, splitType, completionTypes, &wellsIntersectingOtherLgrs );

        RigReservoirGridTools::refreshEclipseCaseDataAndViews( eclipseCase );

        if ( !wellsIntersectingOtherLgrs.empty() )
        {
            RiuMessageDialog::showError( nullptr,
                                         "LGR cells intersected",
                                         "At least one completion intersects with an LGR. No LGR(s) for those cells "
                                         "are produced" );
        }
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateTemporaryLgrFeature::setupActionLook( QAction* actionToSetup )
{
    actionToSetup->setText( "Create Temporary LGR" );
    actionToSetup->setIcon( QIcon( ":/TempLGR16x16.png" ) );
}

//--------------------------------------------------------------------------------------------------
/// Todo: Guarding, caching LGR corner nodes calculations
//--------------------------------------------------------------------------------------------------
RigGridBase* RicCreateTemporaryLgrFeature::createLgrForGrid( const LgrInfo& lgrInfo, RigMainGrid* mainGrid )
{
    int    lgrId          = lgrInfo.id;
    size_t totalCellCount = mainGrid->totalCellCount();
    size_t lgrCellCount   = lgrInfo.lgrCellCount();

    // Create local grid and set properties
    RigLocalGrid* localGrid = new RigLocalGrid( mainGrid );
    localGrid->setAsTempGrid( true );
    localGrid->setAssociatedWellPathName( lgrInfo.associatedWellPathName.toStdString() );
    localGrid->setGridId( lgrId );
    localGrid->setIndexToStartOfCells( totalCellCount );
    localGrid->setGridName( lgrInfo.name.toStdString() );
    localGrid->setCellCounts( lgrInfo.lgrCellCounts() );
    mainGrid->addLocalGrid( localGrid );

    size_t cellStartIndex = mainGrid->totalCellCount();
    size_t nodeStartIndex = mainGrid->nodes().size();

    // Resize global cell and node arrays
    {
        RigCell defaultCell;
        defaultCell.setHostGrid( localGrid );
        mainGrid->reservoirCells().resize( cellStartIndex + lgrCellCount, defaultCell );
        mainGrid->nodes().resize( nodeStartIndex + lgrCellCount * 8, cvf::Vec3d( 0, 0, 0 ) );
    }

    if ( lgrInfo.gradedRefinement )
    {
        localGrid->setRefinement( lgrInfo.gradedRefinement );
        fillGradedLgrCells( lgrInfo, mainGrid, localGrid, cellStartIndex, nodeStartIndex );
        localGrid->setParentGrid( mainGrid );
        return localGrid;
    }

    auto   refinement         = lgrInfo.refinement;
    size_t gridLocalCellIndex = 0;

    // Loop through all new LGR cells
    for ( size_t lgrK = 0; lgrK < lgrInfo.lgrCellCounts().k(); lgrK++ )
    {
        size_t mainK = lgrInfo.mainGridStartCell.k() + lgrK / refinement.k();

        for ( size_t lgrJ = 0; lgrJ < lgrInfo.lgrCellCounts().j(); lgrJ++ )
        {
            size_t mainJ = lgrInfo.mainGridStartCell.j() + lgrJ / refinement.j();

            for ( size_t lgrI = 0; lgrI < lgrInfo.lgrCellCounts().i(); lgrI++, gridLocalCellIndex++ )
            {
                size_t mainI = lgrInfo.mainGridStartCell.i() + lgrI / refinement.i();

                size_t mainCellIndex = mainGrid->cellIndexFromIJK( mainI, mainJ, mainK );
                auto&  mainGridCell  = mainGrid->cell( mainCellIndex );
                mainGridCell.setSubGrid( localGrid );

                RigCell& cell = mainGrid->cell( cellStartIndex + gridLocalCellIndex );
                cell.setGridLocalCellIndex( gridLocalCellIndex );
                cell.setParentCellIndex( mainCellIndex );

                // Corner coordinates
                {
                    size_t                    cIdx;
                    std::array<cvf::Vec3d, 8> vertices = mainGrid->cellCornerVertices( mainCellIndex );

                    std::vector<cvf::Vec3d> lgrCoords =
                        RiaCellDividingTools::createHexCornerCoords( vertices, refinement.i(), refinement.j(), refinement.k() );

                    size_t subI = lgrI % refinement.i();
                    size_t subJ = lgrJ % refinement.j();
                    size_t subK = lgrK % refinement.k();

                    size_t subIndex = subI + subJ * refinement.i() + subK * refinement.i() * refinement.j();

                    for ( cIdx = 0; cIdx < 8; ++cIdx )
                    {
                        auto& node = mainGrid->nodes()[nodeStartIndex + gridLocalCellIndex * 8 + cIdx];
                        node.set( lgrCoords[subIndex * 8 + cIdx] );
                        cell.cornerIndices()[cIdx] = nodeStartIndex + gridLocalCellIndex * 8 + cIdx;
                    }
                }
            }
        }
    }

    localGrid->setParentGrid( mainGrid );

    return localGrid;
}

//--------------------------------------------------------------------------------------------------
/// Cells and nodes of an LGR whose host cells are subdivided by the cumulative fractions of the graded refinement
//--------------------------------------------------------------------------------------------------
void RicCreateTemporaryLgrFeature::fillGradedLgrCells( const LgrInfo& lgrInfo,
                                                       RigMainGrid*   mainGrid,
                                                       RigLocalGrid*  localGrid,
                                                       size_t         cellStartIndex,
                                                       size_t         nodeStartIndex )
{
    const RigRefinement& refinement = *lgrInfo.gradedRefinement;
    const cvf::Vec3st    lgrCounts  = lgrInfo.lgrCellCounts();

    const size_t boxCountI = lgrInfo.mainGridEndCell.i() - lgrInfo.mainGridStartCell.i() + 1;
    const size_t boxCountJ = lgrInfo.mainGridEndCell.j() - lgrInfo.mainGridStartCell.j() + 1;
    const size_t boxCountK = lgrInfo.mainGridEndCell.k() - lgrInfo.mainGridStartCell.k() + 1;

    for ( size_t boxK = 0; boxK < boxCountK; boxK++ )
    {
        const auto& fractionsK = refinement.cumulativeFractions( RigRefinement::DimK, boxK );
        for ( size_t boxJ = 0; boxJ < boxCountJ; boxJ++ )
        {
            const auto& fractionsJ = refinement.cumulativeFractions( RigRefinement::DimJ, boxJ );
            for ( size_t boxI = 0; boxI < boxCountI; boxI++ )
            {
                const auto& fractionsI = refinement.cumulativeFractions( RigRefinement::DimI, boxI );

                size_t mainCellIndex = mainGrid->cellIndexFromIJK( lgrInfo.mainGridStartCell.i() + boxI,
                                                                   lgrInfo.mainGridStartCell.j() + boxJ,
                                                                   lgrInfo.mainGridStartCell.k() + boxK );
                mainGrid->cell( mainCellIndex ).setSubGrid( localGrid );

                const std::vector<cvf::Vec3d> subCellCoords =
                    RiaCellDividingTools::createHexCornerCoords( mainGrid->cellCornerVertices( mainCellIndex ), fractionsI, fractionsJ, fractionsK );

                for ( size_t subK = 0; subK < fractionsK.size(); subK++ )
                {
                    const size_t lgrK = refinement.cumulativeOffset( RigRefinement::DimK, boxK ) + subK;
                    for ( size_t subJ = 0; subJ < fractionsJ.size(); subJ++ )
                    {
                        const size_t lgrJ = refinement.cumulativeOffset( RigRefinement::DimJ, boxJ ) + subJ;
                        for ( size_t subI = 0; subI < fractionsI.size(); subI++ )
                        {
                            const size_t lgrI               = refinement.cumulativeOffset( RigRefinement::DimI, boxI ) + subI;
                            const size_t gridLocalCellIndex = lgrI + lgrJ * lgrCounts.i() + lgrK * lgrCounts.i() * lgrCounts.j();
                            const size_t subIndex = subI + subJ * fractionsI.size() + subK * fractionsI.size() * fractionsJ.size();

                            RigCell& cell = mainGrid->cell( cellStartIndex + gridLocalCellIndex );
                            cell.setGridLocalCellIndex( gridLocalCellIndex );
                            cell.setParentCellIndex( mainCellIndex );

                            for ( size_t cIdx = 0; cIdx < 8; ++cIdx )
                            {
                                const size_t nodeIndex       = nodeStartIndex + gridLocalCellIndex * 8 + cIdx;
                                mainGrid->nodes()[nodeIndex] = subCellCoords[subIndex * 8 + cIdx];
                                cell.cornerIndices()[cIdx]   = nodeIndex;
                            }
                        }
                    }
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RigGridBase* RicCreateTemporaryLgrFeature::createLgr( const LgrInfo& lgrInfo, RigEclipseCaseData* caseData )
{
    auto mainGrid = caseData->mainGrid();

    auto gridCount         = mainGrid->gridCount();
    auto mainGridCellCount = mainGrid->totalCellCount();

    auto localGrid = createLgrForGrid( lgrInfo, mainGrid );

    if ( auto activeInfo = caseData->activeCellInfo( RiaDefines::PorosityModelType::MATRIX_MODEL ) )
    {
        auto currentReservoirCellCount = activeInfo->reservoirCellCount();

        activeInfo->setGridCount( gridCount + 1 );
        activeInfo->setReservoirCellCount( currentReservoirCellCount + localGrid->cellCount() );

        size_t activeCellCount = 0;

        for ( size_t k = 0; k < localGrid->cellCountK(); k++ )
        {
            for ( size_t j = 0; j < localGrid->cellCountJ(); j++ )
            {
                for ( size_t i = 0; i < localGrid->cellCountI(); i++ )
                {
                    size_t cellIndex = localGrid->cellIndexFromIJK( i, j, k );
                    auto   gridCell  = localGrid->cell( cellIndex );
                    if ( gridCell.parentCellIndex() != cvf::UNDEFINED_SIZE_T )
                    {
                        auto resultIndex = activeInfo->cellResultIndex( ReservoirCellIndex( gridCell.parentCellIndex() ) );
                        if ( resultIndex.value() != cvf::UNDEFINED_SIZE_T )
                        {
                            activeInfo->setCellResultIndex( ReservoirCellIndex( cellIndex + mainGridCellCount ), resultIndex );
                            activeCellCount++;
                        }
                    }
                }
            }
        }

        activeInfo->setGridActiveCellCounts( gridCount, activeCellCount );
        activeInfo->computeDerivedData();
    }

    return localGrid;
}
