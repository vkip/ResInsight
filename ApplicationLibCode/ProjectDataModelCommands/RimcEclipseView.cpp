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

#include "RimcEclipseView.h"

#include "ExportCommands/RicEclipseCellResultToFileImpl.h"
#include "ExportCommands/RicFluxSectorModelTools.h"
#include "ExportCommands/RicSaveEclipseInputVisibleCellsFeature.h"
#include "ExportCommands/RicWellRefinementTools.h"

#include "RiaViewRedrawScheduler.h"

#include "RigResultAccessor.h"
#include "RigResultAccessorFactory.h"

#include "RimEclipseCase.h"
#include "RimEclipseCellColors.h"
#include "RimEclipseView.h"
#include "RimFaultDistance.h"
#include "RimFaultDistanceCollection.h"
#include "RimFaultInView.h"
#include "RimFaultInViewCollection.h"
#include "RimWellRefinementRegion.h"

#include "RimcSectorModelResult.h"

#include "cafPdmAbstractFieldScriptingCapability.h"
#include "cafPdmFieldScriptingCapability.h"
#include "cafUtils.h"

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimEclipseView, RimcEclipseView_addFaultDistance, "add_fault_distance" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimcEclipseView_addFaultDistance::RimcEclipseView_addFaultDistance( caf::PdmObjectHandle* self )
    : caf::PdmObjectCreationMethod( self )
{
    CAF_PDM_InitObject( "Add Fault Distance", "", "", "Create a FAULTDIST cell result for a chosen subset of faults" );

    CAF_PDM_InitScriptableFieldNoDefault( &m_resultName, "Name", "Name (default FAULTDIST<n> if empty)" );
    CAF_PDM_InitScriptableFieldNoDefault( &m_faults, "Faults", "Faults to include (empty = all)" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimcEclipseView_addFaultDistance::execute()
{
    auto* eclipseView = self<RimEclipseView>();
    if ( !eclipseView ) return std::unexpected( QString( "No view" ) );

    auto* distanceCollection = eclipseView->faultDistanceCollection();
    if ( !distanceCollection ) return std::unexpected( QString( "No fault distance results collection" ) );

    auto* newResult = distanceCollection->addResult();
    if ( !newResult ) return std::unexpected( QString( "Failed to create fault distance result" ) );

    if ( !m_resultName().isEmpty() ) newResult->setResultName( m_resultName() );

    std::vector<RimFaultInView*> selected = m_faults.ptrReferencedObjectsByType();
    if ( selected.empty() && eclipseView->faultCollection() )
    {
        // Default to all faults, but leave the ResInsight-generated faults unticked.
        selected = eclipseView->faultCollection()->faults();
        std::erase_if( selected, []( RimFaultInView* fault ) { return fault && fault->isGeneratedFault(); } );
    }

    newResult->setSelectedFaults( selected );

    // When created from Python, always trigger the calculation (the UI uses an explicit Generate button instead).
    newResult->compute();

    eclipseView->updateConnectedEditors();

    return newResult;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RimcEclipseView_addFaultDistance::classKeywordReturnedType() const
{
    return RimFaultDistance::classKeywordStatic();
}

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimEclipseView, RimEclipseView_exportCurrentProperty, "exportCurrentProperty" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimEclipseView_exportCurrentProperty::RimEclipseView_exportCurrentProperty( caf::PdmObjectHandle* self )
    : caf::PdmVoidObjectMethod( self )
{
    CAF_PDM_InitObject( "Export Current Property", "", "", "Export the cell result currently shown in the view to a GRDECL style text file" );

    CAF_PDM_InitScriptableField( &m_exportFile, "ExportFile", QString(), "Export File", "", "", "Full path of the file to write" );
    CAF_PDM_InitScriptableField( &m_undefinedValue, "UndefinedValue", 0.0, "Undefined Value", "", "", "Value written for undefined cells" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportCurrentProperty::setExportFile( const QString& exportFile )
{
    m_exportFile = exportFile;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportCurrentProperty::setUndefinedValue( double undefinedValue )
{
    m_undefinedValue = undefinedValue;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimEclipseView_exportCurrentProperty::execute()
{
    auto* view = self<RimEclipseView>();
    if ( !view ) return std::unexpected( "No view is available." );

    if ( m_exportFile().isEmpty() ) return std::unexpected( "No export file specified." );

    RimEclipseCase* eclipseCase = view->eclipseCase();
    if ( !eclipseCase || !eclipseCase->eclipseCaseData() ) return std::unexpected( "The view has no case data." );

    cvf::ref<RigResultAccessor> resultAccessor = currentPropertyAccessor( view );

    const QString propertyName = view->cellResult()->resultVariableUiShortName();

    if ( resultAccessor.isNull() )
    {
        return std::unexpected( QString( "Could not find property '%1' at time step %2 in view '%3'" )
                                    .arg( propertyName )
                                    .arg( view->currentTimeStep() )
                                    .arg( view->name() ) );
    }

    const bool writeEchoKeywords = false;
    QString    errorMessage;
    if ( !RicEclipseCellResultToFileImpl::writeResultToTextFile( m_exportFile(),
                                                                 eclipseCase->eclipseCaseData(),
                                                                 resultAccessor.p(),
                                                                 propertyName,
                                                                 m_undefinedValue(),
                                                                 "exportProperty",
                                                                 writeEchoKeywords,
                                                                 &errorMessage ) )
    {
        return std::unexpected( errorMessage );
    }

    return nullptr;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
cvf::ref<RigResultAccessor> RimEclipseView_exportCurrentProperty::currentPropertyAccessor( const RimEclipseView* view )
{
    if ( !view || !view->eclipseCase() || !view->eclipseCase()->eclipseCaseData() ) return nullptr;

    const int mainGridIndex = 0;
    return RigResultAccessorFactory::createFromResultDefinition( view->eclipseCase()->eclipseCaseData(),
                                                                 mainGridIndex,
                                                                 view->currentTimeStep(),
                                                                 view->cellResult() );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RimEclipseView_exportCurrentProperty::hasCurrentProperty( const RimEclipseView* view )
{
    return currentPropertyAccessor( view ).notNull();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RimEclipseView_exportCurrentProperty::defaultFileBaseName( const RimEclipseView* view )
{
    if ( !view || !view->eclipseCase() ) return {};

    const QString propertyName = view->cellResult()->resultVariableUiShortName();
    const QString fileName =
        QString( "%1-%2-T%3-%4" ).arg( view->eclipseCase()->caseUserDescription() ).arg( view->name() ).arg( view->currentTimeStep() ).arg( propertyName );

    return caf::Utils::makeValidFileBasename( fileName );
}

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimEclipseView, RimEclipseView_exportVisibleCells, "exportVisibleCells" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimEclipseView_exportVisibleCells::RimEclipseView_exportVisibleCells( caf::PdmObjectHandle* self )
    : caf::PdmVoidObjectMethod( self )
{
    CAF_PDM_InitObject( "Export Visible Cells", "", "", "Export a GRDECL keyword with one value per cell based on cell visibility in the view" );

    CAF_PDM_InitScriptableField( &m_exportFile, "ExportFile", QString(), "Export File", "", "", "Full path of the file to write" );
    CAF_PDM_InitScriptableField( &m_exportKeyword,
                                 "ExportKeyword",
                                 RicSaveEclipseInputVisibleCellsUi::FLUXNUM,
                                 "Export Keyword",
                                 "",
                                 "",
                                 "GRDECL keyword to write" );
    CAF_PDM_InitScriptableField( &m_visibleActiveCellsValue,
                                 "VisibleActiveCellsValue",
                                 1,
                                 "Visible Active Cells Value",
                                 "",
                                 "",
                                 "Value for visible active cells" );
    CAF_PDM_InitScriptableField( &m_hiddenActiveCellsValue,
                                 "HiddenActiveCellsValue",
                                 0,
                                 "Hidden Active Cells Value",
                                 "",
                                 "",
                                 "Value for hidden active cells" );
    CAF_PDM_InitScriptableField( &m_inactiveCellsValue, "InactiveCellsValue", 0, "Inactive Cells Value", "", "", "Value for inactive cells" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportVisibleCells::setExportFile( const QString& exportFile )
{
    m_exportFile = exportFile;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportVisibleCells::setExportKeyword( RicSaveEclipseInputVisibleCellsUi::ExportKeyword exportKeyword )
{
    m_exportKeyword = exportKeyword;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportVisibleCells::setVisibleActiveCellsValue( int value )
{
    m_visibleActiveCellsValue = value;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportVisibleCells::setHiddenActiveCellsValue( int value )
{
    m_hiddenActiveCellsValue = value;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimEclipseView_exportVisibleCells::setInactiveCellsValue( int value )
{
    m_inactiveCellsValue = value;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimEclipseView_exportVisibleCells::execute()
{
    auto* view = self<RimEclipseView>();
    if ( !view ) return std::unexpected( "No view is available." );

    if ( m_exportFile().isEmpty() ) return std::unexpected( "No export file specified." );

    RiaViewRedrawScheduler::instance()->clearViewsScheduledForUpdate();

    RicSaveEclipseInputVisibleCellsUi exportSettings;
    exportSettings.exportFilename          = m_exportFile();
    exportSettings.exportKeyword           = m_exportKeyword();
    exportSettings.visibleActiveCellsValue = m_visibleActiveCellsValue();
    exportSettings.hiddenActiveCellsValue  = m_hiddenActiveCellsValue();
    exportSettings.inactiveCellsValue      = m_inactiveCellsValue();

    if ( auto result = RicSaveEclipseInputVisibleCellsFeature::exportVisibleCells( view, exportSettings ); !result )
    {
        return std::unexpected( result.error() );
    }

    return nullptr;
}

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimEclipseView, RimEclipseView_createSectorModelForWell, "createSectorModelForWell" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimEclipseView_createSectorModelForWell::RimEclipseView_createSectorModelForWell( caf::PdmObjectHandle* self )
    : caf::PdmObjectMethod( self, PdmObjectMethod::NullPointerType::NULL_IS_INVALID, PdmObjectMethod::ResultType::PERSISTENT_FALSE )
{
    CAF_PDM_InitObject( "Create Sector Model for Well",
                        "",
                        "",
                        "Create a sector model around a simulation well, and optionally more wells, for an OPM Flow run "
                        "with USEFLUX. The sector is shown as a polygon filter in the view, replacing one of an earlier "
                        "sector around the same wells, and exported as a FLUXNUM region, make_flux builds its boundary "
                        "data, and the deck is copied with USEFLUX, FLUXREG and the region added. The Well LGRs of the "
                        "well are added to the sector deck, with the well moved into them. With export_files false, the "
                        "sector is only shown in the view." );

    CAF_PDM_InitScriptableField( &m_wellName, "WellName", QString(), "Well Name", "", "", "Name of the simulation well" );
    CAF_PDM_InitScriptableField( &m_exportFolder,
                                 "ExportFolder",
                                 QString(),
                                 "Export Folder",
                                 "",
                                 "",
                                 "Folder for the sector files. Default: 'sector_<well>' next to the grid file" );
    CAF_PDM_InitScriptableField( &m_deckFile,
                                 "DeckFile",
                                 QString(),
                                 "Deck File",
                                 "",
                                 "",
                                 "Deck of the case. Default: the DATA file next to the grid file; if there is none, the user is asked" );
    CAF_PDM_InitScriptableField( &m_makeFluxExecutable,
                                 "MakeFluxExecutable",
                                 QString(),
                                 "make_flux Executable",
                                 "",
                                 "",
                                 "Path to make_flux. Default: next to the OPM Flow executable of the preferences, else in PATH" );
    CAF_PDM_InitScriptableField( &m_marginCells, "MarginCells", 3, "Margin Cells", "", "", "Main grid cells around the well connections, laterally" );
    CAF_PDM_InitScriptableField( &m_wellPaddingCells,
                                 "WellPaddingCells",
                                 2,
                                 "Well Padding Cells",
                                 "",
                                 "",
                                 "Cells around another well that the sector boundary would cut, when the sector is enlarged to hold it" );
    CAF_PDM_InitScriptableField( &m_lgrPaddingCells,
                                 "LgrPaddingCells",
                                 2,
                                 "LGR Padding Cells",
                                 "",
                                 "",
                                 "Cells between the LGR around the well and the sector boundary" );
    CAF_PDM_InitScriptableField( &m_runMakeFlux, "RunMakeFlux", true, "Run make_flux", "", "", "Build the FLUX file with make_flux" );
    CAF_PDM_InitScriptableFieldNoDefault( &m_additionalWells,
                                          "AdditionalWells",
                                          "Additional Wells",
                                          "",
                                          "",
                                          "Further simulation wells the sector is built around, with the same margin" );
    CAF_PDM_InitScriptableField( &m_exportFiles,
                                 "ExportFiles",
                                 true,
                                 "Export Files",
                                 "",
                                 "",
                                 "Export the FLUXNUM, FLUX and deck files. False: only show the sector in the view" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimEclipseView_createSectorModelForWell::execute()
{
    auto* view = self<RimEclipseView>();

    RicFluxSectorModelTools::Settings settings;
    settings.wellName           = m_wellName();
    settings.exportFolder       = m_exportFolder();
    settings.deckFile           = m_deckFile();
    settings.makeFluxExecutable = m_makeFluxExecutable();
    settings.marginCells        = m_marginCells();
    settings.wellPaddingCells   = m_wellPaddingCells();
    settings.lgrPaddingCells    = m_lgrPaddingCells();
    settings.runMakeFlux        = m_runMakeFlux();
    settings.exportFiles        = m_exportFiles();
    for ( const auto& well : m_additionalWells() )
    {
        settings.additionalWells.push_back( well );
    }

    auto result = RicFluxSectorModelTools::createSectorModel( view, settings );
    if ( !result ) return std::unexpected( result.error() );

    auto* output              = new RimcSectorModelResult();
    output->m_deckFile        = result->deckFile;
    output->m_fluxFile        = result->fluxFile;
    output->m_fluxnumFile     = result->fluxnumFile;
    output->m_lgrName         = result->lgrNames.join( "," );
    output->m_cellCount       = result->sectorCellCount;
    output->m_activeCellCount = result->activeSectorCellCount;
    output->m_wells           = std::vector<QString>( result->wellsInSector.begin(), result->wellsInSector.end() );
    output->m_warnings        = std::vector<QString>( result->warnings.begin(), result->warnings.end() );
    return output;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RimEclipseView_createSectorModelForWell::classKeywordReturnedType() const
{
    return RimcSectorModelResult::classKeywordStatic();
}

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimEclipseView, RimEclipseView_createWellRefinement, "createWellRefinement" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimEclipseView_createWellRefinement::RimEclipseView_createWellRefinement( caf::PdmObjectHandle* self )
    : caf::PdmObjectCreationMethod( self )
{
    CAF_PDM_InitObject( "Create Well Refinement",
                        "",
                        "",
                        "Create the approximate well path of a simulation well and amalgamated LGRs around its perforated part, as "
                        "Well Refinements of the view. The path cells are refined refinement * layers times laterally, the ring of "
                        "cells around them refinement * (layers - 1) times, and so on; vertically likewise in layers above and below. "
                        "Previous Well Refinements of the well are replaced." );

    CAF_PDM_InitScriptableField( &m_wellName, "WellName", QString(), "Well Name", "", "", "Name of the simulation well" );
    CAF_PDM_InitScriptableField( &m_timeStep,
                                 "TimeStep",
                                 -1,
                                 "Time Step",
                                 "",
                                 "",
                                 "Zero-based time step of the connections. -1 uses the last time step where the well has connections" );
    CAF_PDM_InitScriptableField( &m_lgrName, "LgrName", QString(), "LGR Name", "", "", "Base name of the LGRs. Default: 'L' followed by the well name" );
    CAF_PDM_InitScriptableField( &m_lateralRefinement, "LateralRefinement", 3, "Lateral Refinement", "", "", "nref in X and Y" );
    CAF_PDM_InitScriptableField( &m_lateralLayers, "LateralLayers", 3, "Lateral Layers", "", "", "nlayers in X and Y" );
    CAF_PDM_InitScriptableField( &m_refinementX, "RefinementX", 0, "Refinement X", "", "", "nref in X, if not 0" );
    CAF_PDM_InitScriptableField( &m_layersX, "LayersX", 0, "Layers X", "", "", "nlayers in X, if not 0" );
    CAF_PDM_InitScriptableField( &m_refinementY, "RefinementY", 0, "Refinement Y", "", "", "nref in Y, if not 0" );
    CAF_PDM_InitScriptableField( &m_layersY, "LayersY", 0, "Layers Y", "", "", "nlayers in Y, if not 0" );
    CAF_PDM_InitScriptableField( &m_refinementZ, "RefinementZ", 1, "Refinement Z", "", "", "nref in Z" );
    CAF_PDM_InitScriptableField( &m_layersZ, "LayersZ", 1, "Layers Z", "", "", "nlayers in Z" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimEclipseView_createWellRefinement::execute()
{
    auto* view = self<RimEclipseView>();
    if ( !view ) return std::unexpected( "No view is available." );

    RicWellLgrTools::RingSettings rings;
    rings.refinement = { m_refinementX() > 0 ? m_refinementX() : m_lateralRefinement(),
                         m_refinementY() > 0 ? m_refinementY() : m_lateralRefinement(),
                         m_refinementZ() };
    rings.layers = { m_layersX() > 0 ? m_layersX() : m_lateralLayers(), m_layersY() > 0 ? m_layersY() : m_lateralLayers(), m_layersZ() };

    auto region = RicWellRefinementTools::createWellRefinement( view, m_wellName(), m_timeStep(), m_lgrName(), rings );
    if ( !region ) return std::unexpected( region.error() );
    return *region;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RimEclipseView_createWellRefinement::classKeywordReturnedType() const
{
    return RimWellRefinementRegion::classKeywordStatic();
}
