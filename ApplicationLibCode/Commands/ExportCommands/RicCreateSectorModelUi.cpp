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

#include "RicCreateSectorModelUi.h"

#include "RicWellLgrTools.h"

#include "RigEclipseCaseData.h"
#include "RigMainGrid.h"

#include "RimEclipseCase.h"
#include "RimEclipseView.h"

#include "cafPdmUiFilePathEditor.h"
#include "cafPdmUiOrdering.h"

#include <QFileInfo>

CAF_PDM_SOURCE_INIT( RicCreateSectorModelUi, "RicCreateSectorModelUi" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RicCreateSectorModelUi::RicCreateSectorModelUi()
{
    CAF_PDM_InitObject( "Create Sector Model (USEFLUX)" );

    CAF_PDM_InitFieldNoDefault( &m_wellName, "WellName", "Well" );
    m_wellName.uiCapability()->setUiReadOnly( true );

    CAF_PDM_InitFieldNoDefault( &m_lgrInfo, "LgrInfo", "Well LGRs" );
    m_lgrInfo.uiCapability()->setUiReadOnly( true );
    m_lgrInfo.xmlCapability()->disableIO();

    CAF_PDM_InitFieldNoDefault( &m_exportFolder, "ExportFolder", "Export Folder" );
    m_exportFolder.uiCapability()->setUiEditorTypeName( caf::PdmUiFilePathEditor::uiEditorTypeName() );
    m_exportFolder.uiCapability()->setAttribute( caf::PdmUiFilePathEditor::Keys::SELECT_DIRECTORY, true );

    CAF_PDM_InitFieldNoDefault( &m_deckFile, "DeckFile", "Deck File" );
    m_deckFile.uiCapability()->setUiEditorTypeName( caf::PdmUiFilePathEditor::uiEditorTypeName() );
    m_deckFile.uiCapability()->setAttribute( caf::PdmUiFilePathEditor::Keys::FILE_SELECTION_FILTER,
                                             QString( "Eclipse/OPM Flow deck (*.DATA *.data);;All files (*)" ) );

    CAF_PDM_InitField( &m_runMakeFlux, "RunMakeFlux", true, "Build the FLUX File with make_flux" );

    CAF_PDM_InitFieldNoDefault( &m_makeFluxExecutable, "MakeFluxExecutable", "make_flux" );
    m_makeFluxExecutable.uiCapability()->setUiEditorTypeName( caf::PdmUiFilePathEditor::uiEditorTypeName() );

    RicFluxSectorModelTools::Settings defaults;
    CAF_PDM_InitField( &m_marginCells, "MarginCells", defaults.marginCells, "Margin Cells" );
    m_marginCells.uiCapability()->setUiToolTip( "Main grid cells around the well connections, laterally" );
    CAF_PDM_InitField( &m_wellPaddingCells, "WellPaddingCells", defaults.wellPaddingCells, "Well Padding Cells" );
    m_wellPaddingCells.uiCapability()->setUiToolTip(
        "Cells around another well that the sector boundary would cut, when the sector is enlarged to hold it" );
    CAF_PDM_InitField( &m_lgrPaddingCells, "LgrPaddingCells", defaults.lgrPaddingCells, "LGR Padding Cells" );
    m_lgrPaddingCells.uiCapability()->setUiToolTip( "Cells between the LGRs around the well and the sector boundary" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateSectorModelUi::setSharedSettings( const RicFluxSectorModelTools::Settings& settings )
{
    m_runMakeFlux        = settings.runMakeFlux;
    m_makeFluxExecutable = settings.makeFluxExecutable;
    m_marginCells        = settings.marginCells;
    m_wellPaddingCells   = settings.wellPaddingCells;
    m_lgrPaddingCells    = settings.lgrPaddingCells;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateSectorModelUi::setWell( RimEclipseView* view, const QString& wellName )
{
    m_wellName = wellName;

    RimEclipseCase* eclipseCase = view ? view->eclipseCase() : nullptr;
    if ( !eclipseCase ) return;

    const QString gridFileName = eclipseCase->gridFileName();
    if ( m_exportFolder().isEmpty() ) m_exportFolder = RicFluxSectorModelTools::defaultExportFolder( gridFileName, wellName );

    if ( m_deckFile().isEmpty() )
    {
        const QFileInfo gridInfo( gridFileName );
        for ( const QString& suffix : { ".DATA", ".data", ".Data" } )
        {
            const QString candidate = gridInfo.absolutePath() + "/" + gridInfo.completeBaseName() + suffix;
            if ( QFileInfo::exists( candidate ) )
            {
                m_deckFile = candidate;
                break;
            }
        }
    }

    if ( m_makeFluxExecutable().isEmpty() ) m_makeFluxExecutable = RicFluxSectorModelTools::defaultMakeFluxExecutable();

    const auto* caseData = eclipseCase->eclipseCaseData();
    const auto lgrCount = caseData && caseData->mainGrid() ? RicWellLgrTools::temporaryLgrsForWell( *caseData->mainGrid(), wellName ).size()
                                                           : size_t{ 0 };
    m_lgrInfo           = lgrCount > 0 ? QString( "%1 LGRs, added to the sector deck" ).arg( lgrCount )
                                       : QString( "None. Use Create Well LGR first for a refined sector" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RicFluxSectorModelTools::Settings RicCreateSectorModelUi::settings() const
{
    RicFluxSectorModelTools::Settings settings;
    settings.wellName           = m_wellName();
    settings.exportFolder       = m_exportFolder();
    settings.deckFile           = m_deckFile();
    settings.runMakeFlux        = m_runMakeFlux();
    settings.makeFluxExecutable = m_makeFluxExecutable();
    settings.marginCells        = m_marginCells();
    settings.wellPaddingCells   = m_wellPaddingCells();
    settings.lgrPaddingCells    = m_lgrPaddingCells();
    return settings;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateSectorModelUi::defineUiOrdering( QString uiConfigName, caf::PdmUiOrdering& uiOrdering )
{
    uiOrdering.add( &m_wellName );
    uiOrdering.add( &m_lgrInfo );

    auto* files = uiOrdering.addNewGroup( "Files" );
    files->add( &m_exportFolder );
    files->add( &m_deckFile );
    files->add( &m_runMakeFlux );
    files->add( &m_makeFluxExecutable );
    m_makeFluxExecutable.uiCapability()->setUiReadOnly( !m_runMakeFlux() );

    auto* extent = uiOrdering.addNewGroup( "Sector Extent" );
    extent->add( &m_marginCells );
    extent->add( &m_wellPaddingCells );
    extent->add( &m_lgrPaddingCells );

    uiOrdering.skipRemainingFields( true );
}
