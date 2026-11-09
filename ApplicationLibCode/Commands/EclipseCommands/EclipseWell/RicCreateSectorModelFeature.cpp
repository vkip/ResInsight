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

#include "RicCreateSectorModelFeature.h"

#include "RiaLogging.h"

#include "ExportCommands/RicCreateSectorModelUi.h"
#include "ExportCommands/RicFluxSectorModelTools.h"
#include "RicEclipseWellFeatureImpl.h"

#include "RimEclipseView.h"
#include "RimSimWellInView.h"

#include "Riu3DMainWindowTools.h"
#include "RiuMessageDialog.h"

#include "cafPdmUiPropertyViewDialog.h"

#include <QAction>
#include <QApplication>
#include <QMessageBox>

#include <optional>

CAF_CMD_SOURCE_INIT( RicCreateSectorModelFeature, "RicCreateSectorModelFeature" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateSectorModelFeature::onActionTriggered( bool isChecked )
{
    const auto selection = RicEclipseWellFeatureImpl::selectedWells();
    if ( selection.empty() ) return;
    RimSimWellInView* simWell = selection.front();
    auto*             view    = simWell->firstAncestorOrThisOfType<RimEclipseView>();
    if ( !view ) return;

    // make_flux and the sector extent are kept from one sector model to the next in the session
    static std::optional<RicFluxSectorModelTools::Settings> sharedSettings;

    RicCreateSectorModelUi ui;
    if ( sharedSettings ) ui.setSharedSettings( *sharedSettings );
    ui.setWell( view, simWell->name() );

    caf::PdmUiPropertyViewDialog dialog( Riu3DMainWindowTools::mainWindowWidget(), &ui, "Create Sector Model (USEFLUX)", "" );
    dialog.resize( QSize( 560, 420 ) );
    if ( dialog.exec() != QDialog::Accepted ) return;

    const auto settings = ui.settings();
    sharedSettings      = settings;

    QApplication::setOverrideCursor( QCursor( Qt::WaitCursor ) );
    const auto result = RicFluxSectorModelTools::createSectorModel( view, settings );
    QApplication::restoreOverrideCursor();

    if ( !result )
    {
        RiaLogging::error( result.error().toStdString() );
        RiuMessageDialog::showError( Riu3DMainWindowTools::mainWindowWidget(), "Create Sector Model", result.error() );
        return;
    }

    QString text = QString( "Sector model around well %1: %2 cells (%3 active), wells %4.\n\nDeck: %5" )
                       .arg( settings.wellName )
                       .arg( result->sectorCellCount )
                       .arg( result->activeSectorCellCount )
                       .arg( result->wellsInSector.join( ", " ) )
                       .arg( result->deckFile );
    if ( !result->fluxFile.isEmpty() && settings.runMakeFlux ) text += QString( "\nFLUX file: %1" ).arg( result->fluxFile );
    if ( !result->lgrNames.isEmpty() ) text += QString( "\nLGRs: %1" ).arg( result->lgrNames.size() );
    if ( !result->warnings.isEmpty() ) text += "\n\n" + result->warnings.join( "\n\n" );

    QMessageBox::information( Riu3DMainWindowTools::mainWindowWidget(), "Create Sector Model", text );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateSectorModelFeature::setupActionLook( QAction* actionToSetup )
{
    actionToSetup->setText( "Create Sector Model (USEFLUX)" );
    actionToSetup->setIcon( QIcon( ":/Save.svg" ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RicCreateSectorModelFeature::isCommandEnabled() const
{
    return RicEclipseWellFeatureImpl::isAnyWellSelected();
}
