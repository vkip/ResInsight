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

#include "RicCreateWellLgrFeature.h"

#include "RiaLogging.h"

#include "ExportCommands/RicWellRefinementTools.h"
#include "RicApproximateWellPathTools.h"
#include "RicEclipseWellFeatureImpl.h"

#include "RimEclipseView.h"
#include "RimSimWellInView.h"
#include "RimWellRefinementRegion.h"

#include "Riu3DMainWindowTools.h"
#include "RiuMessageDialog.h"

#include "cafPdmUiPropertyViewDialog.h"

#include <QAction>

CAF_CMD_SOURCE_INIT( RicCreateWellLgrFeature, "RicCreateWellLgrFeature" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateWellLgrFeature::onActionTriggered( bool isChecked )
{
    const auto selection = RicEclipseWellFeatureImpl::selectedWells();
    if ( selection.empty() ) return;
    RimSimWellInView* simWell = selection.front();
    auto*             view    = simWell->firstAncestorOrThisOfType<RimEclipseView>();
    if ( !view ) return;

    // The settings are edited on a region outside the project, preset from the well's current one
    RimWellRefinementRegion settings;
    settings.setSettingsOnly( true );
    settings.setWellName( simWell->name() );
    if ( auto* existing = RicWellRefinementTools::findRegion( view, simWell->name() ) ) settings.copySettingsFrom( *existing );
    settings.setTimeStep( RicApproximateWellPathTools::viewTimeStepWithConnections( *simWell ) );

    caf::PdmUiPropertyViewDialog dialog( Riu3DMainWindowTools::mainWindowWidget(), &settings, "Create Well LGR", "" );
    dialog.resize( QSize( 420, 480 ) );
    if ( dialog.exec() != QDialog::Accepted ) return;

    auto region =
        RicWellRefinementTools::createWellRefinement( view, simWell->name(), settings.timeStep(), settings.lgrName(), settings.ringSettings() );
    if ( !region )
    {
        RiaLogging::error( region.error().toStdString() );
        RiuMessageDialog::showError( Riu3DMainWindowTools::mainWindowWidget(), "Create Well LGR", region.error() );
        return;
    }

    Riu3DMainWindowTools::selectAsCurrentItem( *region );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateWellLgrFeature::setupActionLook( QAction* actionToSetup )
{
    actionToSetup->setText( "Create Well LGR" );
    actionToSetup->setIcon( QIcon( ":/CellFilter_Range.png" ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RicCreateWellLgrFeature::isCommandEnabled() const
{
    return RicEclipseWellFeatureImpl::isAnyWellSelected();
}
