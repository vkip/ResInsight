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

#include "RicCreateApproximateWellPathFeature.h"

#include "RiaGuiApplication.h"
#include "RiaLogging.h"

#include "RicApproximateWellPathTools.h"
#include "RicEclipseWellFeatureImpl.h"

#include "Well/RigSimWellData.h"
#include "Well/RigWellResultFrame.h"
#include "Well/RigWellResultPoint.h"

#include "RimEclipseCase.h"
#include "RimEclipseView.h"
#include "RimSimWellInView.h"
#include "RimWellPath.h"

#include "Riu3DMainWindowTools.h"

#include <QAction>
#include <QMessageBox>

#include <algorithm>

CAF_CMD_SOURCE_INIT( RicCreateApproximateWellPathFeature, "RicCreateApproximateWellPathFeature" );

namespace
{
//--------------------------------------------------------------------------------------------------
/// The time step of the view, if the well has connections then, else the last time step with connections
//--------------------------------------------------------------------------------------------------
int connectionTimeStep( const RimSimWellInView& simWell )
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
} // namespace

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateApproximateWellPathFeature::onActionTriggered( bool isChecked )
{
    RimWellPath* lastCreated = nullptr;
    for ( RimSimWellInView* simWell : RicEclipseWellFeatureImpl::selectedWells() )
    {
        auto* eclipseCase = simWell->firstAncestorOrThisOfType<RimEclipseCase>();
        auto wellPath = RicApproximateWellPathTools::createApproximateWellPath( eclipseCase, simWell->name(), connectionTimeStep( *simWell ) );
        if ( !wellPath )
        {
            RiaLogging::error( wellPath.error().toStdString() );
            QMessageBox::warning( RiaGuiApplication::widgetToUseAsParent(), "Create Approximate Well Path", wellPath.error() );
            continue;
        }
        lastCreated = *wellPath;
    }

    if ( lastCreated ) Riu3DMainWindowTools::selectAsCurrentItem( lastCreated );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RicCreateApproximateWellPathFeature::setupActionLook( QAction* actionToSetup )
{
    actionToSetup->setText( "Create Approximate Well Path" );
    actionToSetup->setIcon( QIcon( ":/Well.svg" ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RicCreateApproximateWellPathFeature::isCommandEnabled() const
{
    return RicEclipseWellFeatureImpl::isAnyWellSelected();
}
