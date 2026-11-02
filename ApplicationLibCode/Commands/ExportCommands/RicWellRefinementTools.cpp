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

#include "RicWellRefinementTools.h"

#include "EclipseCommands/EclipseWell/RicApproximateWellPathTools.h"
#include "RicExportLgrFeature.h"

#include "RimEclipseCase.h"
#include "RimEclipseView.h"
#include "RimRefinementRegionCollection.h"
#include "RimWellPath.h"
#include "RimWellRefinementRegion.h"

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<void, QString> RicWellRefinementTools::applyRegion( RimWellRefinementRegion* region )
{
    if ( !region ) return std::unexpected( "No region." );
    auto* view = region->firstAncestorOrThisOfType<RimEclipseView>();
    if ( !view || !view->eclipseCase() ) return std::unexpected( "The region is not in a view of an Eclipse case." );
    RimEclipseCase* eclipseCase = view->eclipseCase();

    auto wellPath = RicApproximateWellPathTools::createApproximateWellPath( eclipseCase, region->wellName(), region->timeStep() );
    if ( !wellPath ) return std::unexpected( wellPath.error() );

    RicWellLgrTools::WellRefinementSettings settings;
    settings.lgrName = region->lgrName();
    settings.rings   = region->ringSettings();

    auto lgrs = RicWellLgrTools::buildWellRefinementLgrs( *eclipseCase, region->wellName(), **wellPath, settings );
    if ( !lgrs ) return std::unexpected( lgrs.error() );

    auto created = RicWellLgrTools::replaceTemporaryLgrsOfWell( eclipseCase, region->wellName(), std::move( *lgrs ) );
    if ( !created ) return std::unexpected( created.error() );

    std::vector<RicWellLgrTools::RefinedBox> boxes;
    for ( const auto& lgr : *created )
    {
        boxes.push_back(
            { lgr.mainGridStartCell,
              lgr.mainGridEndCell,
              { static_cast<int>( lgr.refinement.x() ), static_cast<int>( lgr.refinement.y() ), static_cast<int>( lgr.refinement.z() ) } } );
    }
    region->setBoxes( boxes );
    region->updateConnectedEditors();
    view->scheduleCreateDisplayModelAndRedraw();
    return {};
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimWellRefinementRegion* RicWellRefinementTools::findRegion( RimEclipseView* view, const QString& wellName )
{
    if ( !view || !view->refinementRegionCollection() ) return nullptr;
    for ( auto* region : view->refinementRegionCollection()->regions() )
    {
        auto* wellRegion = dynamic_cast<RimWellRefinementRegion*>( region );
        if ( wellRegion && wellRegion->wellName() == wellName ) return wellRegion;
    }
    return nullptr;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<RimWellRefinementRegion*, QString> RicWellRefinementTools::createWellRefinement( RimEclipseView* view,
                                                                                               const QString&  wellName,
                                                                                               int             timeStep,
                                                                                               const QString&  lgrName,
                                                                                               const RicWellLgrTools::RingSettings& rings )
{
    if ( !view || !view->refinementRegionCollection() ) return std::unexpected( "No view is available." );
    if ( wellName.isEmpty() ) return std::unexpected( "No well name specified." );

    auto*      region  = findRegion( view, wellName );
    const bool created = region == nullptr;
    if ( created )
    {
        region = new RimWellRefinementRegion();
        region->setWellName( wellName );
        view->refinementRegionCollection()->addRegion( region );
    }
    region->setTimeStep( timeStep );
    region->setLgrName( lgrName );
    region->setRingSettings( rings );

    if ( auto applied = applyRegion( region ); !applied )
    {
        if ( created ) view->refinementRegionCollection()->removeRegion( region );
        return std::unexpected( applied.error() );
    }

    view->refinementRegionCollection()->updateConnectedEditors();
    return region;
}
