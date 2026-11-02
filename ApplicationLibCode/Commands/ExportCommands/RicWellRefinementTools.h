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

#pragma once

#include "RicWellLgrTools.h"

#include <QString>

#include <expected>

class RimEclipseView;
class RimWellRefinementRegion;

//==================================================================================================
///
/// Well Refinements: the approximate well path of a simulation well and amalgamated LGRs graded in rings
/// around its perforated part, kept as a refinement region of the view.
///
//==================================================================================================
namespace RicWellRefinementTools
{
// Creates the approximate well path of the region's well and the LGRs around it; the region keeps their boxes
std::expected<void, QString> applyRegion( RimWellRefinementRegion* region );

// The Well Refinements of the well in the view, or nullptr
RimWellRefinementRegion* findRegion( RimEclipseView* view, const QString& wellName );

// Creates or updates the Well Refinements of the well in the view
std::expected<RimWellRefinementRegion*, QString> createWellRefinement( RimEclipseView*                      view,
                                                                       const QString&                       wellName,
                                                                       int                                  timeStep,
                                                                       const QString&                       lgrName,
                                                                       const RicWellLgrTools::RingSettings& rings );
} // namespace RicWellRefinementTools
