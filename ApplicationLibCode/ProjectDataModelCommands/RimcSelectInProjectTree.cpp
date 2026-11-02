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

#include "RimcSelectInProjectTree.h"

#include "RimSimWellInView.h"
#include "RimWellPath.h"

#include "Riu3DMainWindowTools.h"

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimSimWellInView, RimcSimWellInView_select, "select" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimcSimWellInView_select::RimcSimWellInView_select( caf::PdmObjectHandle* self )
    : caf::PdmVoidObjectMethod( self )
{
    CAF_PDM_InitObject( "Select", "", "", "Select the simulation well in the project tree, making its view the active view" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimcSimWellInView_select::execute()
{
    Riu3DMainWindowTools::selectAsCurrentItem( self<RimSimWellInView>() );
    return nullptr;
}

CAF_PDM_OBJECT_METHOD_SOURCE_INIT( RimWellPath, RimcWellPath_select, "select" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimcWellPath_select::RimcWellPath_select( caf::PdmObjectHandle* self )
    : caf::PdmVoidObjectMethod( self )
{
    CAF_PDM_InitObject( "Select", "", "", "Select the well path in the project tree" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<caf::PdmObjectHandle*, QString> RimcWellPath_select::execute()
{
    Riu3DMainWindowTools::selectAsCurrentItem( self<RimWellPath>() );
    return nullptr;
}
