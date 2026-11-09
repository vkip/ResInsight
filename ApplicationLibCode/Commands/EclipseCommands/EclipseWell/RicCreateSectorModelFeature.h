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

#include "cafCmdFeature.h"

//==================================================================================================
///
/// Creates a sector model for an OPM Flow run with USEFLUX around the selected simulation well: the
/// sector as a polygon filter of the view, its FLUXNUM region, the FLUX file from make_flux, and the
/// sector deck, with the Well LGRs of the well if it has any
///
//==================================================================================================
class RicCreateSectorModelFeature : public caf::CmdFeature
{
    CAF_CMD_HEADER_INIT;

protected:
    void onActionTriggered( bool isChecked ) override;
    void setupActionLook( QAction* actionToSetup ) override;
    bool isCommandEnabled() const override;
};
