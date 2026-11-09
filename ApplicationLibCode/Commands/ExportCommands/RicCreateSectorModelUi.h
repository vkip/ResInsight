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

#include "RicFluxSectorModelTools.h"

#include "cafPdmField.h"
#include "cafPdmObject.h"

class RimEclipseView;

//==================================================================================================
///
/// Settings of a sector model for an OPM Flow run with USEFLUX around a simulation well
///
//==================================================================================================
class RicCreateSectorModelUi : public caf::PdmObject
{
    CAF_PDM_HEADER_INIT;

public:
    RicCreateSectorModelUi();

    // The settings not tied to a well (make_flux and the sector extent), from an earlier sector model
    void setSharedSettings( const RicFluxSectorModelTools::Settings& settings );

    // Presets the folder and the files for the well, and describes the Well LGRs of the well in the view
    void setWell( RimEclipseView* view, const QString& wellName );

    RicFluxSectorModelTools::Settings settings() const;

private:
    void defineUiOrdering( QString uiConfigName, caf::PdmUiOrdering& uiOrdering ) override;

private:
    caf::PdmField<QString> m_wellName;
    caf::PdmField<QString> m_lgrInfo;
    caf::PdmField<QString> m_exportFolder;
    caf::PdmField<QString> m_deckFile;
    caf::PdmField<bool>    m_runMakeFlux;
    caf::PdmField<QString> m_makeFluxExecutable;
    caf::PdmField<int>     m_marginCells;
    caf::PdmField<int>     m_wellPaddingCells;
    caf::PdmField<int>     m_lgrPaddingCells;
};
