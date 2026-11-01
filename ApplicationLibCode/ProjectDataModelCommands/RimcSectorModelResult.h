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

#include "cafPdmField.h"
#include "cafPdmObject.h"

#include <QString>

#include <vector>

//==================================================================================================
/// Files and cells of a sector model created for an OPM Flow run with USEFLUX
//==================================================================================================
class RimcSectorModelResult : public caf::PdmObject
{
    CAF_PDM_HEADER_INIT;

public:
    RimcSectorModelResult();

    caf::PdmField<QString>              m_deckFile;
    caf::PdmField<QString>              m_fluxFile;
    caf::PdmField<QString>              m_fluxnumFile;
    caf::PdmField<QString>              m_lgrName;
    caf::PdmField<int>                  m_cellCount;
    caf::PdmField<int>                  m_activeCellCount;
    caf::PdmField<std::vector<QString>> m_wells;
    caf::PdmField<std::vector<QString>> m_warnings;
};
