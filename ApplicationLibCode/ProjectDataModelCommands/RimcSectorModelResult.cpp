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

#include "RimcSectorModelResult.h"

#include "cafPdmFieldScriptingCapability.h"
#include "cafPdmObjectScriptingCapability.h"

CAF_PDM_SOURCE_INIT( RimcSectorModelResult, "SectorModelResult" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimcSectorModelResult::RimcSectorModelResult()
{
    CAF_PDM_InitScriptableObject( "Sector Model Result", "", "", "Files and cells of a sector model for an OPM Flow run with USEFLUX" );

    CAF_PDM_InitScriptableField( &m_deckFile, "DeckFile", QString(), "Deck File", "", "", "The sector deck, with USEFLUX and FLUXREG" );
    CAF_PDM_InitScriptableField( &m_fluxFile, "FluxFile", QString(), "Flux File", "", "", "The boundary data written by make_flux" );
    CAF_PDM_InitScriptableField( &m_fluxnumFile, "FluxnumFile", QString(), "FLUXNUM File", "", "", "The sector as a FLUXNUM region" );
    CAF_PDM_InitScriptableField( &m_lgrName, "LgrName", QString(), "LGR Name", "", "", "The LGR around the well in the sector deck, if any" );
    CAF_PDM_InitScriptableField( &m_cellCount, "CellCount", 0, "Cell Count", "", "", "Main grid cells in the sector" );
    CAF_PDM_InitScriptableField( &m_activeCellCount, "ActiveCellCount", 0, "Active Cell Count", "", "", "Active main grid cells in the sector" );
    CAF_PDM_InitScriptableFieldNoDefault( &m_wells, "Wells", "Wells", "", "", "Wells with connections in the sector" );
    CAF_PDM_InitScriptableFieldNoDefault( &m_warnings, "Warnings", "Warnings", "", "", "What was adjusted to make the sector model" );
}
