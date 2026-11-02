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

#include "RimRefinementRegion.h"

#include "ExportCommands/RicWellLgrTools.h"

#include "cafPdmField.h"

#include <vector>

//==================================================================================================
///
/// Well Refinements: an amalgamated LGR around a simulation well. The cells along its approximate
/// well path, from the first perforation to the last, are refined refinement * layers times, the
/// ring of cells around them refinement * (layers - 1) times, and so on, per direction. Each box of
/// equal refinement is a temporary LGR; the LGRs are amalgamated.
///
//==================================================================================================
class RimWellRefinementRegion : public RimRefinementRegion
{
    CAF_PDM_HEADER_INIT;

public:
    RimWellRefinementRegion();

    QString wellName() const;
    void    setWellName( const QString& wellName );
    int     timeStep() const;
    void    setTimeStep( int timeStep );
    QString lgrName() const;
    void    setLgrName( const QString& name );

    RicWellLgrTools::RingSettings ringSettings() const;
    void                          setRingSettings( const RicWellLgrTools::RingSettings& settings );
    void                          copySettingsFrom( const RimWellRefinementRegion& other );

    std::vector<RicWellLgrTools::RefinedBox> boxes() const;
    void                                     setBoxes( const std::vector<RicWellLgrTools::RefinedBox>& boxes );

    // Hides what does not apply before the region is created, e.g. in a dialog
    void setSettingsOnly( bool settingsOnly );

    std::vector<PreviewBox> previewBoxes() const override;
    bool                    isSectorGridRefinement() const override;

protected:
    void defineUiOrdering( QString uiConfigName, caf::PdmUiOrdering& uiOrdering ) override;
    void fieldChangedByUi( const caf::PdmFieldHandle* changedField, const QVariant& oldValue, const QVariant& newValue ) override;

private:
    caf::PdmField<QString> m_wellName;
    caf::PdmField<int>     m_timeStep;
    caf::PdmField<QString> m_lgrName;

    caf::PdmField<int>  m_lateralRefinement;
    caf::PdmField<int>  m_lateralLayers;
    caf::PdmField<bool> m_separateXY;
    caf::PdmField<int>  m_refinementX;
    caf::PdmField<int>  m_layersX;
    caf::PdmField<int>  m_refinementY;
    caf::PdmField<int>  m_layersY;
    caf::PdmField<int>  m_refinementZ;
    caf::PdmField<int>  m_layersZ;

    // Nine values per box: I, J, K of the first and the last host cell (zero-based), then the refinement in I, J, K
    caf::PdmField<std::vector<int>> m_boxes;
    caf::PdmField<bool>             m_updateLgrs;

    bool m_settingsOnly = false;
};
