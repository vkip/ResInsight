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

#include "RimWellRefinementRegion.h"

#include "RiaLogging.h"

#include "ExportCommands/RicWellRefinementTools.h"

#include "RigUniformRefinement.h"

#include "Rim3dView.h"

#include "RiuMessageDialog.h"

#include "cafPdmFieldScriptingCapability.h"
#include "cafPdmObjectScriptingCapability.h"
#include "cafPdmUiPushButtonEditor.h"

CAF_PDM_SOURCE_INIT( RimWellRefinementRegion, "WellRefinementRegion" );

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RimWellRefinementRegion::RimWellRefinementRegion()
{
    CAF_PDM_InitScriptableObjectWithNameAndComment( "Well Refinements",
                                                    ":/Well.svg",
                                                    "",
                                                    "",
                                                    "WellRefinementRegion",
                                                    "An amalgamated LGR graded in rings around a simulation well" );

    CAF_PDM_InitScriptableField( &m_wellName, "WellName", QString(), "Well" );
    m_wellName.uiCapability()->setUiReadOnly( true );
    CAF_PDM_InitScriptableField( &m_timeStep,
                                 "TimeStep",
                                 -1,
                                 "Time Step",
                                 "",
                                 "Zero-based time step of the connections. -1 uses the last time step where the well has connections" );
    CAF_PDM_InitScriptableField( &m_lgrName,
                                 "LgrName",
                                 QString(),
                                 "LGR Name",
                                 "",
                                 "Base name of the LGRs, followed by a number. Default: 'L' followed by the well name" );

    CAF_PDM_InitScriptableField( &m_lateralRefinement, "LateralRefinement", 3, "Refinement (nref)" );
    CAF_PDM_InitScriptableField( &m_lateralLayers, "LateralLayers", 3, "Layers (nlayers)" );
    CAF_PDM_InitScriptableField( &m_separateXY, "SeparateXy", false, "Separate X and Y" );
    CAF_PDM_InitScriptableField( &m_refinementX, "RefinementX", 3, "X Refinement" );
    CAF_PDM_InitScriptableField( &m_layersX, "LayersX", 3, "X Layers" );
    CAF_PDM_InitScriptableField( &m_refinementY, "RefinementY", 3, "Y Refinement" );
    CAF_PDM_InitScriptableField( &m_layersY, "LayersY", 3, "Y Layers" );
    CAF_PDM_InitScriptableField( &m_refinementZ, "RefinementZ", 1, "Z Refinement" );
    CAF_PDM_InitScriptableField( &m_layersZ, "LayersZ", 1, "Z Layers" );

    CAF_PDM_InitScriptableField( &m_boxes,
                                 "Boxes",
                                 std::vector<int>(),
                                 "Boxes",
                                 "",
                                 "Nine values per LGR: I, J, K of its first and last host cell (zero-based), then its refinement in I, J, "
                                 "K" );
    m_boxes.uiCapability()->setUiHidden( true );

    CAF_PDM_InitField( &m_updateLgrs, "UpdateLgrs", false, "" );
    caf::PdmUiPushButtonEditor::configureEditorLabelHidden( &m_updateLgrs );
    m_updateLgrs.uiCapability()->setAttribute( caf::PdmUiPushButtonEditor::Keys::BUTTON_TEXT, "Update LGRs" );
    m_updateLgrs.xmlCapability()->disableIO();

    setRegionName( "Well Refinements" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RimWellRefinementRegion::wellName() const
{
    return m_wellName();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::setWellName( const QString& wellName )
{
    m_wellName = wellName;
    setRegionName( "Well Refinements " + wellName );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int RimWellRefinementRegion::timeStep() const
{
    return m_timeStep();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::setTimeStep( int timeStep )
{
    m_timeStep = timeStep;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RimWellRefinementRegion::lgrName() const
{
    return m_lgrName();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::setLgrName( const QString& name )
{
    m_lgrName = name;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
RicWellLgrTools::RingSettings RimWellRefinementRegion::ringSettings() const
{
    RicWellLgrTools::RingSettings settings;
    if ( m_separateXY() )
    {
        settings.refinement = { m_refinementX(), m_refinementY(), m_refinementZ() };
        settings.layers     = { m_layersX(), m_layersY(), m_layersZ() };
    }
    else
    {
        settings.refinement = { m_lateralRefinement(), m_lateralRefinement(), m_refinementZ() };
        settings.layers     = { m_lateralLayers(), m_lateralLayers(), m_layersZ() };
    }
    return settings;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::setRingSettings( const RicWellLgrTools::RingSettings& settings )
{
    m_separateXY        = settings.refinement[0] != settings.refinement[1] || settings.layers[0] != settings.layers[1];
    m_lateralRefinement = settings.refinement[0];
    m_lateralLayers     = settings.layers[0];
    m_refinementX       = settings.refinement[0];
    m_layersX           = settings.layers[0];
    m_refinementY       = settings.refinement[1];
    m_layersY           = settings.layers[1];
    m_refinementZ       = settings.refinement[2];
    m_layersZ           = settings.layers[2];
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::copySettingsFrom( const RimWellRefinementRegion& other )
{
    m_timeStep          = other.m_timeStep();
    m_lgrName           = other.m_lgrName();
    m_lateralRefinement = other.m_lateralRefinement();
    m_lateralLayers     = other.m_lateralLayers();
    m_separateXY        = other.m_separateXY();
    m_refinementX       = other.m_refinementX();
    m_layersX           = other.m_layersX();
    m_refinementY       = other.m_refinementY();
    m_layersY           = other.m_layersY();
    m_refinementZ       = other.m_refinementZ();
    m_layersZ           = other.m_layersZ();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<RicWellLgrTools::RefinedBox> RimWellRefinementRegion::boxes() const
{
    std::vector<RicWellLgrTools::RefinedBox> boxes;
    const auto&                              values = m_boxes();
    for ( size_t n = 0; n + 9 <= values.size(); n += 9 )
    {
        boxes.push_back( { caf::VecIjk0( values[n], values[n + 1], values[n + 2] ),
                           caf::VecIjk0( values[n + 3], values[n + 4], values[n + 5] ),
                           { values[n + 6], values[n + 7], values[n + 8] } } );
    }
    return boxes;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::setBoxes( const std::vector<RicWellLgrTools::RefinedBox>& boxes )
{
    std::vector<int> values;
    for ( const auto& box : boxes )
    {
        for ( size_t d = 0; d < 3; ++d )
            values.push_back( static_cast<int>( box.min[d] ) );
        for ( size_t d = 0; d < 3; ++d )
            values.push_back( static_cast<int>( box.max[d] ) );
        for ( int d = 0; d < 3; ++d )
            values.push_back( box.factors[d] );
    }
    m_boxes = values;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::setSettingsOnly( bool settingsOnly )
{
    m_settingsOnly = settingsOnly;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<RimRefinementRegion::PreviewBox> RimWellRefinementRegion::previewBoxes() const
{
    std::vector<PreviewBox> previews;
    for ( const auto& box : boxes() )
    {
        const cvf::Vec3st size( box.max.i() - box.min.i() + 1, box.max.j() - box.min.j() + 1, box.max.k() - box.min.k() + 1 );
        const cvf::Vec3st factors( box.factors[0], box.factors[1], box.factors[2] );
        previews.push_back( { box.min, box.max, std::make_shared<RigUniformRefinement>( factors, size ) } );
    }
    return previews;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RimWellRefinementRegion::isSectorGridRefinement() const
{
    return false;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::defineUiOrdering( QString uiConfigName, caf::PdmUiOrdering& uiOrdering )
{
    if ( !m_settingsOnly ) uiOrdering.add( userDescriptionField() );
    uiOrdering.add( &m_wellName );
    uiOrdering.add( &m_timeStep );
    uiOrdering.add( &m_lgrName );

    auto* lateral = uiOrdering.addNewGroup( "Lateral Refinement (X, Y)" );
    lateral->add( &m_separateXY );
    if ( m_separateXY() )
    {
        lateral->add( &m_refinementX );
        lateral->add( &m_layersX );
        lateral->add( &m_refinementY );
        lateral->add( &m_layersY );
    }
    else
    {
        lateral->add( &m_lateralRefinement );
        lateral->add( &m_lateralLayers );
    }

    auto* vertical = uiOrdering.addNewGroup( "Vertical Refinement (Z)" );
    vertical->add( &m_refinementZ );
    vertical->add( &m_layersZ );

    const auto  rings = ringSettings();
    QStringList factors;
    for ( int r = 0; r < std::max( rings.layers[0], rings.layers[1] ); ++r )
    {
        const int fx = r < rings.layers[0] ? rings.refinement[0] * ( rings.layers[0] - r ) : 1;
        const int fy = r < rings.layers[1] ? rings.refinement[1] * ( rings.layers[1] - r ) : 1;
        factors.push_back( QString( "%1 x %2" ).arg( fx ).arg( fy ) );
    }
    uiOrdering.addNewLabel( "Path cells and rings: " + factors.join( ", " ) );

    if ( !m_settingsOnly )
    {
        size_t     cellCount = 0;
        const auto current   = boxes();
        for ( const auto& box : current )
        {
            cellCount += ( box.max.i() - box.min.i() + 1 ) * ( box.max.j() - box.min.j() + 1 ) * ( box.max.k() - box.min.k() + 1 ) *
                         box.factors[0] * box.factors[1] * box.factors[2];
        }
        uiOrdering.addNewLabel( QString( "%1 amalgamated LGRs, %2 cells" ).arg( current.size() ).arg( cellCount ) );
        uiOrdering.add( &m_updateLgrs );
        uiOrdering.add( objectToggleField() );
    }

    uiOrdering.skipRemainingFields( true );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void RimWellRefinementRegion::fieldChangedByUi( const caf::PdmFieldHandle* changedField, const QVariant& oldValue, const QVariant& newValue )
{
    if ( changedField == &m_updateLgrs )
    {
        m_updateLgrs = false;
        if ( auto result = RicWellRefinementTools::applyRegion( this ); !result )
        {
            RiaLogging::error( result.error().toStdString() );
            RiuMessageDialog::showError( nullptr, "Well Refinements", result.error() );
        }
        updateConnectedEditors();
        return;
    }

    if ( changedField == &m_separateXY && m_separateXY() )
    {
        m_refinementX = m_lateralRefinement();
        m_refinementY = m_lateralRefinement();
        m_layersX     = m_lateralLayers();
        m_layersY     = m_lateralLayers();
    }

    for ( auto* field :
          { &m_lateralRefinement, &m_lateralLayers, &m_refinementX, &m_layersX, &m_refinementY, &m_layersY, &m_refinementZ, &m_layersZ } )
    {
        if ( changedField == field && ( *field )() < 1 ) *field = 1;
    }

    updateConnectedEditors();
    if ( auto view = firstAncestorOrThisOfType<Rim3dView>() ) view->scheduleCreateDisplayModelAndRedraw();
}
