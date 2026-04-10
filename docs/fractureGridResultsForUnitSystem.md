# `fractureGridResultsForUnitSystem` — How It Is Used

## Overview

`fractureGridResultsForUnitSystem` is a pure virtual method declared in `RimMeshFractureTemplate` and implemented by both `RimStimPlanFractureTemplate` and `RimThermalFractureTemplate`. It returns per-cell result values from the fracture grid, converting them to the unit system of the fracture template (metric or field).

## Declaration

`RimMeshFractureTemplate.h:85`

```cpp
virtual std::vector<double> fractureGridResultsForUnitSystem(
    const QString&                resultName,
    const QString&                unitName,
    size_t                        timeStepIndex,
    RiaDefines::EclipseUnitSystem requiredUnitSystem ) const = 0;
```

| Parameter | Description |
|---|---|
| `resultName` | Name of the fracture property (e.g. width, conductivity) |
| `unitName` | Unit string as stored in the fracture definition file (e.g. `"m"`, `"ft"`, `"/m"`) |
| `timeStepIndex` | Index into the fracture time step data |
| `requiredUnitSystem` | Target unit system (currently always passed as `fractureTemplateUnit()`) |

**Return value:** one `double` per fracture grid cell, converted to the template's unit system.

---

## Implementation: `RimStimPlanFractureTemplate`

`RimStimPlanFractureTemplate.cpp:287`

```cpp
std::vector<double> RimStimPlanFractureTemplate::fractureGridResultsForUnitSystem(
    const QString& resultName, const QString& unitName,
    size_t timeStepIndex, RiaDefines::EclipseUnitSystem requiredUnitSystem ) const
{
    auto resultValues = m_stimPlanFractureDefinitionData->fractureGridResults(
        resultName, unitName, m_activeTimeStepIndex );

    if ( fractureTemplateUnit() == UNITS_METRIC )
        for ( auto& v : resultValues )
            v = RiaEclipseUnitTools::convertToMeter( v, unitName );
    else if ( fractureTemplateUnit() == UNITS_FIELD )
        for ( auto& v : resultValues )
            v = RiaEclipseUnitTools::convertToFeet( v, unitName );

    return resultValues;
}
```

1. Fetches raw cell values from `RigStimPlanFractureDefinition::fractureGridResults`.
2. Converts each value using `RiaEclipseUnitTools::convertToMeter` or `convertToFeet` based on the template's unit system.
3. The `requiredUnitSystem` parameter is accepted but the actual conversion uses `fractureTemplateUnit()` — the two are expected to match at all call sites.

`RimThermalFractureTemplate` follows the same pattern, but sources data from `RigThermalFractureDefinition` instead.

---

## Call Sites in `RimMeshFractureTemplate`

All call sites are in `RimMeshFractureTemplate.cpp`. The method is never called from outside the template hierarchy.

### 1. `wellFractureIntersectionData` — along-well-path orientation

`RimMeshFractureTemplate.cpp:239–247`

Width and conductivity values are fetched for every fracture grid cell, then a length-weighted mean is computed over the cells intersected by the well path:

```cpp
auto nameUnit = widthParameterNameAndUnit();
widthResultValues = fractureGridResultsForUnitSystem(
    nameUnit.first, nameUnit.second, m_activeTimeStepIndex, fractureTemplateUnit() );

auto nameUnit2 = conductivityParameterNameAndUnit();
conductivityResultValues = fractureGridResultsForUnitSystem(
    nameUnit2.first, nameUnit2.second, m_activeTimeStepIndex, fractureTemplateUnit() );
```

The intersection lengths come from `RigWellPathStimplanIntersector`. The weighted means are used to populate `WellFractureIntersectionData::m_width`, `m_conductivity`, and `m_permeability`.

### 2. `wellFractureIntersectionData` — transverse/default orientation

`RimMeshFractureTemplate.cpp:335`

Only width is needed here (conductivity is read directly from the `RigFractureCell`). The value at the well-center cell index is extracted:

```cpp
auto resultValues = fractureGridResultsForUnitSystem(
    nameUnit.first, nameUnit.second, m_activeTimeStepIndex, fractureTemplateUnit() );

if ( wellCellIndex < resultValues.size() )
    widthInRequiredUnit = resultValues[wellCellIndex];
```

Width and permeability are then set on the intersection data struct.

### 3. `widthResultValues`

`RimMeshFractureTemplate.cpp:437`

A thin wrapper that fetches the full width array for the active time step:

```cpp
std::vector<double> RimMeshFractureTemplate::widthResultValues() const
{
    auto nameUnit = widthParameterNameAndUnit();
    return fractureGridResultsForUnitSystem(
        nameUnit.first, nameUnit.second, m_activeTimeStepIndex, fractureTemplateUnit() );
}
```

This is used by callers that need the complete per-cell width grid (e.g. for visualization).

---

## Data Flow

```
Fracture definition file (.xml / thermal)
        |
        v
RigStimPlanFractureDefinition::fractureGridResults()   (raw values, native file units)
        |
        v
RimStimPlanFractureTemplate::fractureGridResultsForUnitSystem()   (unit-converted)
        |
        v
RimMeshFractureTemplate::wellFractureIntersectionData()
    ├── along-well-path: length-weighted mean → WellFractureIntersectionData
    └── transverse:      single cell lookup   → WellFractureIntersectionData
        |
        v
RimMeshFractureTemplate::widthResultValues()   (full grid, for visualization)
```

---

## Key Files

| File | Role |
|---|---|
| `ApplicationLibCode/ProjectDataModel/Completions/RimMeshFractureTemplate.h` | Pure virtual declaration |
| `ApplicationLibCode/ProjectDataModel/Completions/RimMeshFractureTemplate.cpp` | All call sites |
| `ApplicationLibCode/ProjectDataModel/Completions/RimStimPlanFractureTemplate.cpp` | StimPlan implementation |
| `ApplicationLibCode/ProjectDataModel/Completions/RimThermalFractureTemplate.cpp` | Thermal implementation |
| `ApplicationLibCode/Application/Tools/RiaEclipseUnitTools.h` | Unit conversion helpers |
