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

#include "RigRefinement.h"

#include "cafAppEnum.h"
#include "cafVecIjk.h"

#include <QString>

#include <array>
#include <expected>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

class LgrInfo;
class RigEclipseCaseData;
class RigLocalGrid;
class RigMainGrid;
class RimEclipseCase;

//==================================================================================================
///
/// Local grid refinements around simulation wells, i.e. wells known from the simulation results only.
///
//==================================================================================================
namespace RicWellLgrTools
{
struct IjkLess
{
    bool operator()( const caf::VecIjk0& a, const caf::VecIjk0& b ) const
    {
        return std::tie( a.k(), a.j(), a.i() ) < std::tie( b.k(), b.j(), b.i() );
    }
};
using CellSet = std::set<caf::VecIjk0, IjkLess>;

enum class LateralRefinement
{
    UNIFORM,
    LOGARITHMIC
};
using LateralRefinementEnum = caf::AppEnum<LateralRefinement>;

struct Settings
{
    QString           lgrName; // Empty: derived from the well name
    int               timeStep          = -1; // Connections of this time step, or of all time steps if negative
    LateralRefinement lateralRefinement = LateralRefinement::LOGARITHMIC;
    int               refinementI       = 3; // Uniform lateral refinement
    int               refinementJ       = 3;
    int               refinementK       = 2; // Vertical refinement, always uniform
    double            minCellSize       = 1.0; // Logarithmic: lateral cell size at the well [m]
    double            growthFactor      = 1.5; // Logarithmic: size ratio of neighbouring cells away from the well
    int               paddingIJ         = 2; // Main grid cells around the connections, laterally
    int               paddingK          = 0; // and vertically
};

// Main grid cells (0-based) connected to the simulation well. Connections in LGRs map to their host cells.
CellSet simulationWellCells( const RigEclipseCaseData& caseData, const QString& wellName, int timeStep = -1 );

// Cumulative fractions per parent cell for a graded subdivision of one direction. The cell size grows
// geometrically with the distance from the centres of the well cells, h = max(minCellSize, (growthFactor - 1) * d),
// and each parent cell holds at least one refined cell.
std::vector<std::vector<double>> gradedCumulativeFractions( const std::vector<double>& parentCellWidths,
                                                            const std::set<size_t>&    wellParentCells,
                                                            double                     minCellSize,
                                                            double                     growthFactor );

QString defaultLgrName( const QString& wellName );

// Name of LGR number index (zero-based) of count LGRs around one well, at most eight characters
QString lgrName( const QString& baseName, size_t index, size_t count );

using CellBox = std::pair<caf::VecIjk0, caf::VecIjk0>;

// Boxes around the well cells, padded and clipped to the grid, merging boxes that overlap or touch, or whose
// joined box holds at most 1.5 times their cells. A well with connections in separate parts of the grid, e.g.
// in two zones, gets one box per part.
std::vector<CellBox>
    wellCellBoxes( const CellSet& wellCells, const std::array<size_t, 3>& padding, const std::array<size_t, 3>& gridCellCount );

// One LGR per box of wellCellBoxes()
std::expected<std::vector<LgrInfo>, QString>
    buildLgrsForSimulationWell( RimEclipseCase& eclipseCase, const QString& wellName, const Settings& settings );

// Replaces any temporary LGRs of the well, keeping the temporary LGRs of other wells
std::expected<std::vector<LgrInfo>, QString>
    createTemporaryLgrsForSimulationWell( RimEclipseCase* eclipseCase, const QString& wellName, const Settings& settings );

std::vector<const RigLocalGrid*> temporaryLgrsForWell( const RigMainGrid& mainGrid, const QString& wellName );
LgrInfo                          lgrInfoFromTemporaryLgr( const RigMainGrid& mainGrid, const RigLocalGrid& localGrid );

// Number of simulation wells with connections in the box, at least one
int wellCountInBox( const RigEclipseCaseData& caseData, const caf::VecIjk0& boxMin, const caf::VecIjk0& boxMax );

// Subdivision of the LGR box, uniform or graded, relative to the box
std::unique_ptr<RigRefinement> lgrRefinement( const LgrInfo& lgrInfo );

// NXFIN/HXFIN, NYFIN/HYFIN and NZFIN/HZFIN for the directions that are not uniformly refined
QString gradingKeywords( const RigRefinement& refinement );

// CARFIN ... ENDFIN, including the grading keywords
QString carfinBlock( const LgrInfo& lgrInfo, int maxWellCount );
} // namespace RicWellLgrTools
