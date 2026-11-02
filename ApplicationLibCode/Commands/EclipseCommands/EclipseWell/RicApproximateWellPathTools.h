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

#include <QString>

#include <expected>

class RimEclipseCase;
class RimSimWellInView;
class RimWellPath;

//==================================================================================================
///
/// An approximate well path for a simulation well, made from its connections in the simulation results
/// (WELSPECS/COMPDAT level data, segments are ignored). The path runs from above the well head through
/// the centres of the connected cells, in connection order. Each connected cell gets a perforation over the
/// path in the cell, with the diameter and skin of the connection, and a perforated length factor that makes
/// the connection factor of a completion export on the simulation grid equal to the simulated one. On other
/// grids, e.g. with LGRs, the same factor scales the perforated lengths.
///
//==================================================================================================
namespace RicApproximateWellPathTools
{
QString wellPathName( const QString& simulationWellName );

// Time step of the connections, or the last time step with connections if negative. A previous approximate
// well path of the well is replaced. With throughCellFaces (experimental), the path is instead the shortest one
// that passes the connected cells in order through their shared faces, straight inside each cell.
std::expected<RimWellPath*, QString>
    createApproximateWellPath( RimEclipseCase* eclipseCase, const QString& simulationWellName, int timeStep, bool throughCellFaces = false );

// The time step of the well's view if the well has connections then, else -1
int viewTimeStepWithConnections( const RimSimWellInView& simWell );
} // namespace RicApproximateWellPathTools
