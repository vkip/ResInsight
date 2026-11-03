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
/// An approximate well path for a simulation well, made from its connections in the simulation results. The path
/// runs from above the well head through the connected cells in connection order: straight inside each cell, and from
/// one connected cell to the next through their common face. The crossing points on the faces make the path short and
/// smooth. Each connected cell gets a perforation over the path in the cell, with the diameter and skin of the
/// connection, and a perforated length factor that makes the connection factor of a completion export on the
/// simulation grid equal to the simulated one. On other grids, e.g. with LGRs, the same factor scales the perforated
/// lengths.
///
/// For multi-segment wells, the segment data of the restart file (WELSEGS, COMPSEGS, WSEGVALV, WSEGAICD, WSEGSICD) is
/// used too: the measured depths of the path are the lengths along the tubing (from the connection intervals along the
/// tubing where they follow each other, else fitted to the depths of the tubing), and the well path gets the multi-segment
/// settings (top segment as reference measured depth, tubing diameter and roughness, pressure drop). Each perforation of a
/// connection on a valve or device segment gets a valve: an ICD with the flow coefficient and its share of the area of a
/// valve segment, or an AICD/SICD with the parameters of a device segment, and a strength that includes the flow scaling
/// factor of the device and the share of its flow. Each lateral branch with connections becomes a well path lateral
/// ("SIM_<well> L<branch>"), tied in to the path of its parent branch at the length of its outlet along the tubing.
/// Device settings that change with time are not reconstructed, only those at the time step.
///
//==================================================================================================
namespace RicApproximateWellPathTools
{
QString wellPathName( const QString& simulationWellName );

// Time step of the connections, or the last time step with connections if negative. A previous approximate
// well path of the well is replaced.
std::expected<RimWellPath*, QString>
    createApproximateWellPath( RimEclipseCase* eclipseCase, const QString& simulationWellName, int timeStep, bool useMultiSegmentData = true );

// The time step of the well's view if the well has connections then, else -1
int viewTimeStepWithConnections( const RimSimWellInView& simWell );
} // namespace RicApproximateWellPathTools
