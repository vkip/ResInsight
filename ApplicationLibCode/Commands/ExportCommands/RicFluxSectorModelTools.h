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

#include "cvfVector3.h"

#include <QString>
#include <QStringList>

#include <expected>
#include <vector>

class LgrInfo;
class RifOpmFlowDeckFile;
class RimEclipseCase;
class RimEclipseView;
class RimPolygonFilter;

//==================================================================================================
///
/// A sector model around a simulation well for an OPM Flow run with USEFLUX. The sector is a FLUXNUM
/// region defined by a polygon around the well, its boundary data is built by make_flux from the
/// results of the case, and the parent deck is copied with USEFLUX, FLUXREG and the FLUXNUM region
/// added. Temporary LGRs around the well are carried into the sector deck.
///
//==================================================================================================
namespace RicFluxSectorModelTools
{
struct Settings
{
    QString wellName;
    QString exportFolder; // Empty: 'sector_<well>' next to the grid file
    QString deckFile; // Empty: the DATA file next to the grid file, else the user is asked
    QString makeFluxExecutable; // Empty: make_flux from PATH
    int     marginCells        = 3; // Main grid cells around the well connections, laterally
    int     wellPaddingCells   = 2; // Cells around a well that would otherwise be cut by the sector boundary
    int     lgrPaddingCells    = 2; // Cells between an LGR and the sector boundary
    bool    runMakeFlux        = true;
    int     maxExpansionPasses = 50;
};

struct Result
{
    QString     deckFile;
    QString     fluxFile;
    QString     fluxnumFile;
    QStringList lgrNames;
    int         sectorCellCount       = 0;
    int         activeSectorCellCount = 0;
    QStringList wellsInSector;
    QStringList warnings;
};

// Cells of the main grid in the sector, one flag per main grid cell, and the polygon around them
struct SectorRegion
{
    std::vector<bool>       cells;
    std::vector<cvf::Vec3d> polygon;
    QStringList             wellsInSector;
    QStringList             expandedForWells;
};

std::expected<Result, QString> createSectorModel( RimEclipseView* view, const Settings& settings );

std::expected<SectorRegion, QString> computeSectorRegion( RimEclipseCase& eclipseCase, const Settings& settings );

// Convex hull in the XY plane, counter-clockwise, at the given z
std::vector<cvf::Vec3d> convexHullXY( const std::vector<cvf::Vec3d>& points, double z );

// Moves the deck's own FLUXNUM out of the way of the sector region. Region multipliers (MULTREGT) on FLUXNUM are
// moved to MULTNUM or OPERNUM. Returns a description of what was done.
std::expected<QString, QString> relocateDeckFluxnum( RifOpmFlowDeckFile& deckFile );

// Whether the deck defines the well as a multi-segment well (WELSEGS)
bool isMultiSegmentWell( const RifOpmFlowDeckFile& deckFile, const QString& wellName );

// Converts the well's WELSPECS and COMPDAT to WELSPECL and COMPDATL in the LGR
std::expected<int, QString>
    convertWellToLgr( RifOpmFlowDeckFile& deckFile, RimEclipseCase* eclipseCase, const QString& wellName, const LgrInfo& lgrInfo );

std::expected<void, QString> addLgrDimensions( RifOpmFlowDeckFile& deckFile, const std::vector<LgrInfo>& lgrInfos );
} // namespace RicFluxSectorModelTools
