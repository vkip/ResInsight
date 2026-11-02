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

#include <QDate>
#include <QString>
#include <QStringList>

#include <expected>
#include <optional>
#include <vector>

//==================================================================================================
///
/// The well and connection data of a restart file at one report step, i.e. what WELSPECS and COMPDAT
/// defined, from the arrays IWEL, SWEL, ZWEL, ICON, SCON and ZGRP of the global grid.
///
//==================================================================================================
namespace RifRestartWellConnectionReader
{
struct Connection
{
    int    i                = 0; // Zero-based
    int    j                = 0;
    int    k                = 0;
    bool   isOpen           = true;
    int    direction        = 3; // 1: X, 2: Y, 3: Z
    double connectionFactor = 0.0;
    double depth            = 0.0;
    double diameter         = 0.0;
    double kh               = 0.0;
    double skinFactor       = 0.0;
    int    completionNumber = 0;
};

struct Well
{
    QString name;
    QString groupName;
    int     headI          = 0; // Zero-based
    int     headJ          = 0;
    int     wellType       = 0; // 1: producer, 2: oil, 3: water, 4: gas injector
    int     preferredPhase = 0; // 1: oil, 2: water, 3: gas, 4: liquid
    bool    allowCrossFlow = true;
    bool    isMultiSegment = false;
    bool    isInLgr        = false;

    std::optional<double> referenceDepth;
    std::optional<double> drainageRadius;

    std::vector<Connection> connections;
};

// Restart files next to the grid file: the unified restart file, or else the separate ones
QStringList restartFiles( const QString& gridFileName );

// Wells of the report step at the given date
std::expected<std::vector<Well>, QString> readWells( const QStringList& restartFiles, const QDate& date );
} // namespace RifRestartWellConnectionReader
