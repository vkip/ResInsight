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

#include "RifRestartWellConnectionReader.h"

#include "opm/io/eclipse/ERst.hpp"
#include "opm/output/eclipse/VectorItems/connection.hpp"
#include "opm/output/eclipse/VectorItems/intehead.hpp"
#include "opm/output/eclipse/VectorItems/msw.hpp"
#include "opm/output/eclipse/VectorItems/well.hpp"

#include <QDir>
#include <QFileInfo>
#include <QRegularExpression>

namespace
{
namespace VI = Opm::RestartIO::Helpers::VectorItems;

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString trimmed( const std::string& text )
{
    return QString::fromStdString( text ).trimmed();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<RifRestartWellConnectionReader::Well> wellsOfReportStep( Opm::EclIO::ERst& rst, int step )
{
    using Well       = RifRestartWellConnectionReader::Well;
    using Connection = RifRestartWellConnectionReader::Connection;

    std::vector<Well> wells;

    const auto& inteHead = rst.getRestartData<int>( "INTEHEAD", step );
    const int   nWells   = inteHead[VI::intehead::NWELLS];
    if ( nWells <= 0 || !rst.hasArray( "IWEL", step ) ) return wells;

    const int niwelz = inteHead[VI::intehead::NIWELZ];
    const int nswelz = inteHead[VI::intehead::NSWELZ];
    const int nzwelz = inteHead[VI::intehead::NZWELZ];
    const int niconz = inteHead[VI::intehead::NICONZ];
    const int nsconz = inteHead[VI::intehead::NSCONZ];
    const int ncwmax = inteHead[VI::intehead::NCWMAX];
    const int nzgrpz = inteHead[VI::intehead::NZGRPZ];

    const auto& iwel = rst.getRestartData<int>( "IWEL", step );
    const auto& swel = rst.getRestartData<float>( "SWEL", step );
    const auto& zwel = rst.getRestartData<std::string>( "ZWEL", step );

    const bool  hasConnections = rst.hasArray( "ICON", step ) && rst.hasArray( "SCON", step );
    const auto  noInts         = std::vector<int>{};
    const auto  noFloats       = std::vector<float>{};
    const auto& icon           = hasConnections ? rst.getRestartData<int>( "ICON", step ) : noInts;
    const auto& scon           = hasConnections ? rst.getRestartData<float>( "SCON", step ) : noFloats;

    const bool  hasSegments = rst.hasArray( "ISEG", step ) && rst.hasArray( "RSEG", step );
    const auto& iseg        = hasSegments ? rst.getRestartData<int>( "ISEG", step ) : noInts;
    // RSEG is written in double precision, but accept single precision too
    std::vector<double> rseg;
    if ( hasSegments )
    {
        try
        {
            rseg = rst.getRestartData<double>( "RSEG", step );
        }
        catch ( const std::exception& )
        {
            const auto& values = rst.getRestartData<float>( "RSEG", step );
            rseg.assign( values.begin(), values.end() );
        }
    }
    const int nsegmx = inteHead[VI::intehead::NSEGMX];
    const int nisegz = inteHead[VI::intehead::NISEGZ];
    const int nrsegz = inteHead[VI::intehead::NRSEGZ];

    const auto  noStrings = std::vector<std::string>{};
    const auto& zgrp      = rst.hasArray( "ZGRP", step ) ? rst.getRestartData<std::string>( "ZGRP", step ) : noStrings;

    // Values below this are defaulted
    const double undefined = -1.0e19;

    for ( int w = 0; w < nWells; ++w )
    {
        if ( static_cast<size_t>( w + 1 ) * niwelz > iwel.size() || static_cast<size_t>( w + 1 ) * nswelz > swel.size() ) break;
        const int* iw = iwel.data() + static_cast<size_t>( w ) * niwelz;

        Well well;
        well.name           = trimmed( zwel[static_cast<size_t>( w ) * nzwelz] );
        well.headI          = iw[VI::IWell::IHead] - 1;
        well.headJ          = iw[VI::IWell::JHead] - 1;
        well.wellType       = iw[VI::IWell::WType];
        well.preferredPhase = iw[VI::IWell::PreferredPhase];
        well.allowCrossFlow = iw[VI::IWell::XFlow] != 0;
        well.isMultiSegment = iw[VI::IWell::MsWID] > 0;
        well.isInLgr        = niwelz > static_cast<int>( VI::IWell::LGRIndex ) && iw[VI::IWell::LGRIndex] > 0;

        const int groupIndex = iw[VI::IWell::Group] - 1;
        if ( groupIndex >= 0 && nzgrpz > 0 && static_cast<size_t>( groupIndex ) * static_cast<size_t>( nzgrpz ) < zgrp.size() )
        {
            well.groupName = trimmed( zgrp[static_cast<size_t>( groupIndex ) * nzgrpz] );
        }

        const float* sw = swel.data() + static_cast<size_t>( w ) * nswelz;
        if ( nswelz > static_cast<int>( VI::SWell::DatumDepth ) && sw[VI::SWell::DatumDepth] > undefined )
            well.referenceDepth = sw[VI::SWell::DatumDepth];
        if ( nswelz > static_cast<int>( VI::SWell::DrainageRadius ) && sw[VI::SWell::DrainageRadius] > undefined )
        {
            well.drainageRadius = sw[VI::SWell::DrainageRadius];
        }

        const int connectionCount = hasConnections ? iw[VI::IWell::NConn] : 0;
        for ( int c = 0; c < connectionCount; ++c )
        {
            const size_t iOffset = ( static_cast<size_t>( w ) * ncwmax + c ) * niconz;
            const size_t sOffset = ( static_cast<size_t>( w ) * ncwmax + c ) * nsconz;
            if ( iOffset + static_cast<size_t>( niconz ) > icon.size() || sOffset + static_cast<size_t>( nsconz ) > scon.size() ) break;

            const int*   ic = icon.data() + iOffset;
            const float* sc = scon.data() + sOffset;

            Connection connection;
            connection.i                = ic[VI::IConn::CellI] - 1;
            connection.j                = ic[VI::IConn::CellJ] - 1;
            connection.k                = ic[VI::IConn::CellK] - 1;
            connection.isOpen           = ic[VI::IConn::ConnStat] > 0;
            connection.direction        = ic[VI::IConn::ConnDir];
            connection.completionNumber = ic[VI::IConn::ComplNum];

            // The connection transmissibility factor without multipliers (WPIMULT, ROCKTAB), if written
            const double connTrans      = nsconz > static_cast<int>( VI::SConn::ConnTrans ) ? sc[VI::SConn::ConnTrans] : 0.0;
            connection.connectionFactor = connTrans > 0.0 ? connTrans : sc[VI::SConn::EffConnTrans];
            connection.depth            = sc[VI::SConn::Depth];
            connection.diameter         = sc[VI::SConn::Diameter];
            connection.kh               = sc[VI::SConn::EffectiveKH];
            connection.skinFactor       = sc[VI::SConn::SkinFactor];
            connection.segment          = niconz > static_cast<int>( VI::IConn::Segment ) ? ic[VI::IConn::Segment] : 0;
            if ( connection.segment > 0 && nsconz > static_cast<int>( VI::SConn::SegDistStart ) )
            {
                connection.segmentStart = sc[VI::SConn::SegDistStart];
                connection.segmentEnd   = sc[VI::SConn::SegDistEnd];
            }

            well.connections.push_back( connection );
        }

        const int mswIndex = iw[VI::IWell::MsWID] - 1;
        if ( mswIndex >= 0 && hasSegments && nsegmx > 0 )
        {
            well.pressureDropModel = niwelz > static_cast<int>( VI::IWell::MSW_PlossMod ) ? iw[VI::IWell::MSW_PlossMod] : 0;
            const int segmentCount = std::min( iw[VI::IWell::NWseg], nsegmx );
            for ( int slot = 0; slot < segmentCount; ++slot )
            {
                const size_t iOffset = ( static_cast<size_t>( mswIndex ) * nsegmx + slot ) * nisegz;
                const size_t rOffset = ( static_cast<size_t>( mswIndex ) * nsegmx + slot ) * nrsegz;
                if ( iOffset + static_cast<size_t>( nisegz ) > iseg.size() || rOffset + static_cast<size_t>( nrsegz ) > rseg.size() ) break;

                const int*    is = iseg.data() + iOffset;
                const double* rs = rseg.data() + rOffset;

                using Segment = RifRestartWellConnectionReader::Segment;
                Segment segment;
                segment.number    = slot + 1;
                segment.outlet    = is[VI::ISeg::OutSeg];
                segment.branch    = is[VI::ISeg::BranchNo];
                segment.length    = rs[VI::RSeg::DistBHPRef];
                segment.depth     = rs[VI::RSeg::SegNodeDepth];
                segment.diameter  = rs[VI::RSeg::SegDiam];
                segment.roughness = rs[VI::RSeg::SegRough];
                segment.area      = rs[VI::RSeg::SegArea];
                segment.volume    = rs[VI::RSeg::SegVolume];
                switch ( is[VI::ISeg::SegmentType] )
                {
                    case VI::ISeg::Value::Regular:
                        segment.type = Segment::Type::REGULAR;
                        break;
                    case VI::ISeg::Value::Valve:
                        segment.type = Segment::Type::VALVE;
                        break;
                    case VI::ISeg::Value::SpiralICD:
                        segment.type = Segment::Type::SPIRAL_ICD;
                        break;
                    case VI::ISeg::Value::AutoICD:
                        segment.type = Segment::Type::AUTO_ICD;
                        break;
                    default:
                        segment.type = Segment::Type::OTHER;
                }
                if ( segment.type == Segment::Type::VALVE )
                {
                    segment.valveFlowCoefficient = rs[VI::RSeg::ValveFlowCoeff];
                    segment.valveArea            = rs[VI::RSeg::ValveArea];
                    segment.valveMaxArea         = rs[VI::RSeg::ValveMaxArea];
                    segment.valveLength          = rs[VI::RSeg::ValveLength];
                }
                else if ( ( segment.type == Segment::Type::SPIRAL_ICD || segment.type == Segment::Type::AUTO_ICD ) &&
                          nisegz > VI::ISeg::ICDOpenShutFlag && nrsegz > VI::RSeg::flowFractionGasViscosityExponent )
                {
                    segment.deviceStrength        = rs[VI::RSeg::DeviceBaseStrength];
                    segment.deviceScalingMethod   = is[VI::ISeg::ICDScalingMode];
                    segment.deviceOpen            = is[VI::ISeg::ICDOpenShutFlag] == 0;
                    segment.deviceLength          = rs[VI::RSeg::ICDLength];
                    segment.calibrationDensity    = rs[VI::RSeg::CalibrFluidDensity];
                    segment.calibrationViscosity  = rs[VI::RSeg::CalibrFluidViscosity];
                    segment.criticalWaterFraction = rs[VI::RSeg::CriticalWaterFraction];
                    segment.transitionWidth       = rs[VI::RSeg::TransitionRegWidth];
                    segment.maxEmulsionRatio      = rs[VI::RSeg::MaxEmulsionRatio];
                    segment.maxValidFlowRate      = rs[VI::RSeg::MaxValidFlowRate];

                    // The scaling factor is written as a length when it is one (method 1, or a negative length with a
                    // negative method). The simulator uses the length in SI units as the factor.
                    segment.deviceScalingFactor = rs[VI::RSeg::ScalingFactor];
                    const bool scalingIsLength  = segment.deviceScalingMethod == 1 ||
                                                 ( segment.deviceScalingMethod < 0 && segment.deviceLength < 0.0 );
                    if ( scalingIsLength && inteHead[VI::intehead::UNIT] == 2 ) segment.deviceScalingFactor *= 0.3048;
                    if ( scalingIsLength && inteHead[VI::intehead::UNIT] == 3 ) segment.deviceScalingFactor *= 0.01;

                    if ( segment.type == Segment::Type::AUTO_ICD )
                    {
                        segment.flowRateExponent           = rs[VI::RSeg::FlowRateExponent];
                        segment.viscosityExponent          = rs[VI::RSeg::ViscFuncExponent];
                        segment.densityFractionExponents   = { rs[VI::RSeg::flowFractionOilDensityExponent],
                                                               rs[VI::RSeg::flowFractionWaterDensityExponent],
                                                               rs[VI::RSeg::flowFractionGasDensityExponent] };
                        segment.viscosityFractionExponents = { rs[VI::RSeg::flowFractionOilViscosityExponent],
                                                               rs[VI::RSeg::flowFractionWaterViscosityExponent],
                                                               rs[VI::RSeg::flowFractionGasViscosityExponent] };
                    }
                }
                well.segments.push_back( segment );
            }
        }

        wells.push_back( well );
    }

    return wells;
}
} // namespace

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QStringList RifRestartWellConnectionReader::restartFiles( const QString& gridFileName )
{
    const QFileInfo gridInfo( gridFileName );
    const QDir      folder   = gridInfo.absoluteDir();
    const QString   baseName = gridInfo.completeBaseName();

    for ( const QString& extension : { "UNRST", "FUNRST" } )
    {
        for ( const QString& name : { baseName + "." + extension, baseName + "." + extension.toLower() } )
        {
            if ( folder.exists( name ) ) return { folder.absoluteFilePath( name ) };
        }
    }

    QStringList              files;
    const QRegularExpression separate( "^" + QRegularExpression::escape( baseName ) + "\\.[XFxf]\\d{4}$" );
    const QFileInfoList      entries = folder.entryInfoList( QDir::Files, QDir::Name );
    for ( const QFileInfo& entry : entries )
    {
        if ( separate.match( entry.fileName() ).hasMatch() ) files.push_back( entry.absoluteFilePath() );
    }
    return files;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<std::vector<RifRestartWellConnectionReader::Well>, QString>
    RifRestartWellConnectionReader::readWells( const QStringList& restartFiles, const QDate& date )
{
    if ( restartFiles.isEmpty() ) return std::unexpected( "No restart file found." );

    try
    {
        for ( const QString& file : restartFiles )
        {
            Opm::EclIO::ERst rst( file.toStdString() );
            for ( int step : rst.listOfReportStepNumbers() )
            {
                rst.loadReportStepNumber( step );
                const auto& inteHead = rst.getRestartData<int>( "INTEHEAD", step );
                const QDate stepDate( inteHead[VI::intehead::YEAR], inteHead[VI::intehead::MONTH], inteHead[VI::intehead::DAY] );
                if ( stepDate == date ) return wellsOfReportStep( rst, step );
            }
        }
    }
    catch ( const std::exception& e )
    {
        return std::unexpected( QString( "Failed to read the restart data: %1" ).arg( e.what() ) );
    }

    return std::unexpected( QString( "The restart files have no report step at %1." ).arg( date.toString( Qt::ISODate ) ) );
}
