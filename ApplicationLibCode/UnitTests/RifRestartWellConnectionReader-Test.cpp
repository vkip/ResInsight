#include "gtest/gtest.h"

#include "RiaTestDataDirectory.h"

#include "RifRestartWellConnectionReader.h"

#include <QDate>

#include <map>

//--------------------------------------------------------------------------------------------------
/// A horizontal multi-segment well (TestModels/RifRestartWellConnectionReader/MSW_ICD.DATA): tubing segments 1-4 on
/// branch 1, and two valve segments (WSEGVALV) on branches 2 and 3 that the two last connections drain through
//--------------------------------------------------------------------------------------------------
TEST( RifRestartWellConnectionReaderTest, MultiSegmentWellWithValves )
{
    const QString gridFile = QString( "%1/RifRestartWellConnectionReader/MSW_ICD.EGRID" ).arg( TEST_DATA_DIR );
    const auto    files    = RifRestartWellConnectionReader::restartFiles( gridFile );
    ASSERT_FALSE( files.empty() );

    const auto wells = RifRestartWellConnectionReader::readWells( files, QDate( 1983, 1, 2 ) );
    ASSERT_TRUE( wells.has_value() ) << wells.error().toStdString();
    ASSERT_EQ( 1u, wells->size() );

    const auto& well = wells->front();
    EXPECT_EQ( "P1", well.name.toStdString() );
    EXPECT_TRUE( well.isMultiSegment );
    EXPECT_EQ( 1, well.pressureDropModel ); // HF-

    using Type = RifRestartWellConnectionReader::Segment::Type;
    ASSERT_EQ( 6u, well.segments.size() );
    EXPECT_EQ( 0, well.segments[0].outlet );
    EXPECT_DOUBLE_EQ( 2050.0, well.segments[0].length );
    EXPECT_DOUBLE_EQ( 2002.5, well.segments[0].depth );
    EXPECT_EQ( 1, well.segments[3].branch );
    EXPECT_DOUBLE_EQ( 2350.0, well.segments[3].length );
    EXPECT_NEAR( 0.15, well.segments[3].diameter, 1e-9 );
    EXPECT_NEAR( 0.0001, well.segments[3].roughness, 1e-9 );
    EXPECT_EQ( Type::REGULAR, well.segments[3].type );

    const auto& valve = well.segments[5];
    EXPECT_EQ( 6, valve.number );
    EXPECT_EQ( 4, valve.outlet );
    EXPECT_EQ( 3, valve.branch );
    EXPECT_EQ( Type::VALVE, valve.type );
    EXPECT_NEAR( 0.8, valve.valveFlowCoefficient, 1e-6 );
    EXPECT_NEAR( 3.0e-4, valve.valveArea, 1e-9 );

    ASSERT_EQ( 3u, well.connections.size() );
    EXPECT_EQ( 2, well.connections[0].segment );
    EXPECT_EQ( 5, well.connections[1].segment );
    EXPECT_EQ( 6, well.connections[2].segment );
    EXPECT_NEAR( 2200.0, well.connections[1].segmentStart, 1e-3 );
    EXPECT_NEAR( 2300.0, well.connections[1].segmentEnd, 1e-3 );
}

//--------------------------------------------------------------------------------------------------
/// A multi-segment well with a lateral (TestModels/RifRestartWellConnectionReader/MSW_DEVICES.DATA): tubing segments 1-4
/// on the main branch, 5-7 on lateral branch 2 from segment 3, SICD segments 8 and 11, and AICD segments 9 and 10
//--------------------------------------------------------------------------------------------------
TEST( RifRestartWellConnectionReaderTest, MultiSegmentWellWithLateralAndDevices )
{
    const QString gridFile = QString( "%1/RifRestartWellConnectionReader/MSW_DEVICES.EGRID" ).arg( TEST_DATA_DIR );
    const auto    wells =
        RifRestartWellConnectionReader::readWells( RifRestartWellConnectionReader::restartFiles( gridFile ), QDate( 1983, 1, 2 ) );
    ASSERT_TRUE( wells.has_value() ) << wells.error().toStdString();
    ASSERT_EQ( 1u, wells->size() );

    const auto& well = wells->front();
    using Type       = RifRestartWellConnectionReader::Segment::Type;
    ASSERT_EQ( 11u, well.segments.size() );

    const auto& lateralStart = well.segments[4];
    EXPECT_EQ( 2, lateralStart.branch );
    EXPECT_EQ( 3, lateralStart.outlet );
    EXPECT_NEAR( 0.12, lateralStart.diameter, 1e-9 );
    EXPECT_NEAR( 0.0002, lateralStart.roughness, 1e-9 );

    // Method 1: the scaling factor is the length
    const auto& sicd = well.segments[7];
    EXPECT_EQ( Type::SPIRAL_ICD, sicd.type );
    EXPECT_NEAR( 2.0e-5, sicd.deviceStrength, 1e-12 );
    EXPECT_EQ( 1, sicd.deviceScalingMethod );
    EXPECT_NEAR( 1.5, sicd.deviceScalingFactor, 1e-9 );
    EXPECT_NEAR( 1000.0, sicd.calibrationDensity, 1e-6 );
    EXPECT_NEAR( 1.0, sicd.calibrationViscosity, 1e-9 );
    EXPECT_NEAR( 0.5, sicd.criticalWaterFraction, 1e-9 );
    EXPECT_NEAR( 0.05, sicd.transitionWidth, 1e-9 );
    EXPECT_NEAR( 5.0, sicd.maxEmulsionRatio, 1e-9 );
    EXPECT_TRUE( sicd.deviceOpen );

    // Default method: the scaling factor is the device length over the length of the outlet segment
    const auto& lateralSicd = well.segments[10];
    EXPECT_EQ( Type::SPIRAL_ICD, lateralSicd.type );
    EXPECT_EQ( 7, lateralSicd.outlet );
    EXPECT_NEAR( 0.5, lateralSicd.deviceScalingFactor, 1e-9 );

    const auto& aicd = well.segments[9];
    EXPECT_EQ( Type::AUTO_ICD, aicd.type );
    EXPECT_EQ( 5, aicd.outlet );
    EXPECT_NEAR( 2.0e-6, aicd.deviceStrength, 1e-12 );
    EXPECT_NEAR( 0.12, aicd.deviceScalingFactor, 1e-9 );
    EXPECT_NEAR( 2.5, aicd.flowRateExponent, 1e-9 );
    EXPECT_NEAR( 0.7, aicd.viscosityExponent, 1e-9 );
    EXPECT_NEAR( 500.0, aicd.maxValidFlowRate, 1e-6 );
    EXPECT_NEAR( 1.1, aicd.densityFractionExponents[0], 1e-9 );
    EXPECT_NEAR( 1.0, aicd.viscosityFractionExponents[2], 1e-9 );

    ASSERT_EQ( 6u, well.connections.size() );
    std::map<int, int> segmentOfCell; // by j * 10 + i, one-based
    for ( const auto& connection : well.connections )
    {
        segmentOfCell[( connection.j + 1 ) * 10 + connection.i + 1] = connection.segment;
    }
    EXPECT_EQ( 2, segmentOfCell[12] );
    EXPECT_EQ( 8, segmentOfCell[13] );
    EXPECT_EQ( 9, segmentOfCell[14] );
    EXPECT_EQ( 10, segmentOfCell[23] );
    EXPECT_EQ( 6, segmentOfCell[33] );
    EXPECT_EQ( 11, segmentOfCell[43] );
}
