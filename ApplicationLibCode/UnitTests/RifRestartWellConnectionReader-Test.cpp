#include "gtest/gtest.h"

#include "RiaTestDataDirectory.h"

#include "RifRestartWellConnectionReader.h"

#include <QDate>

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
