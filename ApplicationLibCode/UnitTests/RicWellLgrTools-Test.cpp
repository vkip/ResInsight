#include "gtest/gtest.h"

#include "ExportCommands/RicExportLgrFeature.h"
#include "ExportCommands/RicFluxSectorModelTools.h"
#include "ExportCommands/RicWellLgrTools.h"

#include "RifOpmFlowDeckFile.h"

#include "RigNonUniformRefinement.h"

#include <QFile>
#include <QTemporaryDir>
#include <QTextStream>

#include <cmath>

namespace
{
//--------------------------------------------------------------------------------------------------
/// Widths of the refined cells of a parent cell, from its cumulative fractions
//--------------------------------------------------------------------------------------------------
std::vector<double> refinedWidths( const std::vector<double>& fractions, double parentWidth )
{
    std::vector<double> widths;
    double              previous = 0.0;
    for ( double fraction : fractions )
    {
        widths.push_back( ( fraction - previous ) * parentWidth );
        previous = fraction;
    }
    return widths;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString writeDeck( const QTemporaryDir& dir, const QString& gridKeywords )
{
    const QString fileName = dir.path() + "/TEST.DATA";
    QFile         file( fileName );
    if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) ) return {};

    QTextStream( &file ) << "RUNSPEC\nDIMENS\n 2 1 1 /\nOIL\nWATER\nREGDIMS\n 4 4 /\n"
                         << "GRID\nDX\n 2*100 /\nDY\n 2*100 /\nDZ\n 2*10 /\nTOPS\n 2*2000 /\nPORO\n 2*0.2 /\n"
                         << gridKeywords << "PROPS\nSCHEDULE\nTSTEP\n 1 /\nEND\n";
    return fileName;
}
} // namespace

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicWellLgrToolsTest, GradedFractionsAroundVerticalWell )
{
    const std::vector<double> widths    = { 100.0, 100.0, 100.0 };
    const auto                fractions = RicWellLgrTools::gradedCumulativeFractions( widths, { 1 }, 1.0, 1.5 );
    ASSERT_EQ( fractions.size(), 3u );

    for ( const auto& parent : fractions )
    {
        ASSERT_FALSE( parent.empty() );
        EXPECT_DOUBLE_EQ( parent.back(), 1.0 );
        for ( size_t n = 1; n < parent.size(); ++n )
        {
            EXPECT_LT( parent[n - 1], parent[n] );
        }
    }

    // About 1 m at the well in the middle of the well cell, growing away from it, symmetrically
    const auto wellCell = refinedWidths( fractions[1], 100.0 );
    ASSERT_GE( wellCell.size(), 8u );
    const size_t centre = wellCell.size() / 2;
    EXPECT_NEAR( wellCell[centre], 1.0, 0.1 );
    EXPECT_NEAR( wellCell[centre - 1], 1.0, 0.1 );
    for ( size_t n = 0; n + 1 < centre; ++n )
    {
        EXPECT_GE( wellCell[n], wellCell[n + 1] );
        EXPECT_NEAR( wellCell[n], wellCell[wellCell.size() - 1 - n], 1.0e-6 );
    }
    EXPECT_GT( wellCell.front(), 10.0 );

    // The neighbours hold fewer, larger cells
    EXPECT_LT( fractions[0].size(), fractions[1].size() );
    EXPECT_EQ( fractions[0].size(), fractions[2].size() );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicWellLgrToolsTest, GradedFractionsWithoutWellAreUnrefined )
{
    const auto fractions = RicWellLgrTools::gradedCumulativeFractions( { 50.0, 50.0 }, {}, 1.0, 1.5 );
    ASSERT_EQ( fractions.size(), 2u );
    EXPECT_EQ( fractions[0], std::vector<double>{ 1.0 } );
    EXPECT_EQ( fractions[1], std::vector<double>{ 1.0 } );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicWellLgrToolsTest, GradingKeywordsOnlyForGradedDirections )
{
    RigNonUniformRefinement refinement( cvf::Vec3st( 2, 1, 1 ) );
    refinement.setCumulativeFractions( RigRefinement::DimI, 0, { 0.5, 1.0 } );
    refinement.setCumulativeFractions( RigRefinement::DimI, 1, { 0.25, 0.75, 1.0 } );
    refinement.setCumulativeFractions( RigRefinement::DimK, 0, { 0.5, 1.0 } );

    const QString text = RicWellLgrTools::gradingKeywords( refinement );
    EXPECT_TRUE( text.contains( "NXFIN\n  2 3\n/" ) );
    EXPECT_TRUE( text.contains( "HXFIN\n  0.5 0.5 0.25 0.5 0.25\n/" ) );
    EXPECT_FALSE( text.contains( "NYFIN" ) );
    EXPECT_FALSE( text.contains( "NZFIN" ) );

    LgrInfo info( 1, "LGR1", "P1", cvf::Vec3st( 1, 1, 2 ), caf::VecIjk0( 4, 5, 0 ), caf::VecIjk0( 5, 5, 0 ) );
    info.gradedRefinement = std::make_shared<RigNonUniformRefinement>( refinement );

    const QString carfin = RicWellLgrTools::carfinBlock( info, 1 );
    EXPECT_TRUE( carfin.contains( "'LGR1' 5 6 6 6 1 1 5 1 2 1 /" ) );
    EXPECT_TRUE( carfin.trimmed().endsWith( "ENDFIN" ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicWellLgrToolsTest, DefaultLgrName )
{
    EXPECT_EQ( RicWellLgrTools::defaultLgrName( "B-2H" ), QString( "LB-2H" ) );
    EXPECT_EQ( RicWellLgrTools::defaultLgrName( "PRODUCER12" ), QString( "LPRODUCE" ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicFluxSectorModelToolsTest, ConvexHull )
{
    const std::vector<cvf::Vec3d> points = { { 0, 0, 1 }, { 2, 0, 2 }, { 2, 2, 3 }, { 0, 2, 4 }, { 1, 1, 5 }, { 1, 0, 6 } };
    const auto                    hull   = RicFluxSectorModelTools::convexHullXY( points, -10.0 );

    ASSERT_EQ( hull.size(), 4u );
    for ( const auto& p : hull )
    {
        EXPECT_DOUBLE_EQ( p.z(), -10.0 );
        EXPECT_TRUE( ( p.x() == 0.0 || p.x() == 2.0 ) && ( p.y() == 0.0 || p.y() == 2.0 ) );
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicFluxSectorModelToolsTest, RelocateFluxnumUsedByMultregt )
{
    QTemporaryDir dir;
    ASSERT_TRUE( dir.isValid() );
    const QString fileName = writeDeck( dir,
                                        "FLUXNUM\n 1 2 /\n"
                                        "EQUALS\n FLUXNUM 3 1 1 1 1 1 1 /\n PORO 0.25 /\n/\n"
                                        "MULTREGT\n 1 2 0.5 1* 1* 'F' /\n/\n" );

    RifOpmFlowDeckFile deck;
    ASSERT_TRUE( deck.loadDeck( fileName.toStdString() ).has_value() );

    auto result = RicFluxSectorModelTools::relocateDeckFluxnum( deck );
    ASSERT_TRUE( result.has_value() ) << result.error().toStdString();
    EXPECT_TRUE( result->contains( "MULTNUM" ) );

    EXPECT_TRUE( deck.findAllKeywords( "FLUXNUM" ).empty() );
    ASSERT_EQ( deck.findAllKeywords( "MULTNUM" ).size(), 1u );
    EXPECT_EQ( deck.findAllKeywords( "MULTNUM" )[0].getIntData(), ( std::vector<int>{ 1, 2 } ) );

    const auto equals = deck.findAllKeywords( "EQUALS" );
    ASSERT_EQ( equals.size(), 1u );
    ASSERT_EQ( equals[0].size(), 2u );
    EXPECT_EQ( equals[0].getRecord( 0 ).getItem( 0 ).getTrimmedString( 0 ), "MULTNUM" );
    EXPECT_EQ( equals[0].getRecord( 1 ).getItem( 0 ).getTrimmedString( 0 ), "PORO" );

    const auto multregt = deck.findAllKeywords( "MULTREGT" );
    ASSERT_EQ( multregt.size(), 1u );
    EXPECT_EQ( multregt[0].getRecord( 0 ).getItem( "REGION_DEF" ).getTrimmedString( 0 ), "M" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicFluxSectorModelToolsTest, RemoveUnusedFluxnum )
{
    QTemporaryDir dir;
    ASSERT_TRUE( dir.isValid() );
    const QString fileName = writeDeck( dir, "FLUXNUM\n 1 2 /\nEQUALS\n FLUXNUM 3 1 1 1 1 1 1 /\n/\n" );

    RifOpmFlowDeckFile deck;
    ASSERT_TRUE( deck.loadDeck( fileName.toStdString() ).has_value() );

    auto result = RicFluxSectorModelTools::relocateDeckFluxnum( deck );
    ASSERT_TRUE( result.has_value() );
    EXPECT_TRUE( result->contains( "removed" ) );
    EXPECT_TRUE( deck.findAllKeywords( "FLUXNUM" ).empty() );
    EXPECT_TRUE( deck.findAllKeywords( "MULTNUM" ).empty() );
    EXPECT_TRUE( deck.findAllKeywords( "EQUALS" ).empty() );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicFluxSectorModelToolsTest, DeckWithoutFluxnumIsUnchanged )
{
    QTemporaryDir dir;
    ASSERT_TRUE( dir.isValid() );
    const QString fileName = writeDeck( dir, "MULTREGT\n 1 2 0.5 /\n/\n" );

    RifOpmFlowDeckFile deck;
    ASSERT_TRUE( deck.loadDeck( fileName.toStdString() ).has_value() );

    auto result = RicFluxSectorModelTools::relocateDeckFluxnum( deck );
    ASSERT_TRUE( result.has_value() );
    EXPECT_TRUE( result->isEmpty() );
}
