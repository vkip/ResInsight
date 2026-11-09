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

#include "RicFluxSectorModelTools.h"

#include "RiaGuiApplication.h"
#include "RiaLogging.h"
#include "RiaPreferencesOpm.h"

#include "RicExportLgrFeature.h"
#include "RicWellLgrTools.h"

#include "RifOpmDeckTools.h"
#include "RifOpmFlowDeckFile.h"

#include "RigActiveCellInfo.h"
#include "RigCell.h"
#include "RigCellGeometryTools.h"
#include "RigEclipseCaseData.h"
#include "RigLocalGrid.h"
#include "RigMainGrid.h"
#include "SimulationFile/RigSimulationInputConnectionTools.h"
#include "SimulationFile/RigSimulationInputSettings.h"
#include "Well/RigSimWellData.h"

#include "CellFilters/RimCellFilterCollection.h"
#include "CellFilters/RimPolygonFilter.h"
#include "Polygons/RimPolygon.h"
#include "Polygons/RimPolygonCollection.h"
#include "RimEclipseCase.h"
#include "RimEclipseView.h"
#include "RimTools.h"

#include "cvfArray.h"

#include "opm/input/eclipse/Deck/DeckItem.hpp"
#include "opm/input/eclipse/Deck/DeckKeyword.hpp"
#include "opm/input/eclipse/Deck/DeckRecord.hpp"
#include "opm/input/eclipse/Deck/FileDeck.hpp"
#include "opm/input/eclipse/Parser/ParserKeywords/C.hpp"
#include "opm/input/eclipse/Parser/ParserKeywords/I.hpp"
#include "opm/input/eclipse/Parser/ParserKeywords/L.hpp"
#include "opm/input/eclipse/Parser/ParserKeywords/M.hpp"
#include "opm/input/eclipse/Parser/ParserKeywords/O.hpp"
#include "opm/input/eclipse/Parser/ParserKeywords/W.hpp"

#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QProcess>
#include <QRegularExpression>
#include <QStandardPaths>
#include <QTextStream>

#include <algorithm>
#include <map>
#include <set>
#include <tuple>

namespace
{
using CellSet = RicWellLgrTools::CellSet;

//--------------------------------------------------------------------------------------------------
/// Cells within 'padding' cells of the given cells in I and J, inside the grid
//--------------------------------------------------------------------------------------------------
void addLaterallyExpanded( const CellSet& cells, int padding, const RigMainGrid& mainGrid, CellSet& target )
{
    const int nx = static_cast<int>( mainGrid.cellCountI() );
    const int ny = static_cast<int>( mainGrid.cellCountJ() );
    for ( const auto& c : cells )
    {
        for ( int dj = -padding; dj <= padding; ++dj )
        {
            for ( int di = -padding; di <= padding; ++di )
            {
                const int i = static_cast<int>( c.i() ) + di;
                const int j = static_cast<int>( c.j() ) + dj;
                if ( i >= 0 && j >= 0 && i < nx && j < ny ) target.insert( caf::VecIjk0( i, j, c.k() ) );
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
CellSet boxCells( const caf::VecIjk0& boxMin, const caf::VecIjk0& boxMax )
{
    CellSet cells;
    for ( size_t k = boxMin.k(); k <= boxMax.k(); ++k )
        for ( size_t j = boxMin.j(); j <= boxMax.j(); ++j )
            for ( size_t i = boxMin.i(); i <= boxMax.i(); ++i )
                cells.insert( caf::VecIjk0( i, j, k ) );
    return cells;
}

//--------------------------------------------------------------------------------------------------
/// Host cells of the LGRs that are not part of the grid of the case, by name
//--------------------------------------------------------------------------------------------------
std::map<QString, CellSet> lgrHostCells( const RigMainGrid& mainGrid, const QStringList& excludedLgrs )
{
    std::map<QString, CellSet> lgrCells;
    for ( size_t c = 0; c < mainGrid.cellCount(); ++c )
    {
        const RigLocalGrid* subGrid = mainGrid.cell( c ).subGrid();
        if ( !subGrid ) continue;

        const QString name = QString::fromStdString( subGrid->gridName() );
        if ( excludedLgrs.contains( name ) ) continue;

        if ( auto ijk = mainGrid.ijkFromCellIndex( c ) ) lgrCells["LGR " + name].insert( *ijk );
    }
    return lgrCells;
}

//--------------------------------------------------------------------------------------------------
/// Run-length encoded GRDECL keyword of integers
//--------------------------------------------------------------------------------------------------
bool writeIntKeyword( const QString& fileName, const QString& keyword, const std::vector<int>& values, const QString& header )
{
    QFile file( fileName );
    if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) ) return false;

    QTextStream stream( &file );
    stream << header << keyword << "\n";

    int itemsOnLine = 0;
    for ( size_t n = 0; n < values.size(); )
    {
        size_t count = 1;
        while ( n + count < values.size() && values[n + count] == values[n] )
        {
            count++;
        }

        stream << ( itemsOnLine == 0 ? "  " : " " );
        if ( count > 1 )
            stream << count << "*" << values[n];
        else
            stream << values[n];

        if ( ++itemsOnLine == 10 )
        {
            stream << "\n";
            itemsOnLine = 0;
        }
        n += count;
    }
    if ( itemsOnLine > 0 ) stream << "\n";
    stream << "/\n";
    return true;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString sanitizedName( const QString& name )
{
    QString result = name;
    result.replace( QRegularExpression( "[^A-Za-z0-9_\\-]" ), "_" );
    return result;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<QString, QString> resolveDeckFile( const QString& gridFileName, const QString& deckFile )
{
    if ( !deckFile.isEmpty() )
    {
        if ( !QFileInfo::exists( deckFile ) ) return std::unexpected( QString( "The deck file '%1' does not exist." ).arg( deckFile ) );
        return QFileInfo( deckFile ).absoluteFilePath();
    }

    const QFileInfo gridInfo( gridFileName );
    for ( const QString& suffix : { ".DATA", ".data", ".Data" } )
    {
        const QString candidate = gridInfo.absolutePath() + "/" + gridInfo.completeBaseName() + suffix;
        if ( QFileInfo::exists( candidate ) ) return candidate;
    }

    if ( RiaGuiApplication::isRunning() )
    {
        const QString selected = QFileDialog::getOpenFileName( nullptr,
                                                               QString( "Select the deck of %1" ).arg( gridInfo.fileName() ),
                                                               gridInfo.absolutePath(),
                                                               "Eclipse/OPM Flow deck (*.DATA *.data);;All files (*)" );
        if ( !selected.isEmpty() ) return selected;
    }

    return std::unexpected(
        QString( "No deck file found next to '%1'. Specify it with deck_file." ).arg( QDir::toNativeSeparators( gridFileName ) ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<QString, QString> resolveMakeFluxExecutable( const QString& makeFluxExecutable )
{
    if ( !makeFluxExecutable.isEmpty() )
    {
        const QFileInfo info( makeFluxExecutable );
        if ( !info.isFile() || !info.isExecutable() )
            return std::unexpected( QString( "make_flux executable '%1' not found." ).arg( makeFluxExecutable ) );
        return info.absoluteFilePath();
    }

    const QString found = RicFluxSectorModelTools::defaultMakeFluxExecutable();
    if ( found.isEmpty() )
        return std::unexpected( "make_flux was found neither next to the OPM Flow executable of the preferences nor in PATH. "
                                "Specify it with make_flux_executable." );
    return found;
}

//--------------------------------------------------------------------------------------------------
/// The parent case for make_flux: the deck, with the results of the case next to it. A deck that does
/// not lie next to the grid file is copied to a folder holding links to the result files.
//--------------------------------------------------------------------------------------------------
std::expected<QString, QString> makeFluxParent( const QString& gridFileName, const QString& deckFile, const QString& exportFolder )
{
    const QFileInfo gridInfo( gridFileName );
    const QFileInfo deckInfo( deckFile );
    if ( gridInfo.absolutePath() == deckInfo.absolutePath() && gridInfo.completeBaseName() == deckInfo.completeBaseName() )
    {
        return deckInfo.absoluteFilePath();
    }

    const QString parentFolder = exportFolder + "/parent";
    QDir().mkpath( parentFolder );

    RifOpmFlowDeckFile parentDeck;
    if ( auto loaded = parentDeck.loadDeck( deckInfo.absoluteFilePath().toStdString() ); !loaded )
    {
        return std::unexpected( QString( "Failed to read the deck '%1': %2" ).arg( deckFile ).arg( QString::fromStdString( loaded.error() ) ) );
    }
    const QString parentDeckName = gridInfo.completeBaseName() + ".DATA";
    if ( !parentDeck.saveDeck( parentFolder.toStdString(), parentDeckName.toStdString(), /*copyAllIncludes*/ false ) )
    {
        return std::unexpected( QString( "Failed to copy the deck to '%1'." ).arg( parentFolder ) );
    }

    const QDir resultDir = gridInfo.absoluteDir();
    for ( const QString& resultFile : resultDir.entryList( { gridInfo.completeBaseName() + ".*" }, QDir::Files ) )
    {
        if ( resultFile.endsWith( ".DATA", Qt::CaseInsensitive ) ) continue;

        const QString link = parentFolder + "/" + resultFile;
        QFile::remove( link );
        QFile::link( resultDir.absoluteFilePath( resultFile ), link );
    }

    return parentFolder + "/" + parentDeckName;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<void, QString>
    runMakeFlux( const QString& executable, const QStringList& arguments, const QString& workingFolder, const QString& logFile )
{
    RiaLogging::info( QString( "Running %1 %2" ).arg( executable ).arg( arguments.join( " " ) ).toStdString() );

    QProcess process;
    process.setWorkingDirectory( workingFolder );
    process.setProcessChannelMode( QProcess::MergedChannels );
    process.setStandardOutputFile( logFile );
    process.start( executable, arguments );
    if ( !process.waitForStarted() ) return std::unexpected( QString( "Failed to start '%1'." ).arg( executable ) );

    process.waitForFinished( -1 );

    if ( process.exitStatus() != QProcess::NormalExit || process.exitCode() != 0 )
    {
        QString tail;
        QFile   log( logFile );
        if ( log.open( QIODevice::ReadOnly | QIODevice::Text ) )
        {
            QStringList lines = QString::fromUtf8( log.readAll() ).split( '\n', Qt::SkipEmptyParts );
            tail              = lines.mid( std::max( 0, static_cast<int>( lines.size() ) - 5 ) ).join( "\n" );
        }
        return std::unexpected(
            QString( "make_flux failed (exit code %1), see '%2':\n%3" ).arg( process.exitCode() ).arg( logFile ).arg( tail ) );
    }

    return {};
}

//--------------------------------------------------------------------------------------------------
/// A copy of the record with the string item at the given position replaced
//--------------------------------------------------------------------------------------------------
Opm::DeckRecord withStringItem( const Opm::DeckRecord& record, size_t itemIndex, const std::string& value )
{
    std::vector<Opm::DeckItem> items;
    for ( size_t n = 0; n < record.size(); ++n )
    {
        const auto& item = record.getItem( n );
        items.push_back( n == itemIndex ? RifOpmDeckTools::item( item.name(), value ) : item );
    }
    return Opm::DeckRecord{ std::move( items ) };
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::string upper( std::string text )
{
    std::transform( text.begin(), text.end(), text.begin(), ::toupper );
    return text;
}

//--------------------------------------------------------------------------------------------------
/// Record of a keyword with LGR support, e.g. WELSPECL, built from the corresponding global keyword record
//--------------------------------------------------------------------------------------------------
Opm::DeckRecord lgrRecord( const Opm::DeckRecord&                      globalRecord,
                           const std::vector<std::string>&             lgrItemNames,
                           const std::map<std::string, Opm::DeckItem>& replacements )
{
    std::vector<Opm::DeckItem> items;
    for ( const auto& name : lgrItemNames )
    {
        if ( auto it = replacements.find( name ); it != replacements.end() )
        {
            items.push_back( it->second );
        }
        else if ( globalRecord.hasItem( name ) )
        {
            items.push_back( globalRecord.getItem( name ) );
        }
        else
        {
            items.push_back( RifOpmDeckTools::defaultItem( name ) );
        }
    }
    return Opm::DeckRecord{ std::move( items ) };
}

//--------------------------------------------------------------------------------------------------
/// Position of the centre refined cell of a host cell, 1-based in the LGR
//--------------------------------------------------------------------------------------------------
int centreRefinedIndex( const RigRefinement& refinement, RigRefinement::Dimension dim, size_t boxIndex )
{
    return static_cast<int>( refinement.cumulativeOffset( dim, boxIndex ) + refinement.subcellCount( dim, boxIndex ) / 2 ) + 1;
}
} // namespace

//--------------------------------------------------------------------------------------------------
/// Andrew's monotone chain
//--------------------------------------------------------------------------------------------------
std::vector<cvf::Vec3d> RicFluxSectorModelTools::convexHullXY( const std::vector<cvf::Vec3d>& points, double z )
{
    std::vector<cvf::Vec3d> sorted = points;
    std::sort( sorted.begin(),
               sorted.end(),
               []( const cvf::Vec3d& a, const cvf::Vec3d& b ) { return a.x() < b.x() || ( a.x() == b.x() && a.y() < b.y() ); } );
    sorted.erase( std::unique( sorted.begin(),
                               sorted.end(),
                               []( const cvf::Vec3d& a, const cvf::Vec3d& b ) { return a.x() == b.x() && a.y() == b.y(); } ),
                  sorted.end() );
    if ( sorted.size() < 3 ) return sorted;

    auto cross = []( const cvf::Vec3d& o, const cvf::Vec3d& a, const cvf::Vec3d& b )
    { return ( a.x() - o.x() ) * ( b.y() - o.y() ) - ( a.y() - o.y() ) * ( b.x() - o.x() ); };

    std::vector<cvf::Vec3d> hull( 2 * sorted.size() );
    size_t                  k = 0;
    for ( size_t i = 0; i < sorted.size(); ++i )
    {
        while ( k >= 2 && cross( hull[k - 2], hull[k - 1], sorted[i] ) <= 0 )
            k--;
        hull[k++] = sorted[i];
    }
    for ( size_t i = sorted.size() - 1, t = k + 1; i > 0; --i )
    {
        while ( k >= t && cross( hull[k - 2], hull[k - 1], sorted[i - 1] ) <= 0 )
            k--;
        hull[k++] = sorted[i - 1];
    }
    hull.resize( k - 1 );

    for ( auto& p : hull )
    {
        p.z() = z;
    }
    return hull;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<RicFluxSectorModelTools::SectorRegion, QString> RicFluxSectorModelTools::computeSectorRegion( RimEclipseCase& eclipseCase,
                                                                                                            const Settings& settings )
{
    const RigEclipseCaseData* caseData = eclipseCase.eclipseCaseData();
    if ( !caseData || !caseData->mainGrid() ) return std::unexpected( "The case has no grid." );
    const RigMainGrid& mainGrid = *caseData->mainGrid();

    const auto wellCells = RicWellLgrTools::simulationWellCells( *caseData, settings.wellName );
    if ( wellCells.empty() ) return std::unexpected( QString( "Well '%1' has no connections in the main grid." ).arg( settings.wellName ) );

    CellSet required;
    addLaterallyExpanded( wellCells, settings.marginCells, mainGrid, required );

    // LGRs around the well are kept away from the sector boundary
    CellSet     lgrWithPadding;
    QStringList lgrNames;
    for ( const auto* wellLgr : RicWellLgrTools::temporaryLgrsForWell( mainGrid, settings.wellName ) )
    {
        const LgrInfo info = RicWellLgrTools::lgrInfoFromTemporaryLgr( mainGrid, *wellLgr );
        lgrNames.push_back( info.name );

        const int pad = settings.lgrPaddingCells;
        if ( static_cast<int>( info.mainGridStartCell.i() ) < pad || static_cast<int>( info.mainGridStartCell.j() ) < pad ||
             info.mainGridEndCell.i() + pad >= mainGrid.cellCountI() || info.mainGridEndCell.j() + pad >= mainGrid.cellCountJ() )
        {
            return std::unexpected( QString( "The LGR '%1' around well '%2' is within %3 cells of the boundary of the model. Sector models "
                                             "with an LGR on the sector boundary are not supported." )
                                        .arg( info.name )
                                        .arg( settings.wellName )
                                        .arg( pad ) );
        }

        caf::VecIjk0 paddedMin( info.mainGridStartCell.i() - pad,
                                info.mainGridStartCell.j() - pad,
                                info.mainGridStartCell.k() > static_cast<size_t>( pad ) ? info.mainGridStartCell.k() - pad : 0 );
        caf::VecIjk0 paddedMax( info.mainGridEndCell.i() + pad,
                                info.mainGridEndCell.j() + pad,
                                std::min( info.mainGridEndCell.k() + pad, mainGrid.cellCountK() - 1 ) );
        const auto   padded = boxCells( paddedMin, paddedMax );
        lgrWithPadding.insert( padded.begin(), padded.end() );
    }
    required.insert( lgrWithPadding.begin(), lgrWithPadding.end() );

    // Wells and LGRs must lie entirely inside or entirely outside the sector
    std::map<QString, CellSet> objects = lgrHostCells( mainGrid, lgrNames );
    for ( const auto& wellData : caseData->wellResults() )
    {
        if ( wellData->m_wellName == settings.wellName ) continue;
        auto cells = RicWellLgrTools::simulationWellCells( *caseData, wellData->m_wellName );
        if ( !cells.empty() ) objects["well " + wellData->m_wellName] = cells;
    }

    SectorRegion region;
    region.cells.assign( mainGrid.cellCount(), false );

    for ( int pass = 0;; ++pass )
    {
        if ( pass >= settings.maxExpansionPasses )
        {
            return std::unexpected( "The sector polygon kept growing to avoid cutting wells; no sector could be found." );
        }

        std::vector<cvf::Vec3d> corners;
        double                  topZ = -std::numeric_limits<double>::max();
        for ( const auto& c : required )
        {
            const size_t cellIndex = mainGrid.cellIndexFromIJK( c.i(), c.j(), c.k() );
            if ( mainGrid.cell( cellIndex ).isInvalid() ) continue;
            for ( const auto& corner : mainGrid.cellCornerVertices( cellIndex ) )
            {
                corners.push_back( corner );
                topZ = std::max( topZ, corner.z() );
            }
        }

        region.polygon = convexHullXY( corners, topZ );
        if ( region.polygon.size() < 3 ) return std::unexpected( "Could not build a polygon around the well." );

        auto closedPolygon = region.polygon;
        closedPolygon.push_back( closedPolygon.front() );

        for ( size_t c = 0; c < mainGrid.cellCount(); ++c )
        {
            const RigCell& cell = mainGrid.cell( c );
            region.cells[c]     = !cell.isInvalid() && RigCellGeometryTools::pointInsidePolygon2D( cell.center(), closedPolygon );
        }

        bool expanded = false;
        region.wellsInSector.clear();
        for ( const auto& [name, cells] : objects )
        {
            size_t inside = 0;
            for ( const auto& c : cells )
            {
                if ( region.cells[mainGrid.cellIndexFromIJK( c.i(), c.j(), c.k() )] ) inside++;
            }

            if ( inside > 0 && inside < cells.size() )
            {
                addLaterallyExpanded( cells, settings.wellPaddingCells, mainGrid, required );
                if ( !region.expandedForWells.contains( name ) ) region.expandedForWells.push_back( name );
                expanded = true;
            }
            else if ( inside > 0 && name.startsWith( "well " ) )
            {
                region.wellsInSector.push_back( name.mid( 5 ) );
            }
        }

        if ( !expanded ) break;
    }
    region.wellsInSector.push_front( settings.wellName );

    // The whole LGR and the cells around it must be in the sector
    const RigActiveCellInfo* activeInfo = caseData->activeCellInfo( RiaDefines::PorosityModelType::MATRIX_MODEL );
    for ( const auto& c : lgrWithPadding )
    {
        const size_t cellIndex = mainGrid.cellIndexFromIJK( c.i(), c.j(), c.k() );
        const bool   isActive  = activeInfo && activeInfo->isActive( ReservoirCellIndex( cellIndex ) );
        if ( isActive && !region.cells[cellIndex] )
        {
            return std::unexpected( QString( "The LGRs %1 around well '%2' and the %3 cells around them are not all inside the sector, at "
                                             "cell (%4, %5, %6). Sector models with an LGR on the sector boundary are not supported." )
                                        .arg( lgrNames.join( ", " ) )
                                        .arg( settings.wellName )
                                        .arg( settings.lgrPaddingCells )
                                        .arg( c.i() + 1 )
                                        .arg( c.j() + 1 )
                                        .arg( c.k() + 1 ) );
        }
    }

    return region;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<QString, QString> RicFluxSectorModelTools::relocateDeckFluxnum( RifOpmFlowDeckFile& deckFile )
{
    Opm::FileDeck* fileDeck = deckFile.fileDeck();
    if ( !fileDeck ) return std::unexpected( "No deck loaded." );

    const std::vector<std::string> operationKeywords = { "EQUALS", "ADD", "MULTIPLY", "MINVALUE", "MAXVALUE", "COPY" };
    auto                           refersTo = [&]( const std::string& keywordName, const Opm::DeckRecord& record, const std::string& array )
    {
        if ( keywordName == "COPY" )
            return upper( record.getItem( 0 ).getTrimmedString( 0 ) ) == array || upper( record.getItem( 1 ).getTrimmedString( 0 ) ) == array;
        return upper( record.getItem( 0 ).getTrimmedString( 0 ) ) == array;
    };
    auto isDefined = [&]( const std::string& array )
    {
        if ( !deckFile.findAllKeywords( array ).empty() ) return true;
        for ( const auto& keywordName : operationKeywords )
        {
            for ( const auto& keyword : deckFile.findAllKeywords( keywordName ) )
            {
                for ( const auto& record : keyword )
                {
                    if ( refersTo( keywordName, record, array ) ) return true;
                }
            }
        }
        return false;
    };

    for ( const std::string unsupported : { "OPERATE", "OPERATER" } )
    {
        for ( const auto& keyword : deckFile.findAllKeywords( unsupported ) )
        {
            for ( const auto& record : keyword )
            {
                for ( size_t n = 0; n < record.size(); ++n )
                {
                    const auto& item = record.getItem( n );
                    if ( item.getType() == Opm::type_tag::string && item.hasValue( 0 ) && upper( item.getTrimmedString( 0 ) ) == "FLUXNUM" )
                    {
                        return std::unexpected(
                            QString( "The deck uses FLUXNUM in %1, which cannot be moved to make room for the sector region." )
                                .arg( QString::fromStdString( unsupported ) ) );
                    }
                }
            }
        }
    }

    if ( !isDefined( "FLUXNUM" ) ) return QString();

    // Which region arrays the transmissibility multipliers between regions use
    std::set<std::string> multregtRegions;
    for ( const auto& keyword : deckFile.findAllKeywords( "MULTREGT" ) )
    {
        for ( const auto& record : keyword )
        {
            const auto& item = record.getItem( "REGION_DEF" );
            multregtRegions.insert( item.hasValue( 0 ) ? upper( item.getTrimmedString( 0 ) ) : "M" );
        }
    }

    std::string target;
    std::string targetLetter;
    if ( multregtRegions.contains( "F" ) )
    {
        if ( !isDefined( "MULTNUM" ) && !multregtRegions.contains( "M" ) )
        {
            target       = "MULTNUM";
            targetLetter = "M";
        }
        else if ( !isDefined( "OPERNUM" ) && !multregtRegions.contains( "O" ) )
        {
            target       = "OPERNUM";
            targetLetter = "O";
        }
        else
        {
            return std::unexpected( "MULTREGT uses FLUXNUM, and neither MULTNUM nor OPERNUM is free to take its regions." );
        }
    }

    // Renames FLUXNUM, or removes it when nothing refers to its regions
    std::vector<Opm::FileDeck::Index> emptied;
    for ( const auto& [index, keyword] : deckFile.findAllKeywordsWithIndices( "FLUXNUM" ) )
    {
        if ( target.empty() )
        {
            emptied.push_back( index );
            continue;
        }
        const auto&      data    = keyword.getIntData();
        Opm::DeckKeyword renamed = ( target == "MULTNUM" ) ? Opm::DeckKeyword( Opm::ParserKeywords::MULTNUM(), data )
                                                           : Opm::DeckKeyword( Opm::ParserKeywords::OPERNUM(), data );
        deckFile.replaceKeywordAtIndex( index, renamed );
    }

    for ( const auto& keywordName : operationKeywords )
    {
        for ( const auto& [index, keyword] : deckFile.findAllKeywordsWithIndices( keywordName ) )
        {
            bool             changed = false;
            Opm::DeckKeyword updated( keyword.location(), keyword.name() );
            for ( const auto& record : keyword )
            {
                if ( !refersTo( keywordName, record, "FLUXNUM" ) )
                {
                    updated.addRecord( Opm::DeckRecord( record ) );
                    continue;
                }
                changed = true;
                if ( target.empty() ) continue;

                Opm::DeckRecord renamed = record;
                for ( size_t itemIndex : { size_t( 0 ), size_t( 1 ) } )
                {
                    if ( itemIndex == 1 && keywordName != "COPY" ) continue;
                    if ( upper( renamed.getItem( itemIndex ).getTrimmedString( 0 ) ) == "FLUXNUM" )
                    {
                        renamed = withStringItem( renamed, itemIndex, target );
                    }
                }
                updated.addRecord( std::move( renamed ) );
            }

            if ( !changed ) continue;
            if ( updated.size() == 0 )
                emptied.push_back( index );
            else
                deckFile.replaceKeywordAtIndex( index, updated );
        }
    }

    if ( !target.empty() )
    {
        for ( const auto& [index, keyword] : deckFile.findAllKeywordsWithIndices( "MULTREGT" ) )
        {
            Opm::DeckKeyword updated( keyword.location(), keyword.name() );
            for ( const auto& record : keyword )
            {
                const auto& item = record.getItem( "REGION_DEF" );
                if ( item.hasValue( 0 ) && upper( item.getTrimmedString( 0 ) ) == "F" )
                {
                    updated.addRecord( withStringItem( record, record.size() - 1, targetLetter ) );
                }
                else
                {
                    updated.addRecord( Opm::DeckRecord( record ) );
                }
            }
            deckFile.replaceKeywordAtIndex( index, updated );
        }
    }

    // Erase from the back, so the remaining indices stay valid
    std::sort( emptied.begin(), emptied.end() );
    for ( auto it = emptied.rbegin(); it != emptied.rend(); ++it )
    {
        deckFile.removeKeywordAtIndex( *it );
    }

    if ( target.empty() ) return QString( "The deck's own FLUXNUM was removed, nothing refers to its regions." );
    return QString( "The deck's own FLUXNUM was renamed to %1, and MULTREGT uses it as region array '%2'." )
        .arg( QString::fromStdString( target ) )
        .arg( QString::fromStdString( targetLetter ) );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<int>
    RicFluxSectorModelTools::removeAquifersOutsideSector( RifOpmFlowDeckFile& deckFile, const std::vector<int>& fluxnum, size_t nx, size_t ny )
{
    const auto aqucons = deckFile.findAllKeywordsWithIndices( "AQUCON" );
    const auto aqunums = deckFile.findAllKeywordsWithIndices( "AQUNUM" );
    if ( aqunums.empty() ) return {};

    const auto intItem  = []( const Opm::DeckRecord& record, size_t index ) { return record.getItem( index ).get<int>( 0 ); };
    const auto inSector = [&]( int i, int j, int k )
    {
        const size_t cell = static_cast<size_t>( i - 1 ) + nx * ( static_cast<size_t>( j - 1 ) + ny * static_cast<size_t>( k - 1 ) );
        return cell < fluxnum.size() && fluxnum[cell] == 1;
    };

    // Aquifers with a connection (AQUCON: id, I1, I2, J1, J2, K1, K2, ...) to a cell of the sector
    std::set<int> connected;
    for ( const auto& [index, keyword] : aqucons )
    {
        for ( const auto& record : keyword )
        {
            const int id = intItem( record, 0 );
            for ( int k = intItem( record, 5 ); k <= intItem( record, 6 ) && !connected.contains( id ); ++k )
                for ( int j = intItem( record, 3 ); j <= intItem( record, 4 ) && !connected.contains( id ); ++j )
                    for ( int i = intItem( record, 1 ); i <= intItem( record, 2 ); ++i )
                    {
                        if ( inSector( i, j, k ) )
                        {
                            connected.insert( id );
                            break;
                        }
                    }
        }
    }

    std::set<int>                                                  removed;
    std::vector<Opm::FileDeck::Index>                              emptied;
    std::vector<std::pair<Opm::FileDeck::Index, Opm::DeckKeyword>> updates;
    for ( const auto& keywords : { aqunums, aqucons } )
    {
        for ( const auto& [index, keyword] : keywords )
        {
            Opm::DeckKeyword updated( keyword.location(), keyword.name() );
            for ( const auto& record : keyword )
            {
                const int id = intItem( record, 0 );
                if ( connected.contains( id ) )
                    updated.addRecord( Opm::DeckRecord( record ) );
                else
                    removed.insert( id );
            }
            if ( updated.size() == keyword.size() ) continue;
            if ( updated.size() == 0 )
                emptied.push_back( index );
            else
                updates.emplace_back( index, updated );
        }
    }

    // Summary vectors of the removed aquifers (ANQR, ANQT, ANQP: a list of aquifer numbers, or all when empty)
    if ( !removed.empty() )
    {
        for ( const std::string name : { "ANQR", "ANQT", "ANQP" } )
        {
            for ( const auto& [index, keyword] : deckFile.findAllKeywordsWithIndices( name ) )
            {
                if ( !keyword.isDataKeyword() || keyword.size() == 0 ) continue;
                Opm::DeckKeyword updated( keyword );
                auto&            ids   = updated.getIntData();
                const size_t     count = ids.size();
                if ( count == 0 ) continue;
                std::erase_if( ids, [&]( int id ) { return removed.contains( id ); } );
                if ( ids.size() == count ) continue;
                if ( ids.empty() )
                    emptied.push_back( index );
                else
                    updates.emplace_back( index, updated );
            }
        }
    }

    for ( const auto& [index, keyword] : updates )
    {
        deckFile.replaceKeywordAtIndex( index, keyword );
    }
    // Erase from the back, so the remaining indices stay valid
    std::sort( emptied.begin(), emptied.end() );
    for ( auto it = emptied.rbegin(); it != emptied.rend(); ++it )
    {
        deckFile.removeKeywordAtIndex( *it );
    }
    return { removed.begin(), removed.end() };
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool RicFluxSectorModelTools::isMultiSegmentWell( const RifOpmFlowDeckFile& deckFile, const QString& wellName )
{
    const std::string well = wellName.toStdString();
    for ( const std::string multiSegment : { "WELSEGS", "COMPSEGS" } )
    {
        for ( const auto& keyword : const_cast<RifOpmFlowDeckFile&>( deckFile ).findAllKeywords( multiSegment ) )
        {
            if ( keyword.size() > 0 && keyword.getRecord( 0 ).getItem( 0 ).getTrimmedString( 0 ) == well ) return true;
        }
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<int, QString> RicFluxSectorModelTools::convertWellToLgr( RifOpmFlowDeckFile& deckFile,
                                                                       RimEclipseCase*     eclipseCase,
                                                                       const QString&      wellName,
                                                                       const LgrInfo&      lgrInfo )
{
    const std::string well    = wellName.toStdString();
    const std::string lgrName = lgrInfo.name.toStdString();

    if ( isMultiSegmentWell( deckFile, wellName ) )
    {
        return std::unexpected( QString( "Well '%1' is a multi-segment well, which cannot be converted to COMPDATL yet." ).arg( wellName ) );
    }
    for ( const std::string connectionKeyword : { "COMPLUMP", "WPIMULT", "WELOPEN", "CSKIN" } )
    {
        for ( const auto& keyword : deckFile.findAllKeywords( connectionKeyword ) )
        {
            for ( const auto& record : keyword )
            {
                if ( record.getItem( 0 ).getTrimmedString( 0 ) != well ) continue;
                for ( const std::string item : { "I", "J", "K" } )
                {
                    if ( record.hasItem( item ) && record.getItem( item ).hasValue( 0 ) && record.getItem( item ).get<int>( 0 ) > 0 )
                    {
                        RiaLogging::warning(
                            QString( "%1 refers to connections of well '%2' by cell; it is not converted to the LGR '%3'." )
                                .arg( QString::fromStdString( connectionKeyword ) )
                                .arg( wellName )
                                .arg( lgrInfo.name )
                                .toStdString() );
                        break;
                    }
                }
            }
        }
    }

    const auto                 refinement = RicWellLgrTools::lgrRefinement( lgrInfo );
    const auto&                boxMin     = lgrInfo.mainGridStartCell;
    const auto&                boxMax     = lgrInfo.mainGridEndCell;
    const auto                 cellProps  = RigSimulationInputConnectionTools::connectionCellPropertiesProvider( eclipseCase );
    RigSimulationInputSettings settings;
    settings.setMin( boxMin );
    settings.setMax( boxMax );
    settings.setRefinement( refinement->clone() );

    auto insideBox = [&]( int i, int j )
    {
        return i >= static_cast<int>( boxMin.i() ) && i <= static_cast<int>( boxMax.i() ) && j >= static_cast<int>( boxMin.j() ) &&
               j <= static_cast<int>( boxMax.j() );
    };

    const auto wellCells    = RicWellLgrTools::simulationWellCells( *eclipseCase->eclipseCaseData(), wellName );
    const auto fallbackHead = wellCells.empty() ? caf::VecIjk0( boxMin.i(), boxMin.j(), 0 )
                                                : caf::VecIjk0( wellCells.begin()->i(), wellCells.begin()->j(), 0 );

    using WS                                     = Opm::ParserKeywords::WELSPECS;
    using CD                                     = Opm::ParserKeywords::COMPDAT;
    const std::vector<std::string> welspeclItems = { "WELL",
                                                     "GROUP",
                                                     "LGR",
                                                     "HEAD_I",
                                                     "HEAD_J",
                                                     "REF_DEPTH",
                                                     "PHASE",
                                                     "D_RADIUS",
                                                     "INFLOW_EQ",
                                                     "AUTO_SHUTIN",
                                                     "CROSSFLOW",
                                                     "P_TABLE",
                                                     "DENSITY_CALC",
                                                     "FIP_REGION",
                                                     "FRONTSIM1",
                                                     "FRONTSIM2",
                                                     "well_model",
                                                     "POLYMER_TABLE" };
    const std::vector<std::string> compdatlItems = { "WELL",
                                                     "LGR",
                                                     "I",
                                                     "J",
                                                     "K1",
                                                     "K2",
                                                     "STATE",
                                                     "SAT_TABLE",
                                                     "CONNECTION_TRANSMISSIBILITY_FACTOR",
                                                     "DIAMETER",
                                                     "Kh",
                                                     "SKIN",
                                                     "D_FACTOR",
                                                     "DIR",
                                                     "PR" };

    // Keywords of the well in deck order, so a defaulted COMPDAT cell takes the well head in effect
    struct Entry
    {
        Opm::FileDeck::Index index;
        Opm::DeckKeyword     keyword;
    };
    std::vector<Entry> entries;
    for ( const std::string name : { "WELSPECS", "COMPDAT" } )
    {
        for ( const auto& [index, keyword] : deckFile.findAllKeywordsWithIndices( name ) )
        {
            for ( const auto& record : keyword )
            {
                if ( record.getItem( 0 ).getTrimmedString( 0 ) == well )
                {
                    entries.push_back( { index, keyword } );
                    break;
                }
            }
        }
    }
    std::sort( entries.begin(), entries.end(), []( const Entry& a, const Entry& b ) { return a.index < b.index; } );
    if ( entries.empty() ) return std::unexpected( QString( "No WELSPECS or COMPDAT for well '%1' in the deck." ).arg( wellName ) );

    struct Replacement
    {
        Opm::FileDeck::Index            index;
        std::optional<Opm::DeckKeyword> globalKeyword;
        Opm::DeckKeyword                lgrKeyword;
    };
    std::vector<Replacement> replacements;

    std::optional<std::pair<int, int>>  head; // 1-based in the main grid
    int                                 convertedRecords = 0;
    std::set<std::tuple<int, int, int>> refinedConnections;
    for ( const auto& entry : entries )
    {
        const bool       isWelspecs = entry.keyword.name() == "WELSPECS";
        Opm::DeckKeyword globalKeyword( entry.keyword.location(), entry.keyword.name() );
        Opm::DeckKeyword lgrKeyword( entry.keyword.location(), isWelspecs ? "WELSPECL" : "COMPDATL" );

        for ( const auto& record : entry.keyword )
        {
            if ( record.getItem( 0 ).getTrimmedString( 0 ) != well )
            {
                globalKeyword.addRecord( Opm::DeckRecord( record ) );
                continue;
            }

            convertedRecords++;
            if ( isWelspecs )
            {
                const int headI = record.getItem( WS::HEAD_I::itemName ).get<int>( 0 );
                const int headJ = record.getItem( WS::HEAD_J::itemName ).get<int>( 0 );
                head            = std::make_pair( headI, headJ );

                caf::VecIjk0 localHead = insideBox( headI - 1, headJ - 1 ) ? caf::VecIjk0( headI - 1, headJ - 1, 0 ) : fallbackHead;
                if ( !insideBox( headI - 1, headJ - 1 ) )
                {
                    RiaLogging::warning( QString( "The head of well '%1' at (%2, %3) is outside the LGR '%4'; it is moved to (%5, %6)." )
                                             .arg( wellName )
                                             .arg( headI )
                                             .arg( headJ )
                                             .arg( lgrInfo.name )
                                             .arg( localHead.i() + 1 )
                                             .arg( localHead.j() + 1 )
                                             .toStdString() );
                }

                std::map<std::string, Opm::DeckItem> items =
                    { { "LGR", RifOpmDeckTools::item( "LGR", lgrName ) },
                      { "HEAD_I",
                        RifOpmDeckTools::item( "HEAD_I", centreRefinedIndex( *refinement, RigRefinement::DimI, localHead.i() - boxMin.i() ) ) },
                      { "HEAD_J",
                        RifOpmDeckTools::item( "HEAD_J", centreRefinedIndex( *refinement, RigRefinement::DimJ, localHead.j() - boxMin.j() ) ) } };
                lgrKeyword.addRecord( lgrRecord( record, welspeclItems, items ) );
                continue;
            }

            // A defaulted cell of COMPDAT is the well head
            Opm::DeckRecord explicitRecord = record;
            for ( auto [itemName, value] :
                  { std::make_pair( CD::I::itemName, head ? head->first : 0 ), std::make_pair( CD::J::itemName, head ? head->second : 0 ) } )
            {
                const auto& item = explicitRecord.getItem( itemName );
                if ( !item.hasValue( 0 ) || item.get<int>( 0 ) <= 0 )
                {
                    std::vector<Opm::DeckItem> items;
                    for ( size_t n = 0; n < explicitRecord.size(); ++n )
                    {
                        const auto& existing = explicitRecord.getItem( n );
                        items.push_back( existing.name() == itemName ? RifOpmDeckTools::item( itemName, value ) : existing );
                    }
                    explicitRecord = Opm::DeckRecord{ std::move( items ) };
                }
            }

            auto refined = RigSimulationInputConnectionTools::processCompdatRecord( explicitRecord, well, settings, cellProps );
            if ( !refined )
            {
                return std::unexpected( QString( "COMPDAT of well '%1' could not be moved to the LGR '%2': %3" )
                                            .arg( wellName )
                                            .arg( lgrInfo.name )
                                            .arg( refined.error() ) );
            }
            for ( const auto& refinedRecord : *refined )
            {
                using CDI = Opm::ParserKeywords::COMPDAT;
                for ( int k = refinedRecord.getItem( CDI::K1::itemName ).get<int>( 0 );
                      k <= refinedRecord.getItem( CDI::K2::itemName ).get<int>( 0 );
                      ++k )
                {
                    refinedConnections.insert( { refinedRecord.getItem( CDI::I::itemName ).get<int>( 0 ),
                                                 refinedRecord.getItem( CDI::J::itemName ).get<int>( 0 ),
                                                 k } );
                }
                lgrKeyword.addRecord( lgrRecord( refinedRecord, compdatlItems, { { "LGR", RifOpmDeckTools::item( "LGR", lgrName ) } } ) );
            }
        }

        replacements.push_back(
            { entry.index, globalKeyword.size() > 0 ? std::optional<Opm::DeckKeyword>( globalKeyword ) : std::nullopt, lgrKeyword } );
    }

    // From the back, so inserting a keyword does not move those not yet replaced
    Opm::FileDeck* fileDeck = deckFile.fileDeck();
    for ( auto it = replacements.rbegin(); it != replacements.rend(); ++it )
    {
        if ( it->globalKeyword )
        {
            deckFile.replaceKeywordAtIndex( it->index, *it->globalKeyword );
            fileDeck->insert( it->index + 1, it->lgrKeyword );
        }
        else
        {
            deckFile.replaceKeywordAtIndex( it->index, it->lgrKeyword );
        }
    }

    // Splitting connections over refined cells may need room for more connections per well
    if ( auto dims = deckFile.welldims(); dims.size() >= 4 && static_cast<int>( refinedConnections.size() ) > dims[1] )
    {
        deckFile.setWelldims( dims[0], static_cast<int>( refinedConnections.size() ), dims[2], dims[3] );
        RiaLogging::info( QString( "WELLDIMS: maximum connections per well raised from %1 to %2 for well '%3' in the LGR." )
                              .arg( dims[1] )
                              .arg( refinedConnections.size() )
                              .arg( wellName )
                              .toStdString() );
    }

    return convertedRecords;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<void, QString> RicFluxSectorModelTools::addLgrDimensions( RifOpmFlowDeckFile& deckFile, const std::vector<LgrInfo>& lgrInfos )
{
    using L = Opm::ParserKeywords::LGR;

    int lgrCount  = static_cast<int>( lgrInfos.size() );
    int cellCount = 0;
    for ( const auto& lgrInfo : lgrInfos )
    {
        cellCount += static_cast<int>( lgrInfo.lgrCellCount() );
    }
    if ( auto existing = deckFile.findKeyword( "LGR" ) )
    {
        const auto& record = existing->getRecord( 0 );
        if ( record.getItem( L::MAXLGR::itemName ).hasValue( 0 ) ) lgrCount += record.getItem( L::MAXLGR::itemName ).get<int>( 0 );
        if ( record.getItem( L::MAXCLS::itemName ).hasValue( 0 ) ) cellCount += record.getItem( L::MAXCLS::itemName ).get<int>( 0 );
    }

    std::map<QString, int> amalgamations;
    for ( const auto& lgrInfo : lgrInfos )
    {
        if ( !lgrInfo.amalgamationName.isEmpty() ) amalgamations[lgrInfo.amalgamationName]++;
    }
    int largestAmalgamation = 0;
    for ( const auto& [name, count] : amalgamations )
    {
        largestAmalgamation = std::max( largestAmalgamation, count );
    }

    std::vector<Opm::DeckItem> items = { RifOpmDeckTools::item( L::MAXLGR::itemName, lgrCount ),
                                         RifOpmDeckTools::item( L::MAXCLS::itemName, cellCount ) };
    if ( !amalgamations.empty() )
    {
        items.push_back( RifOpmDeckTools::defaultItem( L::MCOARS::itemName ) );
        items.push_back( RifOpmDeckTools::item( L::MAMALG::itemName, static_cast<int>( amalgamations.size() ) ) );
        items.push_back( RifOpmDeckTools::item( L::MXLALG::itemName, largestAmalgamation ) );
    }

    Opm::DeckKeyword keyword( ( L() ) );
    keyword.addRecord( Opm::DeckRecord{ std::move( items ) } );

    if ( !deckFile.replaceKeyword( "RUNSPEC", keyword ) ) return std::unexpected( "The deck has no RUNSPEC section." );
    return {};
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RicFluxSectorModelTools::defaultExportFolder( const QString& gridFileName, const QString& wellName )
{
    return QFileInfo( gridFileName ).absolutePath() + "/sector_" + sanitizedName( wellName );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
QString RicFluxSectorModelTools::defaultMakeFluxExecutable()
{
    if ( const auto* preferences = RiaPreferencesOpm::current() )
    {
        const QFileInfo flow( preferences->opmFlowCommand() );
        if ( !preferences->useWsl() && flow.isAbsolute() )
        {
            const QFileInfo candidate( flow.absolutePath() + "/make_flux" );
            if ( candidate.isFile() && candidate.isExecutable() ) return candidate.absoluteFilePath();
        }
    }
    return QStandardPaths::findExecutable( "make_flux" );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::expected<RicFluxSectorModelTools::Result, QString> RicFluxSectorModelTools::createSectorModel( RimEclipseView* view,
                                                                                                    const Settings& settings )
{
    if ( !view ) return std::unexpected( "No view is available." );
    RimEclipseCase* eclipseCase = view->eclipseCase();
    if ( !eclipseCase || !eclipseCase->eclipseCaseData() || !eclipseCase->eclipseCaseData()->mainGrid() )
        return std::unexpected( "The view has no case with a grid." );
    if ( settings.wellName.isEmpty() ) return std::unexpected( "No well name specified." );
    if ( !eclipseCase->eclipseCaseData()->findSimWellData( settings.wellName ) )
        return std::unexpected( QString( "The case has no well named '%1'." ).arg( settings.wellName ) );

    const RigMainGrid& mainGrid     = *eclipseCase->eclipseCaseData()->mainGrid();
    const QString      gridFileName = eclipseCase->gridFileName();
    const QFileInfo    gridInfo( gridFileName );

    auto deckFile = resolveDeckFile( gridFileName, settings.deckFile );
    if ( !deckFile ) return std::unexpected( deckFile.error() );

    QString makeFlux;
    if ( settings.runMakeFlux )
    {
        auto executable = resolveMakeFluxExecutable( settings.makeFluxExecutable );
        if ( !executable ) return std::unexpected( executable.error() );
        makeFlux = *executable;
    }

    const QString exportFolder = settings.exportFolder.isEmpty() ? defaultExportFolder( gridFileName, settings.wellName )
                                                                 : QFileInfo( settings.exportFolder ).absoluteFilePath();
    if ( !QDir().mkpath( exportFolder ) ) return std::unexpected( QString( "Could not create the folder '%1'." ).arg( exportFolder ) );

    auto region = computeSectorRegion( *eclipseCase, settings );
    if ( !region ) return std::unexpected( region.error() );

    Result result;
    result.wellsInSector = region->wellsInSector;
    for ( const auto& name : region->expandedForWells )
    {
        result.warnings.push_back( QString( "The sector was enlarged to hold all of %1." ).arg( name ) );
    }

    // Polygon and polygon filter in the view
    RimPolygonCollection* polygonCollection = RimTools::polygonCollection();
    if ( !polygonCollection ) return std::unexpected( "The project has no polygon collection." );

    RimPolygon* polygon = polygonCollection->createUserDefinedPolygon();
    polygon->setName( QString( "Sector %1" ).arg( settings.wellName ) );
    polygon->setPointsInDomainCoords( region->polygon );
    polygon->setIsClosed( true );
    polygonCollection->addUserDefinedPolygon( polygon );
    polygonCollection->uiCapability()->updateAllRequiredEditors();

    RimPolygonFilter* filter = view->cellFilterCollection()->addNewPolygonFilter( eclipseCase, polygon );
    filter->setName( QString( "Sector %1" ).arg( settings.wellName ) );
    filter->setPolygonFilterMode( RimPolygonFilter::PolygonFilterModeType::DEPTH_Z, RimPolygonFilter::PolygonIncludeType::CENTER );

    // The FLUXNUM region is the main grid cells the filter shows
    cvf::UByteArray includeCells( mainGrid.cellCount() );
    cvf::UByteArray excludeCells( mainGrid.cellCount() );
    includeCells.setAll( false );
    excludeCells.setAll( true );
    filter->updateCellIndexFilter( &includeCells, &excludeCells, 0 );

    std::vector<int> fluxnum( mainGrid.cellCount(), 0 );
    size_t           mismatches = 0;
    for ( size_t c = 0; c < mainGrid.cellCount(); ++c )
    {
        fluxnum[c] = includeCells[c] ? 1 : 0;
        if ( static_cast<bool>( includeCells[c] ) != region->cells[c] ) mismatches++;
        if ( fluxnum[c] == 1 )
        {
            result.sectorCellCount++;
            if ( eclipseCase->eclipseCaseData()->activeCellInfo( RiaDefines::PorosityModelType::MATRIX_MODEL )->isActive( ReservoirCellIndex( c ) ) )
            {
                result.activeSectorCellCount++;
            }
        }
    }
    if ( mismatches > 0 )
    {
        return std::unexpected( QString( "The polygon filter selects %1 cells that differ from the sector region." ).arg( mismatches ) );
    }

    view->scheduleCreateDisplayModelAndRedraw();

    const QString baseName = sanitizedName( gridInfo.completeBaseName() + "_SECTOR_" + settings.wellName ).toUpper();
    result.fluxnumFile     = exportFolder + "/" + baseName + "_FLUXNUM.grdecl";
    result.fluxFile        = exportFolder + "/" + baseName + ".FLUX";
    result.deckFile        = exportFolder + "/" + baseName + ".DATA";

    const QString header =
        QString( "-- Sector around well %1 in %2, written by ResInsight\n" ).arg( settings.wellName ).arg( gridInfo.fileName() );
    if ( !writeIntKeyword( result.fluxnumFile, "FLUXNUM", fluxnum, header ) )
    {
        return std::unexpected( QString( "Could not write '%1'." ).arg( result.fluxnumFile ) );
    }

    if ( settings.runMakeFlux )
    {
        auto parent = makeFluxParent( gridFileName, *deckFile, exportFolder );
        if ( !parent ) return std::unexpected( parent.error() );

        const QStringList arguments = { "--parent=" + *parent,
                                        "--fluxnum=" + result.fluxnumFile,
                                        "--regions=1",
                                        "--output=" + result.fluxFile,
                                        "--parsing-strictness=low" };
        if ( auto ran = runMakeFlux( makeFlux, arguments, exportFolder, exportFolder + "/" + baseName + "_make_flux.log" ); !ran )
        {
            return std::unexpected( ran.error() );
        }
    }

    // The sector deck
    RifOpmFlowDeckFile deck;
    if ( auto loaded = deck.loadDeck( deckFile->toStdString() ); !loaded )
    {
        return std::unexpected( QString( "Failed to read the deck '%1': %2" ).arg( *deckFile ).arg( QString::fromStdString( loaded.error() ) ) );
    }

    auto relocated = relocateDeckFluxnum( deck );
    if ( !relocated ) return std::unexpected( relocated.error() );
    if ( !relocated->isEmpty() ) result.warnings.push_back( *relocated );

    if ( const auto aquifers = removeAquifersOutsideSector( deck, fluxnum, mainGrid.cellCountI(), mainGrid.cellCountJ() ); !aquifers.empty() )
    {
        QStringList numbers;
        for ( int id : aquifers )
        {
            numbers.push_back( QString::number( id ) );
        }
        result.warnings.push_back(
            QString( "The numerical aquifers %1 connect to no cell of the sector and were removed from the deck." ).arg( numbers.join( ", " ) ) );
    }

    QString gridText = header + QString( "\nINCLUDE\n  '%1' /\n\n"
                                         "USEFLUX\n  '%2' /\n\n"
                                         "FLUXTYPE\n  'PRESSURE' /\n\n"
                                         "FLUXREG\n  1 /\n" )
                                    .arg( QFileInfo( result.fluxnumFile ).fileName() )
                                    .arg( baseName );

    std::vector<LgrInfo> lgrInfos;
    for ( const auto* wellLgr : RicWellLgrTools::temporaryLgrsForWell( mainGrid, settings.wellName ) )
    {
        lgrInfos.push_back( RicWellLgrTools::lgrInfoFromTemporaryLgr( mainGrid, *wellLgr ) );
    }
    if ( !lgrInfos.empty() )
    {
        for ( const auto& lgrInfo : lgrInfos )
        {
            result.lgrNames.push_back( lgrInfo.name );
            gridText += "\n" + RicWellLgrTools::carfinBlock( lgrInfo,
                                                             RicWellLgrTools::wellCountInBox( *eclipseCase->eclipseCaseData(),
                                                                                              lgrInfo.mainGridStartCell,
                                                                                              lgrInfo.mainGridEndCell ) );
        }

        gridText += RicWellLgrTools::amalgamKeyword( lgrInfos );

        if ( auto added = addLgrDimensions( deck, lgrInfos ); !added ) return std::unexpected( added.error() );

        // A well in several LGRs, partly outside its LGR or with segments keeps its COMPDAT (and COMPSEGS)
        const auto& lgr = lgrInfos.front();
        const bool  allInOneLgr =
            lgrInfos.size() == 1 &&
            std::ranges::all_of( RicWellLgrTools::simulationWellCells( *eclipseCase->eclipseCaseData(), settings.wellName ),
                                 [&lgr]( const caf::VecIjk0& c )
                                 {
                                     return c.i() >= lgr.mainGridStartCell.i() && c.i() <= lgr.mainGridEndCell.i() &&
                                            c.j() >= lgr.mainGridStartCell.j() && c.j() <= lgr.mainGridEndCell.j() &&
                                            c.k() >= lgr.mainGridStartCell.k() && c.k() <= lgr.mainGridEndCell.k();
                                 } );
        const bool multiSegment = isMultiSegmentWell( deck, settings.wellName );
        if ( allInOneLgr && !multiSegment )
        {
            auto converted = convertWellToLgr( deck, eclipseCase, settings.wellName, lgrInfos.front() );
            if ( !converted ) return std::unexpected( converted.error() );
        }
        else
        {
            result.warnings.push_back( QString( "Well '%1' %2, so its connections are kept as COMPDAT. OPM Flow moves the connections "
                                                "in the LGRs %3 into them; other simulators need WELSPECL/COMPDATL%4." )
                                           .arg( settings.wellName )
                                           .arg( multiSegment          ? "is a multi-segment well"
                                                 : allInOneLgr         ? ""
                                                 : lgrInfos.size() > 1 ? "has several LGRs"
                                                                       : "has connections outside its LGR" )
                                           .arg( result.lgrNames.join( ", " ) )
                                           .arg( multiSegment ? "/COMPSEGL" : "" ) );
        }
    }

    const QString gridIncludeFile = exportFolder + "/" + baseName + "_GRID.inc";
    {
        QFile file( gridIncludeFile );
        if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) )
            return std::unexpected( QString( "Could not write '%1'." ).arg( gridIncludeFile ) );
        QTextStream( &file ) << gridText;
    }

    Opm::DeckKeyword includeKeyword( ( Opm::ParserKeywords::INCLUDE() ) );
    includeKeyword.addRecord( Opm::DeckRecord{ { RifOpmDeckTools::item( Opm::ParserKeywords::INCLUDE::IncludeFile::itemName,
                                                                        QFileInfo( gridIncludeFile ).fileName().toStdString() ) } } );
    if ( !deck.addKeyword( "GRID", includeKeyword ) ) return std::unexpected( "The deck has no GRID section." );

    if ( !deck.saveDeck( exportFolder.toStdString(), QFileInfo( result.deckFile ).fileName().toStdString(), /*copyAllIncludes*/ false ) )
    {
        return std::unexpected( QString( "Failed to write the sector deck '%1'." ).arg( result.deckFile ) );
    }

    RiaLogging::info( QString( "Sector model around well '%1': %2 cells (%3 active), wells %4. Deck: %5" )
                          .arg( settings.wellName )
                          .arg( result.sectorCellCount )
                          .arg( result.activeSectorCellCount )
                          .arg( result.wellsInSector.join( ", " ) )
                          .arg( result.deckFile )
                          .toStdString() );
    for ( const auto& warning : result.warnings )
    {
        RiaLogging::warning( warning.toStdString() );
    }

    return result;
}
