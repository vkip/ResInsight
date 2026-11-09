#include "gtest/gtest.h"

#include "ExportCommands/RicCreateSectorModelUi.h"
#include "ExportCommands/RicFluxSectorModelTools.h"

#include "cafCmdFeature.h"
#include "cafFactory.h"

#include <algorithm>

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicCreateSectorModelUiTest, DefaultsAreThoseOfTheSectorModel )
{
    RicCreateSectorModelUi                  ui;
    const RicFluxSectorModelTools::Settings defaults;
    const auto                              settings = ui.settings();

    EXPECT_EQ( defaults.marginCells, settings.marginCells );
    EXPECT_EQ( defaults.wellPaddingCells, settings.wellPaddingCells );
    EXPECT_EQ( defaults.lgrPaddingCells, settings.lgrPaddingCells );
    EXPECT_EQ( defaults.runMakeFlux, settings.runMakeFlux );
    EXPECT_TRUE( settings.exportFolder.isEmpty() );
    EXPECT_TRUE( settings.deckFile.isEmpty() );
}

//--------------------------------------------------------------------------------------------------
/// The settings not tied to a well carry over from one sector model to the next; the well's own do not
//--------------------------------------------------------------------------------------------------
TEST( RicCreateSectorModelUiTest, SharedSettingsCarryOver )
{
    RicFluxSectorModelTools::Settings previous;
    previous.wellName           = "A4";
    previous.exportFolder       = "/data/sector_A4";
    previous.deckFile           = "/data/CASE.DATA";
    previous.makeFluxExecutable = "/opt/opm/bin/make_flux";
    previous.runMakeFlux        = false;
    previous.marginCells        = 5;
    previous.wellPaddingCells   = 4;
    previous.lgrPaddingCells    = 1;

    RicCreateSectorModelUi ui;
    ui.setSharedSettings( previous );
    const auto settings = ui.settings();

    EXPECT_EQ( previous.makeFluxExecutable, settings.makeFluxExecutable );
    EXPECT_EQ( previous.runMakeFlux, settings.runMakeFlux );
    EXPECT_EQ( previous.marginCells, settings.marginCells );
    EXPECT_EQ( previous.wellPaddingCells, settings.wellPaddingCells );
    EXPECT_EQ( previous.lgrPaddingCells, settings.lgrPaddingCells );
    EXPECT_TRUE( settings.wellName.isEmpty() );
    EXPECT_TRUE( settings.exportFolder.isEmpty() );
    EXPECT_TRUE( settings.deckFile.isEmpty() );
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST( RicCreateSectorModelUiTest, DefaultExportFolderIsNextToTheGrid )
{
    EXPECT_EQ( QString( "/data/model/sector_A_4" ), RicFluxSectorModelTools::defaultExportFolder( "/data/model/CASE.EGRID", "A 4" ) );
}

//--------------------------------------------------------------------------------------------------
/// The context menus refer to the feature by name
//--------------------------------------------------------------------------------------------------
TEST( RicCreateSectorModelUiTest, FeatureIsRegistered )
{
    const auto keys = caf::Factory<caf::CmdFeature, std::string>::instance()->allKeys();
    EXPECT_NE( keys.end(), std::find( keys.begin(), keys.end(), std::string( "RicCreateSectorModelFeature" ) ) );
}
