###################################################################################
# This example creates a local grid refinement (LGR) around a simulation well and a
# sector model around the same well for an OPM Flow run with USEFLUX:
# 1. Load a case whose deck lies next to the grid file
# 2. Refine around the well: twice vertically and logarithmically laterally, with
#    cells of about 1 m in the connection cells
# 3. Create the sector model. The sector is shown as a polygon filter, written as a
#    FLUXNUM region, make_flux builds its boundary data from the case's results,
#    and the deck is copied with USEFLUX, FLUXREG, the LGR and the well in the LGR
#
# The parent run needs FLUXALL in its SUMMARY section, and preferably a restart at
# every report step (RPTRST BASIC=2), see the sector modelling documentation of OPM
# Flow. Graded LGRs (NXFIN/HXFIN) need OPM Flow's Conforming LGR backend. A well with
# connections far apart, e.g. in two zones, gets one LGR per part.
###################################################################################

import os

import rips

resinsight = rips.Instance.find()

case_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "../../../../../TestModels/flux_sector_model/SECTOR_TEST.EGRID",
)
case = resinsight.project.load_case(path=os.path.abspath(case_path))
view = case.create_view()
well_name = "P3"

# Logarithmic lateral refinement, 1 m at the well, growing by about 1.5 per cell
case.create_lgr_for_simulation_well(
    well_name=well_name,
    lateral_refinement="LOGARITHMIC",
    min_cell_size=1.0,
    growth_factor=1.5,
    refinement_k=2,
    lateral_padding=1,
)
for grid in case.grids():
    dims = grid.dimensions()
    print(f"Grid {grid.index}: {dims.i} x {dims.j} x {dims.k}")

# make_flux is taken from PATH unless given with make_flux_executable
result = view.create_sector_model_for_well(well_name=well_name, run_make_flux=False)
print(f"Sector deck: {result.deck_file}")
print(f"FLUXNUM:     {result.fluxnum_file}")
print(f"LGRs:        {result.lgr_name}")
print(f"Cells:       {result.cell_count} ({result.active_cell_count} active)")
print(f"Wells:       {', '.join(result.wells)}")
for warning in result.warnings:
    print(f"Note: {warning}")
