import os
import re
import shutil
import sys
import tempfile

import pytest

sys.path.insert(1, os.path.join(sys.path[0], "../../"))
import dataroot
import rips

CASE_PATH = dataroot.PATH + "/flux_sector_model/SECTOR_TEST.EGRID"


def _load_case(rips_instance):
    case = rips_instance.project.load_case(path=CASE_PATH)
    view = case.views()[0] if case.views() else case.create_view()
    return case, view


def test_create_lgr_for_simulation_well(rips_instance, initialize_test):
    case, _ = _load_case(rips_instance)
    assert len(case.grids()) == 1

    case.create_lgr_for_simulation_well(well_name="P3")
    grids = case.grids()
    assert len(grids) == 2
    graded = grids[1].dimensions()
    # Three by three host cells laterally, graded, and two layers per host cell over the three layers
    assert graded.i > 9 and graded.j > 9
    assert graded.k == 6

    # A uniform LGR around another well is added, and the LGR of P3 is kept
    case.create_lgr_for_simulation_well(
        well_name="P1",
        lateral_refinement="UNIFORM",
        refinement_i=3,
        refinement_j=3,
        refinement_k=1,
        lateral_padding=1,
    )
    dims = [
        (g.dimensions().i, g.dimensions().j, g.dimensions().k) for g in case.grids()
    ]
    assert len(dims) == 3
    assert (graded.i, graded.j, graded.k) in dims
    assert (9, 9, 2) in dims

    # Replacing the LGR of P3 keeps the one of P1
    case.create_lgr_for_simulation_well(
        well_name="P3", lateral_refinement="UNIFORM", refinement_k=2
    )
    dims = [
        (g.dimensions().i, g.dimensions().j, g.dimensions().k) for g in case.grids()
    ]
    assert len(dims) == 3
    assert (9, 9, 2) in dims
    assert (15, 15, 6) in dims


def test_create_lgr_for_simulation_well_errors(rips_instance, initialize_test):
    case, _ = _load_case(rips_instance)

    with pytest.raises(rips.RipsError, match="no connections"):
        case.create_lgr_for_simulation_well(well_name="NO_SUCH_WELL")

    with pytest.raises(rips.RipsError, match="growth factor"):
        case.create_lgr_for_simulation_well(well_name="P3", growth_factor=1.0)

    assert len(case.grids()) == 1


def test_create_lgr_holding_other_wells(rips_instance, initialize_test):
    case, _ = _load_case(rips_instance)

    # A padding that reaches the connections of P1 from P3: OPM Flow moves P1's connections
    # into the LGR, so the LGR is made, with a warning
    case.create_lgr_for_simulation_well(
        well_name="P3", lateral_refinement="UNIFORM", lateral_padding=8
    )
    assert len(case.grids()) == 2


def _read_fluxnum(file_name):
    with open(file_name) as f:
        text = f.read()
    body = text.split("FLUXNUM", 1)[1].split("/", 1)[0]
    values = []
    for token in body.split():
        if "*" in token:
            count, value = token.split("*")
            values += [int(value)] * int(count)
        else:
            values.append(int(token))
    return values


def test_create_sector_model_for_well(rips_instance, initialize_test):
    case, view = _load_case(rips_instance)
    case.create_lgr_for_simulation_well(
        well_name="P3", lateral_refinement="UNIFORM", refinement_k=2, lateral_padding=1
    )

    with tempfile.TemporaryDirectory(prefix="rips_sector") as export_folder:
        result = view.create_sector_model_for_well(
            well_name="P3", export_folder=export_folder, run_make_flux=False
        )

        assert result.lgr_name == "LP3"
        assert result.wells[0] == "P3"
        assert result.cell_count > 0
        for file_name in (result.deck_file, result.fluxnum_file):
            assert os.path.exists(file_name)

        fluxnum = _read_fluxnum(result.fluxnum_file)
        assert len(fluxnum) == 15 * 15 * 3
        assert sum(fluxnum) == result.cell_count

        # The LGR and two cells around it are in the sector: host cells I,J 10-12, all layers
        for k in range(3):
            for j in range(7, 14):
                for i in range(7, 14):
                    assert fluxnum[i + 15 * j + 225 * k] == 1

        with open(result.deck_file) as f:
            deck = f.read()
        assert re.search(r"INCLUDE\s+'\S+_GRID\.inc'", deck)
        assert re.search(r"WELSPECL\s+'P3'\s+'G1'\s+'LP3'\s+5\s+5", deck)
        assert re.search(r"COMPDATL\s+'P3'\s+'LP3'\s+5\s+5\s+1\s+6", deck)
        assert re.search(r"^LGR\s+1\s+486", deck, re.MULTILINE)

        grid_include = os.path.join(
            export_folder,
            os.path.basename(result.deck_file)[: -len(".DATA")] + "_GRID.inc",
        )
        with open(grid_include) as f:
            grid_text = f.read()
        for keyword in ("USEFLUX", "FLUXREG", "CARFIN", "ENDFIN"):
            assert keyword in grid_text

    # The sector polygon, used by the polygon filter in the view
    polygons = rips_instance.project.descendants(rips.PolygonCollection)[0].polygons()
    assert [p.name for p in polygons] == ["Sector P3"]
    assert len(polygons[0].coordinates) >= 4


def test_create_sector_model_for_well_needs_deck(rips_instance, initialize_test):
    with tempfile.TemporaryDirectory(prefix="rips_sector") as folder:
        source = os.path.dirname(os.path.abspath(CASE_PATH))
        for name in os.listdir(source):
            if not name.endswith(".DATA"):
                shutil.copy(os.path.join(source, name), folder)
        case = rips_instance.project.load_case(
            path=os.path.join(folder, "SECTOR_TEST.EGRID")
        )
        view = case.create_view()

        with pytest.raises(rips.RipsError, match="does not exist"):
            view.create_sector_model_for_well(
                well_name="P3",
                export_folder=os.path.join(folder, "sector"),
                deck_file=os.path.join(folder, "MISSING.DATA"),
                run_make_flux=False,
            )

        result = view.create_sector_model_for_well(
            well_name="P3",
            export_folder=os.path.join(folder, "sector"),
            deck_file=os.path.join(source, "SECTOR_TEST.DATA"),
            run_make_flux=False,
        )
        assert os.path.exists(result.deck_file)
        assert result.lgr_name == ""


def test_create_well_refinement(rips_instance, initialize_test):
    case, view = _load_case(rips_instance)

    region = view.create_well_refinement(well_name="P3")
    assert region.well_name == "P3"
    values = list(region.boxes)
    assert len(values) % 9 == 0
    boxes = [values[n : n + 9] for n in range(0, len(values), 9)]

    # 9 x 9 at the vertical well, 6 x 6 and 3 x 3 in the rings around it, unrefined vertically
    refined = {}
    for box in boxes:
        hosts = (box[3] - box[0] + 1) * (box[4] - box[1] + 1) * (box[5] - box[2] + 1)
        key = tuple(box[6:9])
        refined[key] = refined.get(key, 0) + hosts
    layers = refined[(9, 9, 1)]
    assert refined == {(9, 9, 1): layers, (6, 6, 1): 8 * layers, (3, 3, 1): 16 * layers}
    assert len(case.grids()) == 1 + len(boxes)

    # Again with other settings: the region and its LGRs are replaced
    region = view.create_well_refinement(
        well_name="P3", lateral_refinement=2, lateral_layers=2, refinement_z=2
    )
    values = list(region.boxes)
    factors = {tuple(values[n + 6 : n + 9]) for n in range(0, len(values), 9)}
    assert factors == {(4, 4, 2), (2, 2, 2)}
    assert len(case.grids()) == 1 + len(values) // 9

    with tempfile.TemporaryDirectory(prefix="rips_sector") as export_folder:
        result = view.create_sector_model_for_well(
            well_name="P3", export_folder=export_folder, run_make_flux=False
        )
        grid_include = os.path.splitext(result.deck_file)[0] + "_GRID.inc"
        with open(grid_include) as f:
            text = f.read()
        assert text.count("CARFIN") == len(values) // 9
        assert "AMALGAM" in text


def test_create_well_refinement_errors(rips_instance, initialize_test):
    case, view = _load_case(rips_instance)

    with pytest.raises(rips.RipsError, match="no well named"):
        view.create_well_refinement(well_name="NO_SUCH_WELL")
    with pytest.raises(rips.RipsError, match="at least 1"):
        view.create_well_refinement(well_name="P3", lateral_layers=0)
