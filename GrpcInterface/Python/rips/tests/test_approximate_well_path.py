import glob
import os
import sys
import tempfile

import pytest

sys.path.insert(1, os.path.join(sys.path[0], "../../"))
import dataroot
import rips
from rips.generated.generated_classes import Perforation

CASE_PATH = dataroot.PATH + "/flux_sector_model/SECTOR_TEST.EGRID"

# Connection factors of well P3 in the restart file at the last time step
P3_CONNECTION_FACTORS = {
    (11, 11, 1): 12.157634,
    (11, 11, 2): 4.863053,
    (11, 11, 3): 30.394085,
}


def _load_case(rips_instance):
    case = rips_instance.project.load_case(path=CASE_PATH)
    if not case.views():
        case.create_view()
    return case


def _compdat(folder):
    rows = {}
    for file_name in glob.glob(os.path.join(folder, "*")):
        keyword = None
        with open(file_name) as f:
            lines = f.readlines()
        for line in lines:
            text = line.split("--")[0].strip()
            if not text:
                continue
            if text in ("WELSPECS", "COMPDAT", "WELSEGS", "COMPSEGS"):
                keyword = text
                continue
            if keyword == "COMPDAT" and not text.startswith("/"):
                tokens = text.replace("'", "").replace("/", "").split()
                rows[(tokens[0], int(tokens[1]), int(tokens[2]), int(tokens[3]))] = (
                    float(tokens[7])
                )
    return rows


def test_create_approximate_well_path(rips_instance, initialize_test):
    case = _load_case(rips_instance)

    well_path = case.create_approximate_well_path(well_name="P3")
    assert well_path.name == "SIM_P3"

    perforations = well_path.descendants(Perforation)
    assert len(perforations) == 3
    for perforation in perforations:
        assert perforation.end_measured_depth > perforation.start_measured_depth
        assert perforation.diameter == pytest.approx(0.2)
        assert perforation.perforated_length_factor > 0.0

    # The connection factors of a completion export are those of the simulation
    with tempfile.TemporaryDirectory(prefix="rips_approx") as folder:
        case.export_completions(well_paths=[well_path], export_folder=folder)
        compdat = _compdat(folder)
        assert len(compdat) == len(P3_CONNECTION_FACTORS)
        for (i, j, k), connection_factor in P3_CONNECTION_FACTORS.items():
            assert compdat[("P3", i, j, k)] == pytest.approx(
                connection_factor, rel=1e-5
            )

    # Creating it again replaces the well path
    case.create_approximate_well_path(well_name="P3")
    names = [w.name for w in rips_instance.project.well_paths()]
    assert names.count("SIM_P3") == 1


def test_create_approximate_well_path_errors(rips_instance, initialize_test):
    case = _load_case(rips_instance)

    with pytest.raises(rips.RipsError, match="no well named"):
        case.create_approximate_well_path(well_name="NO_SUCH_WELL")

    with pytest.raises(rips.RipsError, match="out of range"):
        case.create_approximate_well_path(well_name="P3", time_step=1000)
