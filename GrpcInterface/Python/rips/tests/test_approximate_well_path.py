import glob
import os
import sys
import tempfile

import pytest

sys.path.insert(1, os.path.join(sys.path[0], "../../"))
import dataroot
import rips
from rips.generated.generated_classes import Perforation, WellPathValve

CASE_PATH = dataroot.PATH + "/flux_sector_model/SECTOR_TEST.EGRID"
MSW_CASE_PATH = (
    dataroot.PATH
    + "/../ApplicationLibCode/UnitTests/TestData/RifRestartWellConnectionReader/MSW_ICD.EGRID"
)
MSW_DEVICES_CASE_PATH = (
    dataroot.PATH
    + "/../ApplicationLibCode/UnitTests/TestData/RifRestartWellConnectionReader/MSW_DEVICES.EGRID"
)

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


def _keyword_rows(folder, keyword):
    rows = []
    for file_name in glob.glob(os.path.join(folder, "*")):
        current = None
        with open(file_name) as f:
            for line in f:
                text = line.split("--")[0].strip()
                if not text:
                    continue
                if text.isupper() and text.isalpha():
                    current = text
                    continue
                if current == keyword and not text.startswith("/"):
                    rows.append(text.replace("'", "").replace("/", "").split())
    return rows


def _reconstructed_valves(well_path):
    return [valve for valve in well_path.descendants(WellPathValve) if valve.name]


def test_create_approximate_well_path_multi_segment(rips_instance, initialize_test):
    # A horizontal multi-segment well with two valves (WSEGVALV), the last two connections on valve segments
    case = rips_instance.project.load_case(path=MSW_CASE_PATH)
    case.create_view()

    well_path = case.create_approximate_well_path(well_name="P1")
    perforations = sorted(
        well_path.descendants(Perforation), key=lambda p: p.start_measured_depth
    )
    assert len(perforations) == 3
    # The measured depths are those along the tubing (COMPSEGS)
    for perforation, (start, end) in zip(
        perforations, ((2100, 2200), (2200, 2300), (2300, 2400))
    ):
        assert perforation.start_measured_depth == pytest.approx(start, abs=0.5)
        assert perforation.end_measured_depth == pytest.approx(end, abs=0.5)
    # A valve over each perforation of the two connections on valve segments
    assert len(_reconstructed_valves(well_path)) >= 2

    with tempfile.TemporaryDirectory(prefix="rips_approx_msw") as folder:
        case.export_completions(
            well_paths=[well_path], export_folder=folder, include_msw=True
        )
        # The top segment (WELSEGS) and the valves (WSEGVALV) of the simulation
        top = [row for row in _keyword_rows(folder, "WELSEGS") if row[0] == "P1"]
        assert len(top) == 1
        assert float(top[0][1]) == pytest.approx(2002.5, abs=0.01)
        assert float(top[0][2]) == pytest.approx(2050.0, abs=0.01)

        valves = _keyword_rows(folder, "WSEGVALV")
        assert sorted(round(float(row[2]), 3) for row in valves) == [0.7, 0.8]
        assert sum(float(row[3]) for row in valves) == pytest.approx(5.0e-4, rel=1e-3)

    # Without the segment data: no valves
    well_path = case.create_approximate_well_path(
        well_name="P1", use_multi_segment_data=False
    )
    assert len(_reconstructed_valves(well_path)) == 0


def test_export_valves_sharing_a_cell(rips_instance, initialize_test):
    # A valve that reaches into the cell of the next perforation, which has its own valve: the cell gets one ICD
    # segment with both shares of the area, as a cell connects to one segment only
    case = rips_instance.project.load_case(path=MSW_CASE_PATH)
    case.create_view()
    well_path = case.create_approximate_well_path(well_name="P1")
    first = min(
        well_path.descendants(Perforation), key=lambda p: p.start_measured_depth
    )
    template = rips_instance.project.valve_templates().add_template(
        completion_type="ICD", orifice_diameter=8.0, flow_coefficient=0.7
    )
    first.add_valve(start_md=2150.0, end_md=2250.0, valve_count=1, template=template)

    with tempfile.TemporaryDirectory(prefix="rips_shared_cell") as folder:
        case.export_completions(
            well_paths=[well_path], export_folder=folder, include_msw=True
        )
        cells = [tuple(row[:3]) for row in _keyword_rows(folder, "COMPSEGS")[1:]]
        assert len(cells) == len(set(cells))
        area = 3.14159265 * 0.004**2
        valves = _keyword_rows(folder, "WSEGVALV")
        assert sum(float(row[3]) for row in valves) == pytest.approx(
            5.0e-4 + area, rel=1e-3
        )


def test_create_approximate_well_path_lateral_and_devices(
    rips_instance, initialize_test
):
    # A multi-segment well with a lateral branch, two SICDs and two AICDs (WSEGSICD, WSEGAICD)
    case = rips_instance.project.load_case(path=MSW_DEVICES_CASE_PATH)
    case.create_view()

    well_path = case.create_approximate_well_path(well_name="P1")
    lateral = rips_instance.project.well_path_by_name("SIM_P1 L2")
    assert lateral is not None
    assert len(well_path.descendants(Perforation)) == 3
    lateral_perforations = sorted(
        lateral.descendants(Perforation), key=lambda p: p.start_measured_depth
    )
    assert len(lateral_perforations) == 3
    # The measured depths of the lateral are those along its tubing (COMPSEGS)
    for perforation, (start, end) in zip(
        lateral_perforations, ((2300, 2400), (2400, 2500), (2500, 2600))
    ):
        assert perforation.start_measured_depth == pytest.approx(start, abs=0.5)
        assert perforation.end_measured_depth == pytest.approx(end, abs=0.5)

    with tempfile.TemporaryDirectory(prefix="rips_approx_devices") as folder:
        case.export_completions(
            well_paths=[well_path], export_folder=folder, include_msw=True
        )
        # The lateral is a branch of the well, with its own tubing diameter
        segments = [row for row in _keyword_rows(folder, "WELSEGS") if row[0] != "P1"]
        assert {row[2] for row in segments} >= {"1", "2"}
        assert {round(float(row[6]), 3) for row in segments if row[2] == "2"} == {0.12}

        # The strengths of the devices include their flow scaling factors
        sicds = _keyword_rows(folder, "WSEGSICD")
        assert sorted(float(row[3]) for row in sicds) == pytest.approx(
            [3.0e-5 * 0.5**2, 2.0e-5 * 1.5**2], rel=1e-4
        )
        aicds = _keyword_rows(folder, "WSEGAICD")
        assert sorted(float(row[3]) for row in aicds) == pytest.approx(
            [2.0e-6 * 0.12**2.5, 2.0e-8 * 2.0**2.5], rel=1e-4
        )
        assert all(float(row[12]) == pytest.approx(2.5) for row in aicds)

    # Recreating replaces the lateral
    case.create_approximate_well_path(well_name="P1")
    names = [w.name for w in rips_instance.project.well_paths()]
    assert names.count("SIM_P1 L2") == 1


def test_create_approximate_well_path_errors(rips_instance, initialize_test):
    case = _load_case(rips_instance)

    with pytest.raises(rips.RipsError, match="no well named"):
        case.create_approximate_well_path(well_name="NO_SUCH_WELL")

    with pytest.raises(rips.RipsError, match="out of range"):
        case.create_approximate_well_path(well_name="P3", time_step=1000)


def test_select_in_project_tree(rips_instance, initialize_test):
    case = _load_case(rips_instance)

    # Without a main window (console mode) selecting does nothing, but must succeed
    simulation_well = [w for w in case.simulation_wells() if w.name == "P3"][0]
    simulation_well.select()
    case.create_approximate_well_path(well_name="P3").select()
