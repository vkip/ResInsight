###################################################################################
# This example creates an approximate well path for a simulation well, from the
# connections in the simulation results, and exports its completions:
# 1. The well path 'SIM_<well>' runs through the centres of the connected cells
# 2. Each connected cell gets a perforation with the connection's diameter and skin
#    factor, and a perforated length factor that makes the exported connection
#    factors equal to those of the simulation
# 3. The completion export (WELSPECS/COMPDAT) uses the default options
###################################################################################

import os
import tempfile

import rips

resinsight = rips.Instance.find()

case_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "../../../../../TestModels/flux_sector_model/SECTOR_TEST.EGRID",
)
case = resinsight.project.load_case(path=os.path.abspath(case_path))
case.create_view()

# The connections of the last time step where the well has any
well_path = case.create_approximate_well_path(well_name="P3")
print(f"Well path: {well_path.name}")

export_folder = tempfile.mkdtemp(prefix="approximate_well_path_")
case.export_completions(well_paths=[well_path], export_folder=export_folder)
print(f"Completions exported to {export_folder}")
