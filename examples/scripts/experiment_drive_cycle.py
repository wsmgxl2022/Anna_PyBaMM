#
# Constant-current constant-voltage charge with US06 Drive Cycle using Experiment Class.
#
import pybamm
import pandas as pd
import os

os.chdir(pybamm.__path__[0] + "/..")

pybamm.set_logging_level("INFO")

# import drive cycle from file
drive_cycle_current = pd.read_csv(
    "pybamm/input/drive_cycles/US06.csv", comment="#", header=None
).to_numpy()

print(drive_cycle_current)


# # Map Drive Cycle
# def map_drive_cycle(x, min_op_value, max_op_value):
#     min_ip_value = x[:, 1].min()
#     max_ip_value = x[:, 1].max()
#     x[:, 1] = (x[:, 1] - min_ip_value) / (max_ip_value - min_ip_value) * (
#         max_op_value - min_op_value
#     ) + min_op_value
#     return x


# # Map current drive cycle to voltage and power
# drive_cycle_power = map_drive_cycle(drive_cycle_current, 1.5, 3.5)

# experiment = pybamm.Experiment(
#     [
#         "Charge at 1 A until 4.0 V",
#         "Hold at 4.0 V until 50 mA",
#         "Rest for 30 minutes",
#         "Run US06_A (A)",
#         "Rest for 30 minutes",
#         "Run US06_W (W)",
#         "Rest for 30 minutes",
#     ],
#     drive_cycles={
#         "US06_A": drive_cycle_current,
#         "US06_W": drive_cycle_power,
#     },
# )

# model = pybamm.lithium_ion.DFN()
# sim = pybamm.Simulation(model, experiment=experiment, solver=pybamm.CasadiSolver())
# sim.solve()

# # Show all plots
# sim.plot()
