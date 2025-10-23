import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, sqrt, radians
from scipy.interpolate import griddata
import os

from gs_energy_finder import estimate_ground_state_energy, Atom

USE_EXISTING_DATA = True

EQUILIBRIUM_DISTANCE = 0.957
EQUILIBRIUM_DISTANCE_FOCUS_BUFFER = 0.1
EQUILIBRIUM_ANGLE = radians(104.5)
EQUILIBRIUM_ANGLE_FOCUS_BUFFER = 0.2

DATA_FOLDER = "data"
MAX_BOND_LENGTH = 2
MIN_BOND_LENGTH = 0.7
BOND_LENGTH_RINGS = 50
TEST_ANGLES_COUNT = 50
ANGLE_BUFFER = 0.5


if not USE_EXISTING_DATA:

    # Test bond lengths and angles, but focus around known equilibrium point and angle
    bond_lengths = np.linspace(MIN_BOND_LENGTH, EQUILIBRIUM_DISTANCE-EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, int(BOND_LENGTH_RINGS // 3))
    bond_lengths = np.concatenate([bond_lengths, np.linspace(EQUILIBRIUM_DISTANCE-EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, EQUILIBRIUM_DISTANCE+EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, int(BOND_LENGTH_RINGS // 3))])
    bond_lengths = np.concatenate([bond_lengths, np.linspace(EQUILIBRIUM_DISTANCE+EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, MAX_BOND_LENGTH, int(BOND_LENGTH_RINGS // 3))])
    angles = np.linspace(ANGLE_BUFFER, EQUILIBRIUM_ANGLE-EQUILIBRIUM_ANGLE_FOCUS_BUFFER, int(TEST_ANGLES_COUNT // 3))
    angles = np.concatenate([angles, np.linspace(EQUILIBRIUM_ANGLE-EQUILIBRIUM_ANGLE_FOCUS_BUFFER, EQUILIBRIUM_ANGLE+EQUILIBRIUM_ANGLE_FOCUS_BUFFER, int(TEST_ANGLES_COUNT // 3))])
    angles = np.concatenate([angles, np.linspace(EQUILIBRIUM_ANGLE+EQUILIBRIUM_ANGLE_FOCUS_BUFFER, pi, int(TEST_ANGLES_COUNT // 3))])

    hydrogen_coordinates: list[tuple[float, float, float]] = []
    distance_angle_combinations = []
    gs_energies = []
    for length in bond_lengths:
        for angle in angles:
            second_hydrogen_position = (length * sin(angle), length * cos(angle), 0)
            atoms = [
                Atom("O", (0, 0, 0), 8),
                Atom("H", (0, length, 0), 1),
                Atom("H", (second_hydrogen_position[0], second_hydrogen_position[1], second_hydrogen_position[2]), 1)
                ]
            distance_angle_combinations.append((length, angle))
            hydrogen_coordinates.append(second_hydrogen_position)
            gs_energies.append(estimate_ground_state_energy(atoms))

    min_index, _ = min(enumerate(gs_energies), key=lambda x: x[1])
    min_energy_hydrogen_coordinates = hydrogen_coordinates[min_index]
    min_gs_energy = gs_energies[min_index]

    # Plot the surface
    x_plot_data = [c[0] for c in hydrogen_coordinates]
    y_plot_data = [c[1] for c in hydrogen_coordinates]
    z_plot_data = [e.real for e in gs_energies]

else:
    x_plot_data = np.load(os.path.join(DATA_FOLDER, "X_DATA.npy"))
    y_plot_data = np.load(os.path.join(DATA_FOLDER, "Y_DATA.npy"))
    z_plot_data = np.load(os.path.join(DATA_FOLDER, "Z_DATA.npy"))

np.save(os.path.join(DATA_FOLDER, "X_DATA"), np.array(x_plot_data))
np.save(os.path.join(DATA_FOLDER, "Y_DATA"), np.array(y_plot_data))
np.save(os.path.join(DATA_FOLDER, "Z_DATA"), np.array(z_plot_data))

X, Y = np.meshgrid(x_plot_data, y_plot_data)
Z = griddata((x_plot_data, y_plot_data), z_plot_data, (X, Y), method='cubic')

# Get minimum position
min_index, _ = min(enumerate(z_plot_data), key=lambda x: x[1])
min_x = x_plot_data[min_index]
min_y = y_plot_data[min_index]
min_z = z_plot_data[min_index]
min_dist = sqrt(min_x **2 + min_y **2)

# Plot the surface
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# ax.scatter(x_plot_data, y_plot_data, z_plot_data)
# ax.plot([min_x], [min_y], [min_z], marker="H", c="r", markersize=4)
# ax.plot([0], [0], [min_z], marker="o", c="g", markersize=4)
# ax.plot([0], [sqrt(min_x **2 + min_y **2)], [min_z], marker="o", c="r", markersize=4)
ax.plot_trisurf(x_plot_data, y_plot_data, z_plot_data)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_xlim(-MIN_BOND_LENGTH, MAX_BOND_LENGTH)
ax.set_ylim(-MAX_BOND_LENGTH, MAX_BOND_LENGTH)
ax.set_zlim(min_z, min_z+0.25)

# Plot the H2O configuration as another plot

plt.show()