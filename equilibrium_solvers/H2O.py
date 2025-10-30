import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, sqrt, radians, acos, degrees
import os

from gs_energy_finder import estimate_ground_state_energy, Atom

USE_EXISTING_DATA = False
CONCENTRATE_AROUND_KNOWN_EQUILIBRIUM = False

EQUILIBRIUM_DISTANCE = 0.957
EQUILIBRIUM_DISTANCE_FOCUS_BUFFER = 0.1
EQUILIBRIUM_ANGLE = radians(104.5)
EQUILIBRIUM_ANGLE_FOCUS_BUFFER = 0.2

DATA_FOLDER = "data"
MAX_BOND_LENGTH = 2
MIN_BOND_LENGTH = 0.1
BOND_LENGTH_RINGS = 5
TEST_ANGLES_COUNT = 20
ANGLE_BUFFER = radians(60)

if not USE_EXISTING_DATA:

    # Test bond lengths and angles, but focus around known equilibrium point and angle
    if CONCENTRATE_AROUND_KNOWN_EQUILIBRIUM:
        bond_lengths = np.linspace(MIN_BOND_LENGTH, EQUILIBRIUM_DISTANCE-EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, int(BOND_LENGTH_RINGS // 3))
        bond_lengths = np.concatenate([bond_lengths, np.linspace(EQUILIBRIUM_DISTANCE-EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, EQUILIBRIUM_DISTANCE+EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, int(BOND_LENGTH_RINGS // 3))])
        bond_lengths = np.concatenate([bond_lengths, np.linspace(EQUILIBRIUM_DISTANCE+EQUILIBRIUM_DISTANCE_FOCUS_BUFFER, MAX_BOND_LENGTH, int(BOND_LENGTH_RINGS // 3))])
        angles = np.linspace(ANGLE_BUFFER, EQUILIBRIUM_ANGLE-EQUILIBRIUM_ANGLE_FOCUS_BUFFER, int(TEST_ANGLES_COUNT // 3))
        angles = np.concatenate([angles, np.linspace(EQUILIBRIUM_ANGLE-EQUILIBRIUM_ANGLE_FOCUS_BUFFER, EQUILIBRIUM_ANGLE+EQUILIBRIUM_ANGLE_FOCUS_BUFFER, int(TEST_ANGLES_COUNT // 3))])
        angles = np.concatenate([angles, np.linspace(EQUILIBRIUM_ANGLE+EQUILIBRIUM_ANGLE_FOCUS_BUFFER, pi, int(TEST_ANGLES_COUNT // 3))])
    else:
        bond_lengths = np.linspace(MIN_BOND_LENGTH, MAX_BOND_LENGTH, BOND_LENGTH_RINGS)
        angles = np.linspace(ANGLE_BUFFER, pi-ANGLE_BUFFER, TEST_ANGLES_COUNT)

    print([degrees(a) for a in angles])

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

# X, Y = np.meshgrid(x_plot_data, y_plot_data)
# Z = griddata((x_plot_data, y_plot_data), z_plot_data, (X, Y), method='cubic')

# Get minimum position
min_index, _ = min(enumerate(z_plot_data), key=lambda x: x[1])
min_x = x_plot_data[min_index]
min_y = y_plot_data[min_index]
min_z = z_plot_data[min_index]
min_dist = sqrt(min_x **2 + min_y **2)
min_angle = acos(min_y/min_dist) * 180 / pi

print("Min Angle:", min_angle)
print("Min Distance H-O:", min_dist)

# Plot the surface
fig = plt.figure()

# plot_3d
ax = fig.add_subplot(projection="3d")
ax.scatter(x_plot_data, y_plot_data, z_plot_data, s=5)
# ax.plot_trisurf(x_plot_data, y_plot_data, z_plot_data)

ax.set_xlim(-MIN_BOND_LENGTH, MAX_BOND_LENGTH)
ax.set_ylim(-MAX_BOND_LENGTH, MAX_BOND_LENGTH)
ax.set_zlim(min_z, min_z+0.25)
ax.plot([min_x], [min_y], [min_z], marker="o", c="r", markersize=10)
ax.text(min_x, min_y, min_z, f"({round(min_x, 2)}, {round(min_y, 2)}, {round(min_z, 2)})")
ax.set_title("Ground State Energy (E$_h$) vs. Hydrogen Position")
ax.set_xlabel("Second Hydrogen x Position (Å)")
ax.set_ylabel("Second Hydrogen y Position (Å)")
ax.set_zlabel("Ground State Energy (E$_h$)")

# Plot the H2O configuration as another plot
# Plot configuration
# ax2 = fig.add_subplot()
# plt.plot([min_x, 0], [min_y, 0], 'ro-')
# plt.plot([0, 0], [0, min_dist], 'ro-')
# # ax2.plot([min_x], [min_y], marker="o", c="r", markersize=4)
# ax2.plot([0], [0], marker="o", c="g")
# arc = Arc((0,0), 0.25, 0.25, angle=90-min_angle, theta1=0, theta2=min_angle)
# ax2.add_patch(arc)
# ax2.text(0.1, 0.1, f"{round(min_angle, 2)}$\degree$")
# ax2.text(0, -0.06, f"O (0, 0)")
# ax2.text(min_x-0.2, min_y+0.05, f"H ({round(min_x, 2)}, {round(min_y, 2)})")
# ax2.text(0.025, min_dist, f"H (0, {round(min_dist, 2)})")

# ax2.set_xlabel("x Position Relative to Oxygen (Å)")
# ax2.set_ylabel("y Position Relative to Oxygen (Å)")
# ax2.set_title("H$_2$O Equilibrium Configuration Approximation")
# ax2.plot([0], [sqrt(min_x **2 + min_y **2)], marker="o", c="r", markersize=4)

plt.show()