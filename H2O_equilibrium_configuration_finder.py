import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos

from gs_energy_finder import estimate_ground_state_energy

bond_lengths = np.arange(0.2, 1.5, 0.005)
angles = np.arange(0, pi, 0.01)
molecule_template = "H {0} {1} {2}; H {3} {4} {5}; O 0 0 0"

hydrogen_coordinates: list[tuple[float, float, float, float, float, float]] = []
distance_angle_combinations = []
gs_energies = []
for length in bond_lengths:
    for angle in angles:
        distance_angle_combinations.append((length, angle))
        positions = (0, length, 0, -length * sin(angle), length * cos(angle), 0)
        hydrogen_coordinates.append(positions)
        molecule = molecule_template.format(positions)
        gs_energies.append(estimate_ground_state_energy(molecule))

gs_energies = np.array(list(map(lambda dist: estimate_ground_state_energy(molecule_template.format(dist), dist), distances)))
print(gs_energies)
min_index, _ = min(enumerate(gs_energies), key=lambda x: x[1])
equilibrium_distance = bond_lengths[min_index]
min_gs_energy = gs_energies[min_index]

plt.plot(bond_lengths, gs_energies)
ax = plt.subplot(1, 1, 1)
ax.plot(equilibrium_distance, min_gs_energy, "or")
ax.annotate(f"({round(equilibrium_distance, 3)}, {round(min_gs_energy, 3)})", (equilibrium_distance-0.05, min_gs_energy+0.05))
plt.title("Potential Energy vs. Interatomic Distance")
plt.xlabel("Interatomic Distance (Å)")
plt.ylabel("Potential Energy ($E_h$)")
plt.show()