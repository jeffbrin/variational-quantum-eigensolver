import numpy as np
import matplotlib.pyplot as plt

from gs_energy_finder import estimate_ground_state_energy, Atom

distances = np.arange(0.2, 1.5, 0.005)

gs_energies = np.array(list(map(lambda dist: estimate_ground_state_energy([Atom("H", (0, 0, 0), 1), Atom("H", (0, 0, dist), 1)]), distances)))
print(gs_energies)
min_index, _ = min(enumerate(gs_energies), key=lambda x: x[1])
equilibrium_distance = distances[min_index]
min_gs_energy = gs_energies[min_index]

plt.plot(distances, gs_energies)
ax = plt.subplot(1, 1, 1)
ax.plot(equilibrium_distance, min_gs_energy, "or")
ax.annotate(f"({round(equilibrium_distance, 3)}, {round(min_gs_energy, 3)})", (equilibrium_distance-0.05, min_gs_energy+0.05))
plt.title("$H_2$ Potential Energy vs. Interatomic Distance")
plt.xlabel("Interatomic Distance (Å)")
plt.ylabel("Potential Energy ($E_h$)")
plt.show()