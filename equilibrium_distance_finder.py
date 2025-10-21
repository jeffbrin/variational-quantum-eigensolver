from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.mappers import JordanWignerMapper
import numpy as np
import matplotlib.pyplot as plt

ANGSTROM_TO_BOHR = 1.8897259886

def estimate_ground_state_energy(molecule: str, distance: float) -> float:
    """Estimates the ground state energy for a particular molecule using VQE

    Parameters
    ----------
    molecule : str
        The molecule to estimate ground state energy for.
    distance : float
        The distance between the two atoms in the molecule in Angstrom units.

    Returns
    -------
    float
        The estimated ground state energy in Hartrees (electronic and nuclear repulsion energy).
    """

    print(f"Calculating for {molecule}")

    # Define the PySCF driver to simulate an H_2 molecule 
    # with the hydrogen atoms distanced 0.735 ANGSTROM units apart
    driver = PySCFDriver(
        atom=molecule,
        basis="sto3g",
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
    )

    problem = driver.run()

    solver = GroundStateEigensolver(
        JordanWignerMapper(),
        NumPyMinimumEigensolver(),
    )

    result = solver.solve(problem)
    electronic_energy = result.eigenvalues[0]

    # The Qiskit ground state eigensolver doesn't automatically dd nuclear repulsion energy.
    # Manually add it to the result
    dist_bohr = distance * ANGSTROM_TO_BOHR
    nuclear_repulsion_energy = 1/dist_bohr

    return electronic_energy + nuclear_repulsion_energy

EQUILIBRIUM = 0.735
distances = np.arange(0.2, 1.5, 0.1)
molecule_template = "H 0 0 0; H 0 0 {0}"

gs_energies = np.array(list(map(lambda dist: estimate_ground_state_energy(molecule_template.format(dist), dist), distances)))
print(gs_energies)
min_index, _ = min(enumerate(gs_energies), key=lambda x: x[1])
equilibrium_distance = distances[min_index]
min_gs_energy = gs_energies[min_index]

plt.plot(distances, gs_energies)
plt.axvline(x = equilibrium_distance, color = 'g', label = 'Equilibrium Distance')
ax = plt.subplot(1, 1, 1)
ax.plot(equilibrium_distance, min_gs_energy, "or")
plt.show()