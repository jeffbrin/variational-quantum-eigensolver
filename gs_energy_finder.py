from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.mappers import JordanWignerMapper

from math import sqrt

ANGSTROM_TO_BOHR = 1.8897259886

class Atom:
    def __init__(self, symbol: str, position: tuple[float, float, float], atomic_number: int):
        self.symbol = symbol
        self.position = position
        self.atomic_number = atomic_number

    def to_scf_string(self) -> str:
        return f"{self.symbol} {self.position[0]} {self.position[1]} {self.position[2]}"

def estimate_ground_state_energy(atoms: list[Atom]) -> float:
    """Estimates the ground state energy for a particular molecule using VQE

    Parameters
    ----------
    molecule : list[Atom]
        The atoms which make up the molecule

    Returns
    -------
    float
        The estimated ground state energy in Hartrees (electronic and nuclear repulsion energy).
    """


    molecule = "; ".join(a.to_scf_string() for a in atoms)
    print(f"Estimating ground state energy for {molecule}")

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
    # Calculate nuclear repulsion for each pair of atoms
    nuclear_repulsion_energy = 0
    for first_atom_index in range(len(atoms)):
        for second_atom_index in range(first_atom_index+1, len(atoms)):
            atom1, atom2 = atoms[first_atom_index], atoms[second_atom_index]
            distance = sqrt((atom1.position[0] - atom2.position[0]) ** 2 + 
                            (atom1.position[1] - atom2.position[1]) ** 2 +
                            (atom1.position[2] - atom2.position[2]) ** 2)
            dist_in_bohr = distance * ANGSTROM_TO_BOHR
            nuclear_repulsion_energy += atom1.atomic_number * atom2.atomic_number / dist_in_bohr

    return electronic_energy + nuclear_repulsion_energy
