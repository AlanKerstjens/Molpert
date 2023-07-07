import time
import molpert as mpt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from typing import Union
from numpy.random import rand, uniform

OptionalConstraints = Union[mpt.MolecularConstraints, None]

def RefreshConstraints(
    constraints: OptionalConstraints,
    molecule: Chem.Mol) -> OptionalConstraints:
    mpt.TagMolecule(molecule, skip_if_tagged=True)
    if constraints is not None:
        constraints.ClearConstraints()
        constraints.GenerateConstraints(molecule)
    return constraints


class MoleculePerturber:
    def __init__(self,
        constraints: OptionalConstraints = None,
        start_molecule_id: int = 0):
        self.perturber = mpt.MoleculePerturber()
        self.prng = mpt.PRNG()
        self.constraints = constraints
        self.max_molecule_id = start_molecule_id
        self.n_molecules = 0
        self.design_time = 0.0

    def __call__(self, molecule: Chem.Mol) -> Chem.Mol:
        t0 = time.perf_counter()
        RefreshConstraints(self.constraints, molecule)
        perturbation = self.perturber.Perturbation(
            molecule, self.prng, self.constraints)
        if not perturbation:
            self.design_time += time.perf_counter() - t0
            return molecule
        perturbed_molecule = Chem.Mol(molecule)
        mpt.ExecutePerturbation(perturbation, perturbed_molecule)
        self.design_time += time.perf_counter() - t0
        self.n_molecules += 1
        self.max_molecule_id += 1
        perturbed_molecule.SetIntProp("ID", self.max_molecule_id)
        return perturbed_molecule

    def RandomMolecule(self,
        seed: Chem.Mol, n_atoms: int, n_bonds: int) -> Chem.Mol:
        molecule = Chem.Mol(seed)
        for _ in range(n_atoms - molecule.GetNumAtoms()):
            RefreshConstraints(self.constraints, molecule)
            perturbation = self.perturber.InsertAtom(
                molecule, self.prng, self.constraints)
            if perturbation:
                mpt.ExecutePerturbation(perturbation, molecule)
        for _ in range(n_bonds - molecule.GetNumBonds()):
            RefreshConstraints(self.constraints, molecule)
            perturbation = self.perturber.InsertBond(
                molecule, self.prng, self.constraints)
            if perturbation:
                mpt.ExecutePerturbation(perturbation, molecule)
        self.max_molecule_id += 1
        molecule.SetIntProp("ID", self.max_molecule_id)
        return molecule


class MoleculeEnumerator:
    def __init__(self,
        constraints: OptionalConstraints = None):
        self.perturber = mpt.MoleculePerturber()
        self.perturber.assess_connectivity_with_sssr = True
        self.perturber.assess_distances_with_distance_matrix = True
        self.queue = mpt.MolecularPerturbationQueue()
        self.constraints = constraints
        self.n_molecules = 0
        self.design_time = 0.0

    def __call__(self, molecule: Chem.Mol) -> list[Chem.Mol]:
        neighbors = []
        neighbor_hashes = set()
        t0 = time.perf_counter()
        RefreshConstraints(self.constraints, molecule)
        self.perturber.Perturbations(self.queue, molecule, self.constraints)
        self.design_time += time.perf_counter() - t0
        while not self.queue.empty():
            perturbation = self.queue.front()
            perturbed_molecule = Chem.Mol(molecule)
            t0 = time.perf_counter()
            mpt.ExecutePerturbation(perturbation, perturbed_molecule)
            self.design_time += time.perf_counter() - t0
            self.n_molecules += 1
            self.queue.pop()
            perturbed_hash = mpt.HashMolecule(perturbed_molecule)
            if perturbed_hash in neighbor_hashes:
                continue
            neighbors.append(perturbed_molecule)
            neighbor_hashes.add(perturbed_hash)
        return neighbors


class Mimicry:
    def __init__(self,
        constraints: OptionalConstraints = None,
        min_fraction_atoms: float = 0.1,
        max_fraction_atoms: float = 0.5,
        start_molecule_id: int = 0):
        self.constraints = constraints
        self.min_fraction_atoms = min_fraction_atoms
        self.max_fraction_atoms = max_fraction_atoms
        self.requests = mpt.SubgraphRequests()
        self.prng = mpt.PRNG()
        self.max_molecule_id = start_molecule_id
        self.n_molecules = 0
        self.design_time = 0.0

    def __call__(self, molecule1: Chem.Mol, molecule2: Chem.Mol) -> Chem.Mol:
        mpt.TagMolecule(molecule1, skip_if_tagged=True)
        mpt.TagMolecule(molecule2, skip_if_tagged=True)
        parent1 = molecule1
        parent2 = molecule2
        if rand() > 0.5:
            parent1, parent2 = parent2, parent1
        subgraph_size = int(
            uniform(self.min_fraction_atoms, self.max_fraction_atoms) * \
            parent1.GetNumAtoms())
        self.requests.size = max(1, subgraph_size) # Mimic at least 1 atom.
        t0 = time.perf_counter()
        perturbation = mpt.Mimicry(parent1, parent2, self.requests, self.prng)
        self.design_time += time.perf_counter() - t0
        if not perturbation:
            return parent1
        RefreshConstraints(self.constraints, molecule1)
        if self.constraints:
            if not self.constraints.IsPerturbationAllowed(perturbation, parent1):
                return parent1
        t0 = time.perf_counter()
        child = Chem.Mol(parent1)
        mpt.ExecutePerturbation(perturbation, child)
        self.design_time += time.perf_counter() - t0
        self.n_molecules += 1
        self.max_molecule_id += 1
        child.SetIntProp("ID", self.max_molecule_id)
        return child


def AtomKeyFingerprint(
    molecule: Chem.Mol,
    radius: int = 2) -> DataStructs.UIntSparseIntVect:
    max_unsigned = 2**32 - 1 # We return uint64, but the RDKit expects uint32.
    invariants = [mpt.AtomKeyHash(a) % max_unsigned for a in molecule.GetAtoms()]
    return AllChem.GetMorganFingerprint(molecule, radius, invariants=invariants)

class MoleculeSimilarity:
    def __init__(self, cache_cleaning_frequency: int = 100000):
        self.cache: dict[Chem.Mol, list[DataStructs.UIntSparseIntVect, int]] = {}
        self.n_calls = 0
        self.cache_hits = 0
        self.cache_misses = 0
        self.cache_cleaning_frequency = cache_cleaning_frequency

    def get(self, molecule: Chem.Mol) -> DataStructs.UIntSparseIntVect:
        molecule_id = molecule.GetIntProp("ID")
        try:
            fp = self.cache[molecule_id][0]
            self.cache[molecule_id][1] += 1
            self.cache_hits += 1
        except KeyError:
            fp = AtomKeyFingerprint(molecule)
            self.cache[molecule_id] = [fp, 1]
            self.cache_misses += 1
        return fp

    def clean(self):
        self.cache = {k:v for k,v in self.cache.items() if not v[1]}

    def clear(self):
        self.cache = {}

    def __call__(self, molecule1: Chem.Mol, molecule2: Chem.Mol) -> float:
        fp1 = self.get(molecule1)
        fp2 = self.get(molecule2)
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        self.n_calls += 1
        if not self.n_calls % self.cache_cleaning_frequency:
            self.clean()
        return similarity