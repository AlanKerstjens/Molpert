import molpert as mpt
from rdkit import Chem
from typing import Callable, Union


Key = tuple[int, ...]
AtomBond = Union[Chem.Atom, Chem.Bond]
AtomKeyGenerator = Callable[[Chem.Atom], Key]
BondKeyGenerator = Callable[[Chem.Bond], Key]
EnvironmentKeyGenerator = Callable[[Chem.Atom], Key]
KeyGenerator = Callable[[AtomBond], Key]
AtomHasher = Callable[[Chem.Atom], int]


def ValenceConstraint(atom: Chem.Atom) -> bool:
    valence = mpt.ExplicitValence(atom)
    return mpt.IsValenceBelowMax(atom.GetAtomicNum(), int(valence))

def ValenceConstraintGenerator(prior_atom: Chem.Atom, posterior_atom: Chem.Atom):
    return ValenceConstraint


def LocalAtomKey(atom: Chem.Atom) -> Key:
    return (
        atom.GetDegree(),
        mpt.ExplicitValence(atom), 
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        atom.GetNumExplicitHs())

def RingAwareAtomKey(atom: Chem.Atom) -> Key:
    molecule = atom.GetOwningMol()
    mpt.FindSSSRIfNotInitialized(molecule)
    ring_info = molecule.GetRingInfo()
    atom_idx = atom.GetIdx()
    n_rings = ring_info.NumAtomRings(atom_idx)
    min_atom_ring_size = 0
    max_atom_ring_size = 0
    if n_rings:
        atom_ring_sizes = ring_info.AtomRingSizes(atom_idx)
        min_atom_ring_size = min(atom_ring_sizes)
        max_atom_ring_size = max(atom_ring_sizes)
    return (
        n_rings,
        min_atom_ring_size,
        max_atom_ring_size,
        int(atom.GetIsAromatic()),
        atom.GetDegree(),
        mpt.ExplicitValence(atom),
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        atom.GetNumExplicitHs())


def LocalBondKey(bond: Chem.Bond) -> Key:
    begin_atom_key = LocalAtomKey(bond.GetBeginAtom())
    end_atom_key = LocalAtomKey(bond.GetEndAtom())
    if end_atom_key < begin_atom_key:
        begin_atom_key, end_atom_key = end_atom_key, begin_atom_key
    return begin_atom_key + end_atom_key + (bond.GetBondTypeAsDouble(),)

def RingAwareBondKey(bond: Chem.Bond) -> Key:
    begin_atom_key = RingAwareAtomKey(bond.GetBeginAtom())
    end_atom_key = RingAwareAtomKey(bond.GetEndAtom())
    if end_atom_key < begin_atom_key:
        begin_atom_key, end_atom_key = end_atom_key, begin_atom_key
    return begin_atom_key + end_atom_key + (bond.GetBondTypeAsDouble(),) 


def HashAtomKey(atom: Chem.Atom, atom_key_generator: AtomKeyGenerator) -> int:
    atom_key = atom_key_generator(atom)
    atom_hash = hash(atom_key) & 0xFFFFFFFFFFFFFFFF # Convert int64 to uint64
    return atom_hash


class KeyConstraint:
    def __init__(self,
        key_generator: KeyGenerator,
        key_frequency_dictionary: dict[Key, int],
        key_frequency_threshold: int = 1):
        self.key_generator = key_generator
        self.key_frequency_dictionary = key_frequency_dictionary
        self.key_frequency_threshold = key_frequency_threshold

    def __call__(self, atom_bond: AtomBond) -> bool:
        key = self.key_generator(atom_bond)
        key_frequency = self.key_frequency_dictionary.get(key, 0)
        return key_frequency >= self.key_frequency_threshold


class AtomKeyDictionary:
    def __init__(self, 
        atom_key_generator: AtomKeyGenerator,
        atom_key_frequency_threshold: int = 1):
        self.atom_key_generator = atom_key_generator
        self.atom_key_frequency_dictionary = {}
        self.atom_key_frequency_threshold = atom_key_frequency_threshold
    
    def __call__(self,
        prior_atom: Chem.Atom,
        posterior_atom: Chem.Atom) -> Union[KeyConstraint, None]:
        if not posterior_atom:
            return None
        return KeyConstraint(
            key_generator=self.atom_key_generator,
            key_frequency_dictionary=self.atom_key_frequency_dictionary,
            key_frequency_threshold=self.atom_key_frequency_threshold)

    def AddMolecule(self, molecule: Chem.Mol):
        for atom in molecule.GetAtoms():
            atom_key = self.atom_key_generator(atom)
            try:
                self.atom_key_frequency_dictionary[atom_key] += 1
            except KeyError:
                self.atom_key_frequency_dictionary[atom_key] = 1

class BondKeyDictionary:
    def __init__(self, 
        bond_key_generator: BondKeyGenerator,
        bond_key_frequency_threshold: int = 1):
        self.bond_key_generator = bond_key_generator
        self.bond_key_frequency_dictionary = {}
        self.bond_key_frequency_threshold = bond_key_frequency_threshold
    
    def __call__(self,
        prior_bond: Chem.Bond,
        posterior_bond: Chem.Bond) -> Union[KeyConstraint, None]:
        if not posterior_bond:
            return None
        return KeyConstraint(
            key_generator=self.bond_key_generator,
            key_frequency_dictionary=self.bond_key_frequency_dictionary,
            key_frequency_threshold=self.bond_key_frequency_threshold)

    def AddMolecule(self, molecule: Chem.Mol):
        for bond in molecule.GetBonds():
            bond_key = self.bond_key_generator(bond)
            try:
                self.bond_key_frequency_dictionary[bond_key] += 1
            except KeyError:
                self.bond_key_frequency_dictionary[bond_key] = 1


class EnvironmentKeyDictionary:
    def __init__(self,
        atom_hasher: AtomHasher,
        env_radius: int,
        env_key_frequency_threshold: int = 1):
        self.atom_hasher = atom_hasher
        self.env_radius = env_radius
        self.env_key_frequency_dictionary = {}
        self.env_key_frequency_threshold = env_key_frequency_threshold
        self.cached_molecule_id = 0
        self.atom_hashes = [] # Cache
    
    def __call__(self,
        prior_atom: Chem.Atom,
        posterior_atom: Chem.Atom) -> Union[KeyConstraint, None]:
        if not posterior_atom:
            return None
        return KeyConstraint(
            key_generator=self.EnvironmentKey,
            key_frequency_dictionary=self.env_key_frequency_dictionary,
            key_frequency_threshold=self.env_key_frequency_threshold)

    def UpdateAtomHashes(self, molecule: Chem.Mol):
        self.cached_molecule_id += 1
        molecule.SetIntProp("HID", self.cached_molecule_id)
        self.atom_hashes = []
        for atom in molecule.GetAtoms():
            self.atom_hashes.append(self.atom_hasher(atom))

    def EnvironmentKey(self, atom: Chem.Atom) -> int:
        molecule = atom.GetOwningMol()
        try:
            molecule_id = molecule.GetIntProp("HID")
            if molecule_id != self.cached_molecule_id:
                self.UpdateAtomHashes(molecule)
        except KeyError:
            self.UpdateAtomHashes(molecule)
        environment = mpt.CircularAtomicEnvironment(atom, self.env_radius)
        return environment.Key(self.atom_hashes)
        
    def AddMolecule(self, molecule: Chem.Mol):
        for atom in molecule.GetAtoms():
            env_key = self.EnvironmentKey(atom)
            try:
                self.env_key_frequency_dictionary[env_key] += 1
            except KeyError:
                self.env_key_frequency_dictionary[env_key] = 1