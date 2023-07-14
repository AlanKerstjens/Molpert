# About

Molpert is a library for graph-based perturbation of RDKit molecules. It supports atom and bond property modifications, insertions and deletions and combinations thereof. A plethora of settings provide fine grained control over how perturbations are executed. The library can be used stochastically, applying a random eligible perturbation, or deterministically, enumerating all eligible perturbations. Perturbation eligibility is determined by arbitrary user-provided constraints.

# Installation

## Installation from source

### Prerequisites

Ensure the following dependencies are installed:

* [RDKit](https://rdkit.org/)
* [Boost](https://www.boost.org/). You already have this if you installed the RDKit. If you'd like to build the Python bindings make sure Boost.Python is installed.
* [CMake](https://cmake.org/). Only necessary to build the Python bindings.

### Instructions

The following instructions are for GNU+Linux. For alternative operating systems you'll have to adapt these commands slightly.

```shell
git clone https://github.com/AlanKerstjens/Molpert.git
```

Molpert is header-only. If all you care about is the C++ API you are done. If you want the Python bindings too you must build them.

```shell
export MOLPERT="$(pwd)/Molpert"
mkdir ${MOLPERT}/build && cd ${MOLPERT}/build
cmake ..
make install
```

To be able to import the library from Python add `${MOLPERT}/lib` to your `${PYTHONPATH}`. Consider doing so in your `bash_profile` file. Otherwise you'll have to manually extend `${PYTHONPATH}` everytime you open a new shell.

```shell
export PYTHONPATH="${PYTHONPATH}:${MOLPERT}/lib"
```

### Troubleshooting

CMake will try to find the dependencies for you. To avoid problems ensure you build Molpert with the same Boost and Python versions that you used to build the RDKit. If CMake finds a different Boost or Python installation you'll need to point it to the correct one, as described [here](https://cmake.org/cmake/help/latest/module/FindBoost.html) and [here](https://cmake.org/cmake/help/latest/module/FindPython.html).

CMake will search for the RDKit in the active Anaconda environment (if you have one) and at `${RDBASE}` if set. If neither of these are the case you need to specify the path to the RDKit yourself. Replace the above CMake command with the one below, substituting the `<placeholder/path>` with your path.

```shell
cmake -DRDKit_ROOT=<path/to/rdkit> ..
```

# Getting started

This section provides a high level overview of how to use Molpert and its capabilities. It's not intended to serve as comprehensive documentation, and boilerplate code is omitted. Example code snippets are provided for C++ and Python (in that order).

To start using Molpert include or import the relevant headers or package.

```cpp
#include "MoleculePerturber.hpp"
```
```python
import molpert as mpt
```

## What are perturbations?

Perturbations are specifications of how a molecule ought to be modified. They are represented as classes that inherit from the `MolecularPerturbation` base class. There are 8 perturbation types of interest to the average user. Each derived class has a corresponding entry in the `MolecularPerturbation::Type` enum.

| Derived class             | `MolecularPerturbation::Type` | Classification |
| :-----------------------: | :---------------------------: | :------------: |
| `AtomicNumberChange`      | `AtomicNumberChange_t`        | Decorative     |
| `FormalChargeChange`      | `FormalChargeChange_t`        | Decorative     |
| `ExplicitHydrogensChange` | `ExplicitHydrogensChange_t`   | Decorative     |
| `BondTypeChange`          | `BondTypeChange_t`            | Decorative     |
| `AtomInsertion`           | `AtomInsertion_t`             | Topological    |
| `AtomDeletion`            | `AtomDeletion_t`              | Topological    |
| `BondInsertion`           | `BondInsertion_t`             | Topological    |
| `BondDeletion`            | `BondDeletion_t`              | Topological    |

After construction a perturbation can be executed by calling its `operator()`. Said operator has two overloads: one to perturb mutable `RDKit::RWMol`s in place, and one to perturb copies of immutable `RDKit::ROMol`s. In Python the RDKit makes no distinction between `RDKit::ROMol` and `RDKit::RWMol`, exposing a single mutable `Chem.Mol` class instead. Perturbations are always executed in place by calling the `mpt.ExecutePerturbation()` free function. If you want to create a copy you will have to do so manually.

```cpp
RDKit::RWMol rwmol (...);
RDKit::ROMol romol (...);

std::shared_ptr<MolecularPerturbation> perturbation (...);

perturbation(rwmol); // Perturb in place.
RDKit::RWMol perturbed_romol = perturbation(romol); // Copy and then perturb.
```
```python
molecule = Chem.Mol(...)
molecule_copy = Chem.Mol(molecule)

perturbation = mpt.MolecularPerturbation(...)

mpt.ExecutePerturbation(perturbation, molecule)
mpt.ExecutePerturbation(perturbation, molecule_copy)
```

You can characterize a perturbation by retrieving its type or calculating its ID (i.e. hash).

```cpp
MolecularPerturbation::Type perturbation_type = perturbation->GetType();
std::size_t perturbation_id = perturbation->ID();
```
```python
perturbation_type = perturbation.type
perturbation_id = perturbation.ID()
```

Compared to standard functions, the added value of `MolecularPerturbation`s is:
- Deferred execution
- Deduplication
- Constraint enforcement

## Tagging molecules
In the RDKit specific atoms and bonds can be referenced through their pointers (`RDKit::Atom*` and `RDKit::Bond*`) or indices. Since atoms and bonds are stored in contiguous memory both their addresses and indices may be invalidated when you modify the topology of a molecule. This makes keeping track of atoms and bonds in ever-changing molecules challenging.

To keep track of an atom or bond throughout a molecule's lifespan we use tags. Tags are integers attached to atoms and bonds. They remain associated with the atom or bond regardless of how we modify the molecule. Tags shouldn't be re-used within the same molecule.

**Before constructing perturbations you must tag your molecules! Perturbations that insert atoms or bonds will automatically tag the new atoms and bonds. You only need to tag the initial molecule.**

```cpp
RDKit::ROMol molecule (...);
TagMolecule(molecule);
```
```python
molecule = Chem.Mol(...)
mpt.TagMolecule(molecule)
```

To get the tag of an atom or bond:

```cpp
const RDKit::Atom* atom = molecule.getAtomWithIdx(0);
const RDKit::Bond* bond = molecule.getBondWithIdx(0);
std::size_t atom_tag = GetTag(atom);
std::size_t bond_tag = GetTag(bond);
```
```python
atom = molecule.GetAtomWithIdx(0)
bond = molecule.GetBondWithIdx(0)
atom_tag = mpt.GetAtomTag(atom)
bond_tag = mpt.GetBondTag(bond)
```

## Manually constructing perturbations

You can construct perturbations by manually calling their constructor. Multiple constructor overloads are exposed. Since Python isn't strongly typed overload selection on the Python side is a bit finnicky. You may have to specify all arguments for overload resolution to work properly.

```cpp
AtomicNumberChange perturbation1 (
  0, // Atom index of the affected atom
  8  // New atomic number
)

AtomDeletion perturbation2 (
  molecule, // Subject molecule
  0,        // Atom index of the atom to delete
  1         // Atom index of the reconnection atom that takes its place
)
```
```python
perturbation1 = mpt.AtomicNumberChange(
  atom_idx=0,
  atomic_number=8)

perturbation2 = mpt.AtomDeletion(
  molecule=molecule,
  atom_idx=0,
  reconnection_atom_idx=1,
  preserve_bond_types=True, # Default
  default_bond_type=Chem.BondType.SINGLE # Default
)
```

## Automatic perturbation generation

Molpert's real power rests in its `MolecularPerturbation` factory: `MoleculePerturber`.

```cpp
MoleculePerturber perturber;
```
```python
perturber = mpt.MoleculePerturber()
```

`MoleculePerturber` has a plethora settings to regulate which perturbations are generated. Reasonable default values are set for you. If you don't like them you can change them by accessing the corresponding public member variables. Broadly speaking settings can be classified into two categories:
- Perturbation specific settings
- Sampling values and weights

Each perturbation type has *perturbation specific settings*. Setting names are prefaced with the perturbation type's name. Most of these settings are unsigned integers or booleans.

```cpp
perturber.atom_insertion_max_n_neighbors = 4;   // Default is 3
perturber.bond_deletion_allow_reroutes = false; // Default is true
```
```python
perturber.atom_insertion_max_n_neighbors = 3
perturber.bond_deletion_allow_reroutes = False
```

*Sampling values* are vectors (`std::vector`) from which atom and bond properties are sampled. Each value array has a corresponding size-matched vector of *sampling weights*. The default values and weights are derived from ChEMBL. If you manually modify the value vector make sure you modify the weights vector accordingly. If you are a Python user you'll have to use the corresponding setter instead. In the following example you would set a 60% probability of sampling carbons and a 20% probability of sampling nitrogens and oxygens. Atomic numbers in particular are sampled during `AtomicNumberChange`s and `AtomInsertion`s.

```cpp
perturber.atomic_numbers = {6, 7, 8};
perturber.atomic_numbers_weights = {0.6, 0.2, 0.2};
```
```python
perturber.SetAtomicNumbers(
  atomic_numbers=[6, 7, 8],
  weights=[0.6, 0.2, 0.2])
```

You can conveniently modify some of these arrays by passing arguments to the `MoleculePerturber` constructor. This is less error prone than doing it manually.

```cpp
MoleculePerturber perturber (
  false, // Sample uniformly instead of using the ChEMBL sampling weights.
  true,  // Include aromatic bonds in the eligible bond types.
  true   // Allow sampling of aromatic bond type for acyclic bonds.
)
```
```python
perturber = MoleculePerturber(
  use_chembl_distribution=False,
  use_aromatic_bonds=True,
  acyclic_can_be_aromatic=True)
```

You can also adjust which perturbation types are sampled and with which probabilities. The main difference with the above is that perturbation sampling settings are fixed-size variables (`std::bitset` and `std::array`). It's recommended to modify them using `operator[]` instead of `operator=`. In Python you must use a setter function instead. By default all perturbations may be sampled, with sampling weights being derived from ChEMBL.

```cpp
perturber.perturbation_types.reset() // Reset all perturbation types to false.
perturber.perturbation_types[MolecularPerturbation::Type::AtomicNumberChange] = true;
perturber.perturbation_types[MolecularPerturbation::Type::BondTypeChange] = true;
perturber.perturbation_types_weights[MolecularPerturbation::Type::AtomicNumberChange] = 0.5;
perturber.perturbation_types_weights[MolecularPerturbation::Type::BondTypeChange] = 0.5;
```
```python
perturber.SetPerturbationTypes(
    perturbation_types=[
        mpt.MolecularPerturbation.Type.AtomicNumberChange,
        mpt.MolecularPerturbation.Type.BondTypeChange],
    weights=[0.5, 0.5],
    reset=True # Reset everything to 0 before setting the specified values (default)
)
```

### Stochastic perturbation construction

`MoleculePerturber` has member functions to randomly generate applicable perturbations. These functions rely on a random number generator which must be passed as argument to the perturbation generators.

```cpp
std::random_device rd;
std::mt19937 prng (rd()); // Pseudo-random number generator
```
```python
prng = mpt.PRNG()
```

If you don't mind the type of the generated perturbations you can call `MoleculePerturber::operator()`. In Python the equivalent member function is `MoleculePerturber.Perturbation()`.

```cpp
// Perturb the molecule in some random way
std::shared_ptr<MolecularPerturbation> perturbation = perturber(
  molecule, prng);
```
```python
# Perturb the molecule in some random way
perturbation = perturber.Perturbation(
  molecule=molecule,
  prng=prng)
```

To limit the generation to a specific type of perturbation call the appropiate member function. Stochastic perturbation generators are named using the perturbation type in imperative tense:

```cpp
// Insert an atom anywhere in molecule
std::shared_ptr<AtomInsertion> atom_insertion = perturber.InsertAtom(
  molecule, prng);
// Delete any bond of molecule
std::shared_ptr<BondDeletion> bond_deletion = perturber.DeleteBond(
  molecule, prng);
```
```python
# Insert an atom anywhere in molecule
atom_insertion = perturber.InsertAtom(
  molecule=molecule,
  prng=prng)
# Delete any bond of molecule
bond_deletion = perturber.DeleteBond(
  molecule=molecule,
  prng=prng)
```

Generators have multiple overloads allowing you to specify, among other things, which part of the molecule (i.e. which atom or bond) to target.

```cpp
// Insert an atom adjacent to atom with index 5
std::shared_ptr<AtomInsertion> atom_insertion = perturber.InsertAtom(
  molecule,
  5, // Central atom index
  prng);
// Delete the bond with index 3
std::shared_ptr<BondDeletion> bond_deletion = perturber.DeleteBond(
  molecule,
  3, // Deleted bond index
  prng);
```
```python
# Insert an atom adjacent to atom with index 5
atom_insertion = perturber.InsertAtom(
    molecule=molecule,
    central_atom_idx=5,
    prng=prng)
# Delete the bond with index 3
bond_deletion = perturber.DeleteBond(
    molecule=molecule,
    bond_idx=5,
    prng=prng)
```

### Deterministic perturbation construction

`MoleculePerturber` has member functions to enumerate all perturbations that are applicable to a molecule. Enumerated perturbations are stored in a queue, which must be passed as argument to the perturbation generators.

```cpp
MolecularPerturbationQueue queue;
```
```python
queue = mpt.MolecularPerturbationQueue()
```

`MolecularPerturbationQueue` automatically deduplicates enumerated perturbations, ensuring that the same perturbation can't be stored more than once. Beware that a perturbation is considered a duplicate if it ever entered the queue, regardless of whether it's in it at present time. If you don't like this behaviour call `MolecularPerturbationQueue::clear()` to erase the record of previously queued perturbations.

If you'd like to enumerate all perturbations of all kinds that can be applied to a molecule call `MoleculePerturber::operator()`. The equivalent function in Python is `MoleculePerturber.Perturbations()`.

```cpp
// Enumerate all perturbations that can be applied to molecule
perturber(queue, molecule);
```
```python
# Enumerate all perturbations that can be applied to molecule
perturber.Perturbations(
  queue=queue,
  molecule=molecule)
```

To limit the generation to a specific type of perturbation call the corresponding member function. Deterministic perturbation generators are named after the perturbation type in plural:

```cpp
// Enumerate all atom insertions that can be applied to molecule
perturber.AtomInsertions(queue, molecule);
// Enumerate all bond deletions that can be applied to molecule
perturber.BondDeletions(queue, molecule);
```
```python
# Enumerate all atom insertions that can be applied to molecule
perturber.AtomInsertions(queue=queue, molecule=molecule)
# Enumerate all bond deletions that can be applied to molecule
perturber.BondDeletions(queue=queue, molecule=molecule)
```

Generators have multiple overloads allowing you to specify, among other things, which part of the molecule (i.e. which atom or bond) to target.

```cpp
// Enumerate all atom insertions that insert an atom adjacent to atom 5
perturber.AtomInsertions(
  queue,
  molecule,
  5); // Central atom index
// Enumerate all bond deletions that deleted bond 3
perturber.BondDeletions(
  queue,
  molecule,
  3); // Deleted bond index
```
```python
# Enumerate all atom insertions that insert an atom adjacent to atom 5
perturber.AtomInsertions(
  queue=queue,
  molecule=molecule,
  central_atom_idx=5)
# Enumerate all bond deletions that deleted bond 3
perturber.BondDeletions(
  queue=queue,
  molecule=molecule,
  bond_idx=3)
```

### Perturbation constraints

`MoleculePerturber`s perturbation generation functions accept `MolecularConstraints` as optional arguments. Constraints allow you to specify arbitrary requirements that molecules resulting from a perturbation must fulfill.

A constraint is a callback function that takes as input an atom, bond or molecule and returns a boolean, evaluating to `true` if the corresponding atom, bond or molecule is allowed, and `false` otherwise. Atom and bond constraints apply to a specific atom or bond. Molecule constraints apply to the whole molecule. 

The C++ signatures of constraints are:

```cpp
typedef std::function<bool(const RDKit::Atom*)> AtomConstraint;
typedef std::function<bool(const RDKit::Bond*)> BondConstraint;
typedef std::function<bool(const RDKit::ROMol*)> MoleculeConstraint;
```

Let us define some constraints. Constraints can be free functions or functors:

```cpp
bool HasPositiveCharge(const RDKit::Atom* atom) {
  return atom->getFormalCharge() > 0;
};

bool IsHeteroatom(const RDKit::Atom* atom) {
  int z = atom->getAtomicNum();
  return z != 1 && z != 6;
};

bool NotHeteroBond(const RDKit::Bond* bond) {
  return !(IsHeteroatom(bond->getBeginAtom()) && 
           IsHeteroatom(bond->getEndAtom()));
};

struct IsRigid {
  unsigned max_n_rotatable_bonds = 0;

  IsRigid(unsigned max_n_rotatable_bonds) :
    max_n_rotatable_bonds(max_n_rotatable_bonds) {};

  bool operator()(const RDKit::ROMol* molecule) const {
    unsigned n_rotatable_bonds = RDKit::Descriptors::calcNumRotatableBonds(
      *molecule);
    return n_rotatable_bonds <= max_n_rotatable_bonds;
  };
};
```
```python
def HasPositiveCharge(atom: Chem.Atom) -> bool:
    return atom.GetFormalCharge() > 0

def IsHeteroatom(atom: Chem.Atom) -> bool:
    z = atom.GetAtomicNum()
    return z != 1 and z != 6

def NotHeteroBond(bond: Chem.Bond) -> bool:
    return not (IsHeteroatom(bond.GetBeginAtom()) and 
                IsHeteroatom(bond.GetEndAtom()))

class IsRigid:
    def __init__(self, max_n_rotatable_bonds: int):
        self.max_n_rotatable_bonds = max_n_rotatable_bonds
    
    def __call__(self, molecule: Chem.Mol) -> bool:
        n_rotatable_bonds = Chem.Descriptors.NumRotatableBonds(molecule)
        return n_rotatable_bonds <= self.max_n_rotatable_bonds
```

You can store these constraints in a `MolecularConstraints` instance.

```cpp
MolecularConstraints constraints;
// Atom with tag 3 should have a positive charge.
constraints.SetAtomConstraint(3, HasPositiveCharge);
// Bond with tag 2 shouldn't bond two heteroatoms.
constraints.SetBondConstraint(2, NotHeteroBond);
// Molecule may not have more than 4 rotatable bonds.
constraints.SetMoleculeConstraint(IsRigid(4));
```
```python
constraints = mpt.MolecularConstraints()
# Atom with tag 3 should have a positive charge.
constraints.SetAtomConstraint(
    atom_tag=3,
    atom_constraint=HasPositiveCharge)
# Bond with tag 2 shouldn't bond two heteroatoms.
constraints.SetBondConstraint(
    bond_tag=2,
    bond_constraint=NotHeteroBond)
# Molecule may not have more than 4 rotatable bonds.
constraints.SetMoleculeConstraint(
    molecule_constraint=IsRigid(4))
```

Pass the `MolecularConstraints` to the perturbation generator as the (optional) last argument:

```cpp
// Stochastic perturbation generation
std::shared_ptr<AtomInsertion> perturbation = perturber.InsertAtom(
  molecule, prng, &constraints);
// Deterministic perturbation generation
perturber.AtomInsertions(queue, molecule, &constraints);
```
```python
# Stochastic perturbation generation
perturbation = perturber.InsertAtom(
    molecule=molecule,
    prng=prng,
    constraints=constraints)
# Deterministic perturbation generation
perturber.AtomInsertions(
    queue=queue,
    molecule=molecule,
    constraints=constraints)
```

You can store as many `MoleculeConstraint`s as you please, but only up to one `AtomConstraint` or `BondConstraint` per atom or bond. If you need to check for multiple conditions merge them into a single function.

```cpp
constraints.SetAtomConstraint(3, HasPositiveCharge); // OK
constraints.SetAtomConstraint(3, IsHeteroatom); // Careful! Overwrites constraint.
```

If you need all of your atoms and/or bonds to fulfill some condition you have two options:
- Option A (recommended): Use a `MoleculeConstraint` that iterates over all atoms and/or bonds and checks each one for the condition.
- Option B: Define constraint generators and construct `MolecularConstraints` with them. A constraint generator is a function that takes as input a pair of atoms or bonds (prior and posterior), compares both and generates or updates the constraint for the atom or bond whenever a molecule is perturbed. Use constraint generators only if your constraints depend on both the prior and posterior states. The signature of constraint generators is as follows:

```cpp
typedef std::function<AtomConstraint(const RDKit::Atom*, const RDKit::Atom*)
> AtomConstraintGenerator;
typedef std::function<BondConstraint(const RDKit::Bond*, const RDKit::Bond*)
> BondConstraintGenerator;
```

### Performance advice

During topological perturbation generation the molecular graph may have to be traversed to determine if a perturbation is applicable or not. If you plan on perturbing a molecule just once it's advisable to perform these graph traversals lazily (default).

```cpp
perturber.assess_connectivity_with_sssr = false;
perturber.assess_distances_with_distance_matrix = false;
```

If you will be enumerating many neighbors of a molecule it's worthwhile to spend some extra time calculating some topological information once to be used in subsequent perturbation generation:

```cpp
perturber.assess_connectivity_with_sssr = true;
perturber.assess_distances_with_distance_matrix = true;
```

If you'd like to remove an upper limit for a setting set it to `std::numeric_limits<unsigned>::max()`. On most systems the Python equivalent of this is `0xFFFFFFFF`. The code will check against this value and skip unnecessary calculations. DO NOT set it to an arbitrarily high value, as the graph traversals will still be executed.

```cpp
perturber.bond_insertion_max_distance_partner = std::numeric_limits<unsigned>::max(); // OK
perturber.bond_insertion_max_distance_partner = perturber.max_unsigned; // OK. Same as above.
perturber.bond_insertion_max_distance_partner = 9999; // BAD. Will trigger graph traversal.
```
```python
perturber.bond_insertion_max_distance_partner = 0xFFFFFFFF # OK
```


Be conservative in your use of constraints. It may be possible to achieve what you want without them by changing the `MoleculePerturber` settings. Even if changing settings doesn't get you all the way there it can greatly reduce the number of generated perturbations and with it the number of times your constraints are evaluated. For example, if you'd like to avoid generating macrocycles consider using the below settings. They won't guarantee it doesn't happen, but they will reduce the number of macrocycles that are generated in a computationally more efficient way.

```cpp
// 5 consecutive bonds span 6 atoms.
// Bonding the peripheral atoms results in a 6-membered ring.
perturber.bond_insertion_max_distance_partner = 5;
// 4 consecutive bonds span 5 atoms. 
// Inserting an atom between the peripheral atoms results in a 6-membered ring.
perturber.atom_insertion_max_distance_neighbor = 4;
// 6 consecutive bonds span 7 atoms. However, we are removing one of them.
// Bonding the peripheral atoms results in a 6-membered ring.
perturber.atom_deletion_max_distance_reconnection_atom = 6;
```