#pragma once
#ifndef _MOLECULAR_PERTURBATIONS_HPP_
#define _MOLECULAR_PERTURBATIONS_HPP_

#include "MolecularGraphProjection.hpp"
#include <GraphMol/RWMol.h>

// typedefs so we can use less bits if we really need to
typedef std::size_t AtomIdx;
typedef std::size_t BondIdx;

// Abstract base class for molecular perturbations
class MolecularPerturbation {
public:
  enum Type : std::size_t {
    AtomicNumberChange_t = 0,
    FormalChargeChange_t = 1,
    ExplicitHydrogensChange_t = 2,
    BondTypeChange_t = 3,
    AtomConstruction_t = 4,
    AtomDestruction_t = 5,
    BondConstruction_t = 6,
    BondDestruction_t = 7,
    TopologicalPerturbation_t = 8,
    AtomInsertion_t = 9,
    AtomDeletion_t = 10,
    BondInsertion_t = 11,
    BondDeletion_t = 12,
    SubgraphConstruction_t = 13,
    SubgraphDestruction_t = 14,
    SubgraphPerturbation_t = 15
  };
  static const std::size_t n_types = 16;
  typedef std::bitset<n_types> TypeMask;

  enum TargetType {Atom, Bond};
  typedef std::pair<TargetType, std::size_t> Target;

public:
  // In-place perturbation
  virtual void operator()(RDKit::RWMol&) const = 0;
  // Copy perturbation
  RDKit::RWMol operator()(const RDKit::ROMol& molecule) const {
    RDKit::RWMol perturbed_molecule (molecule);
    operator()(perturbed_molecule);
    return perturbed_molecule;
  };
  // Simulate the changes on the molecular graph associated with the perturbation
  virtual void ProjectMolecularGraph(MolecularGraphProjection&) const = 0;
  // Identifier (e.g. hash) for deduplication purposes
  virtual std::size_t ID() const = 0;
  // Perturbation type getter
  virtual Type GetType() const = 0;
  friend auto operator<=>(
    const MolecularPerturbation&, const MolecularPerturbation&) = default;
};


class AtomicNumberChange : public MolecularPerturbation {
  static const Type type = Type::AtomicNumberChange_t;
  AtomIdx atom_idx;
  std::uint8_t atomic_number;

public:
  AtomicNumberChange(AtomIdx atom_idx, std::uint8_t atomic_number) :
    atom_idx(atom_idx), atomic_number(atomic_number) {};

  void operator()(RDKit::RWMol& molecule) const {
    RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    atom->setAtomicNum(atomic_number);
  };

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    MolecularGraphProjection::Atom& atom = projection.EditAtom(atom_idx);
    atom.atomic_number = atomic_number;
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, atom_idx);
    boost::hash_combine(id, atomic_number);
    return id;
  };

  Type GetType() const { return type; };
};


class FormalChargeChange : public MolecularPerturbation {
  static const Type type = Type::FormalChargeChange_t;
  AtomIdx atom_idx;
  std::int8_t formal_charge;

public:
  FormalChargeChange(AtomIdx atom_idx, std::int8_t formal_charge) :
    atom_idx(atom_idx), formal_charge(formal_charge) {};

  void operator()(RDKit::RWMol& molecule) const {
    RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    atom->setFormalCharge(formal_charge);
  };

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    MolecularGraphProjection::Atom& atom = projection.EditAtom(atom_idx);
    atom.formal_charge = formal_charge;
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, atom_idx);
    boost::hash_combine(id, formal_charge);
    return id;
  };

  Type GetType() const { return type; };
};


class ExplicitHydrogensChange : public MolecularPerturbation {
  static const Type type = Type::ExplicitHydrogensChange_t;
  AtomIdx atom_idx;
  std::uint8_t n_explicit_hydrogens;

public:
  ExplicitHydrogensChange(AtomIdx atom_idx, std::uint8_t n_explicit_hydrogens) :
    atom_idx(atom_idx), n_explicit_hydrogens(n_explicit_hydrogens) {};

  void operator()(RDKit::RWMol& molecule) const {
    RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    atom->setNumExplicitHs(n_explicit_hydrogens);
  };

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    MolecularGraphProjection::Atom& atom = projection.EditAtom(atom_idx);
    atom.n_explicit_hydrogens = n_explicit_hydrogens;
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, atom_idx);
    boost::hash_combine(id, n_explicit_hydrogens);
    return id;
  };

  Type GetType() const { return type; };
};


class BondTypeChange : public MolecularPerturbation {
  static const Type type = Type::BondTypeChange_t;
  AtomIdx begin_atom_idx, end_atom_idx;
  RDKit::Bond::BondType bond_type;

public:
  BondTypeChange(
    AtomIdx begin_atom_idx,
    AtomIdx end_atom_idx,
    RDKit::Bond::BondType bond_type) :
    begin_atom_idx(begin_atom_idx),
    end_atom_idx(end_atom_idx),
    bond_type(bond_type) {
    if (begin_atom_idx == end_atom_idx) {
      throw std::invalid_argument("Self-bond can't exist");
    };
    if (end_atom_idx < begin_atom_idx) {
      std::swap(begin_atom_idx, end_atom_idx);
    };
  };

  BondTypeChange(
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    RDKit::Bond::BondType bond_type) :
    bond_type(bond_type) {
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    begin_atom_idx = bond->getBeginAtomIdx();
    end_atom_idx = bond->getEndAtomIdx();
    if (end_atom_idx < begin_atom_idx) {
      std::swap(begin_atom_idx, end_atom_idx);
    };
  };

  void operator()(RDKit::RWMol& molecule) const {
    RDKit::Bond* bond =
      molecule.getBondBetweenAtoms(begin_atom_idx, end_atom_idx);
    bond->setBondType(bond_type);
  };

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    MolecularGraphProjection::Bond& bond = projection.EditBond(
      begin_atom_idx, end_atom_idx);
    bond.bond_type = bond_type;
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, begin_atom_idx);
    boost::hash_combine(id, end_atom_idx);
    boost::hash_combine(id, bond_type);
    return id;
  };

  Type GetType() const { return type; };
};


class AtomConstruction : public MolecularPerturbation {
  static const Type type = Type::AtomConstruction_t;
  Tag atom_tag;
  std::uint8_t atomic_number;
  std::int8_t formal_charge;
  std::uint8_t n_explicit_hydrogens;

public:
  AtomConstruction(
    Tag atom_tag,
    std::uint8_t atomic_number = 6,
    std::int8_t formal_charge = 0,
    std::uint8_t n_explicit_hydrogens = 0) :
    atom_tag(atom_tag),
    atomic_number(atomic_number),
    formal_charge(formal_charge),
    n_explicit_hydrogens(n_explicit_hydrogens) {};

  AtomConstruction(
    const RDKit::ROMol& molecule,
    std::uint8_t atomic_number = 6,
    std::int8_t formal_charge = 0,
    std::uint8_t n_explicit_hydrogens = 0) :
    atom_tag(GetMaxAtomTag(molecule) + 1),
    atomic_number(atomic_number),
    formal_charge(formal_charge),
    n_explicit_hydrogens(n_explicit_hydrogens) {};

  void operator()(RDKit::RWMol& molecule) const {
    AtomIdx atom_idx = molecule.addAtom();
    RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    TagAtom(molecule, atom, atom_tag);
    atom->setAtomicNum(atomic_number);
    atom->setFormalCharge(formal_charge);
    atom->setNumExplicitHs(n_explicit_hydrogens);
  };

  friend auto operator<=>(
    const AtomConstruction&, const AtomConstruction&) = default;

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    projection.AddAtom(
      atom_tag, atomic_number, formal_charge, n_explicit_hydrogens);
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, atom_tag);
    boost::hash_combine(id, atomic_number);
    boost::hash_combine(id, formal_charge);
    boost::hash_combine(id, n_explicit_hydrogens);
    return id;
  };

  Type GetType() const { return type; };
};


class AtomDestruction : public MolecularPerturbation {
  static const Type type = Type::AtomDestruction_t;
  AtomIdx atom_idx;

public:
  AtomDestruction(AtomIdx atom_idx) : atom_idx(atom_idx) {};

  void operator()(RDKit::RWMol& molecule) const {
    molecule.removeAtom(atom_idx);
  };

  friend auto operator<=>(
    const AtomDestruction&, const AtomDestruction&) = default;

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    projection.RemoveAtom(atom_idx);
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, atom_idx);
    return id;
  };

  Type GetType() const { return type; };
};


class BondConstruction : public MolecularPerturbation {
  static const Type type = Type::BondConstruction_t;
  Tag bond_tag;
  AtomIdx begin_atom_idx;
  AtomIdx end_atom_idx;
  RDKit::Bond::BondType bond_type;

public:
  BondConstruction(
    Tag bond_tag,
    AtomIdx begin_atom_idx,
    AtomIdx end_atom_idx,
    RDKit::Bond::BondType bond_type = RDKit::Bond::SINGLE) :
    bond_tag(bond_tag),
    begin_atom_idx(begin_atom_idx),
    end_atom_idx(end_atom_idx),
    bond_type(bond_type) {
    if (begin_atom_idx == end_atom_idx) {
      throw std::invalid_argument("Can't create self-bond");
    };
    if (end_atom_idx < begin_atom_idx) {
      std::swap(begin_atom_idx, end_atom_idx);
    };
  };

  BondConstruction(
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx,
    AtomIdx end_atom_idx,
    RDKit::Bond::BondType bond_type = RDKit::Bond::SINGLE) :
    BondConstruction(
      GetMaxBondTag(molecule) + 1,
      begin_atom_idx,
      end_atom_idx,
      bond_type) {};

  void operator()(RDKit::RWMol& molecule) const {
    if (molecule.getBondBetweenAtoms(begin_atom_idx, end_atom_idx)) {
      return;
    };
    molecule.addBond(begin_atom_idx, end_atom_idx, bond_type);
    RDKit::Bond* bond = molecule.getBondBetweenAtoms(
      begin_atom_idx, end_atom_idx);
    TagBond(molecule, bond, bond_tag);
  };

  // Note that the bond tag is included in the comparison. This means that two 
  // instances differing solely in their bond tag aren't considered equal, and 
  // may be present simultaneously in containers enforcing uniqueness.
  friend auto operator<=>(
    const BondConstruction&, const BondConstruction&) = default;

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    projection.AddBond(bond_tag, begin_atom_idx, end_atom_idx, bond_type);
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, bond_tag);
    boost::hash_combine(id, begin_atom_idx);
    boost::hash_combine(id, end_atom_idx);
    boost::hash_combine(id, bond_type);
    return id;
  };

  Type GetType() const { return type; };
};


class BondDestruction : public MolecularPerturbation {
  static const Type type = Type::BondDestruction_t;
  AtomIdx begin_atom_idx;
  AtomIdx end_atom_idx;

public:
  BondDestruction(
    AtomIdx begin_atom_idx,
    AtomIdx end_atom_idx) :
    begin_atom_idx(begin_atom_idx),
    end_atom_idx(end_atom_idx) {
    if (begin_atom_idx == end_atom_idx) {
      throw std::invalid_argument("Self-bond can't exist");
    };
    if (end_atom_idx < begin_atom_idx) {
      std::swap(begin_atom_idx, end_atom_idx);
    };
  };

  void operator()(RDKit::RWMol& molecule) const {
    molecule.removeBond(begin_atom_idx, end_atom_idx);
  };

  friend auto operator<=>(
    const BondDestruction&, const BondDestruction&) = default;

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    projection.RemoveBond(begin_atom_idx, end_atom_idx);
  };

  std::size_t ID() const {
    std::size_t id = type;
    boost::hash_combine(id, begin_atom_idx);
    boost::hash_combine(id, end_atom_idx);
    return id;
  };

  AtomIdx GetBeginAtomIdx() const {
    return begin_atom_idx;
  };

  AtomIdx GetEndAtomIdx() const {
    return end_atom_idx;
  };

  Type GetType() const { return type; };
};


class TopologicalPerturbation : public MolecularPerturbation {
  static const Type type = Type::TopologicalPerturbation_t;

protected:
  // Sets are used for de-duplication purposes and to standardize the ID.
  std::set<AtomConstruction> atom_constructions;
  std::set<BondConstruction> bond_constructions;
  std::set<BondDestruction> bond_destructions;
  // Sort atoms in descending order to avoid index invalidation during deletion.
  std::set<AtomDestruction, std::greater<AtomDestruction>> atom_destructions;

public:
  template <class... Args>
  bool AddAtomConstruction(Args&&... args) {
    auto [it, emplaced] = atom_constructions.emplace(
      std::forward<Args>(args)...);
    return emplaced;
  };

  template <class... Args>
  bool AddBondConstruction(Args&&... args) {
    auto [it, emplaced] = bond_constructions.emplace(
      std::forward<Args>(args)...);
    return emplaced;
  };

  template <class... Args>
  bool AddBondDestruction(Args&&... args) {
    auto [it, emplaced] = bond_destructions.emplace(
      std::forward<Args>(args)...);
    return emplaced;
  };

  template <class... Args>
  bool AddAtomDestruction(Args&&... args) {
    auto [it, emplaced] = atom_destructions.emplace(
      std::forward<Args>(args)...);
    return emplaced;
  };

  void operator()(RDKit::RWMol& molecule) const {
    for (const AtomConstruction& atom_construction : atom_constructions) {
      atom_construction(molecule);
    };
    for (const BondConstruction& bond_construction : bond_constructions) {
      bond_construction(molecule);
    };
    for (const BondDestruction& bond_destruction : bond_destructions) {
      bond_destruction(molecule);
    };
    for (const AtomDestruction& atom_destruction : atom_destructions) {
      atom_destruction(molecule);
    };
  };

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    for (const AtomConstruction& atom_construction : atom_constructions) {
      atom_construction.ProjectMolecularGraph(projection);
    };
    for (const BondConstruction& bond_construction : bond_constructions) {
      bond_construction.ProjectMolecularGraph(projection);
    };
    for (const BondDestruction& bond_destruction : bond_destructions) {
      bond_destruction.ProjectMolecularGraph(projection);
    };
    for (const AtomDestruction& atom_destruction : atom_destructions) {
      atom_destruction.ProjectMolecularGraph(projection);
    };
  };

  std::size_t ID() const {
    std::size_t id = type;
    for (const AtomConstruction& atom_construction : atom_constructions) {
      boost::hash_combine(id, atom_construction.ID());
    };
    for (const BondConstruction& bond_construction : bond_constructions) {
      boost::hash_combine(id, bond_construction.ID());
    };
    for (const BondDestruction& bond_destruction : bond_destructions) {
      boost::hash_combine(id, bond_destruction.ID());
    };
    for (const AtomDestruction& atom_destruction : atom_destructions) {
      boost::hash_combine(id, atom_destruction.ID());
    };
    return id;
  };

  const std::set<AtomConstruction>& GetAtomConstructions() const {
    return atom_constructions;
  };

  const std::set<BondConstruction>& GetBondConstructions() const {
    return bond_constructions;
  };

  const std::set<BondDestruction>& GetBondDestructions() const {
    return bond_destructions;
  };

  const std::set<AtomDestruction, std::greater<AtomDestruction>>&
  GetAtomDestructions() const {
    return atom_destructions;
  };

  Type GetType() const { return type; };
};


class AtomInsertion : public TopologicalPerturbation {
  static const Type type = Type::AtomInsertion_t;

public:
  AtomInsertion(
    const RDKit::ROMol& molecule,
    std::uint8_t atomic_number = 6,
    std::int8_t formal_charge = 0,
    std::uint8_t n_explicit_hydrogens = 0) {
    atom_constructions.emplace(
      molecule, atomic_number, formal_charge, n_explicit_hydrogens);
  };

  AtomInsertion(
    const RDKit::ROMol& molecule,
    const std::vector<AtomIdx>& neighbor_atom_indices,
    const std::vector<RDKit::Bond::BondType>& bond_types,
    std::uint8_t atomic_number = 6,
    std::int8_t formal_charge = 0,
    std::uint8_t n_explicit_hydrogens = 0,
    int dropped_atom_idx = -1) {
    std::size_t n = neighbor_atom_indices.size();
    if (n != bond_types.size()) {
      throw std::length_error("Size mismatch between input vectors");
    };
    AtomIdx atom_idx = molecule.getNumAtoms();
    atom_constructions.emplace(
      molecule, atomic_number, formal_charge, n_explicit_hydrogens);
    Tag max_bond_tag = GetMaxBondTag(molecule);
    for (std::size_t i = 0; i < n; ++i) {
      AtomIdx neighbor_idx = neighbor_atom_indices[i];
      bond_constructions.emplace(
        ++max_bond_tag, atom_idx, neighbor_idx, bond_types[i]);
      // Inserting an atom may involve creating multiple bonds (between the
      // inserted atom and the specified neighboring atoms), which can create
      // very complex cyclic molecular topologies if no other bonds are
      // destroyed. The user may specify an atom whose bonds with the specified
      // neighboring atoms are destroyed. This yields more relaxed topologies,
      // which may look more natural to the user.
      if (dropped_atom_idx < 0) {
        continue;
      };
      const RDKit::Bond* dropped_bond =
        molecule.getBondBetweenAtoms(dropped_atom_idx, neighbor_idx);
      if (dropped_bond) {
        bond_destructions.emplace(dropped_atom_idx, neighbor_idx);
      };
    };
  };

  Type GetType() const { return type; };
};


class AtomDeletion : public TopologicalPerturbation {
  static const Type type = Type::AtomDeletion_t;

public:
  AtomDeletion(AtomIdx atom_idx) {
    atom_destructions.emplace(atom_idx);
  };

  AtomDeletion(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    AtomIdx reconnection_atom_idx,
    bool preserve_bond_types = true,
    RDKit::Bond::BondType default_bond_type = RDKit::Bond::SINGLE) {
    if (atom_idx == reconnection_atom_idx) {
      throw std::invalid_argument("Can't use deleted atom for reconnection");
    };
    // Atom deletions can disconnect the molecular graph. If we want to prevent
    // this while still enabling atom deletions we can follow up the deletion
    // with some new bond formations (i.e. "reconnections").
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    Tag max_bond_tag = GetMaxBondTag(molecule);
    for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
      AtomIdx neighbor_idx = bond->getOtherAtomIdx(atom_idx);
      if (reconnection_atom_idx == neighbor_idx) {
        continue;
      };
      if (molecule.getBondBetweenAtoms(reconnection_atom_idx, neighbor_idx)) {
        continue;
      };
      RDKit::Bond::BondType bond_type =
        preserve_bond_types ? bond->getBondType() : default_bond_type;
      bond_constructions.emplace(
        ++max_bond_tag, reconnection_atom_idx, neighbor_idx, bond_type);
    };
    atom_destructions.emplace(atom_idx);
  };

  AtomDeletion(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    AtomIdx reconnection_center_atom_idx,
    const std::vector<AtomIdx>& reconnection_partner_atom_indices,
    const std::vector<RDKit::Bond::BondType> reconnection_bond_types) {
    if (atom_idx == reconnection_center_atom_idx) {
      throw std::invalid_argument("Can't use deleted atom for reconnection");
    };
    std::size_t n = reconnection_partner_atom_indices.size();
    if (n != reconnection_bond_types.size()) {
      throw std::length_error("Size mismatch between input vectors");
    };
    Tag max_bond_tag = GetMaxBondTag(molecule);
    for (std::size_t i = 0; i < n; ++i) {
      AtomIdx partner_atom_idx = reconnection_partner_atom_indices[i];
      if (reconnection_center_atom_idx == partner_atom_idx) {
        continue;
      };
      if (molecule.getBondBetweenAtoms(
        reconnection_center_atom_idx, partner_atom_idx)) {
        continue;
      };
      bond_constructions.emplace(++max_bond_tag, reconnection_center_atom_idx,
        partner_atom_idx, reconnection_bond_types[i]);
    };
    atom_destructions.emplace(atom_idx);
  };

  Type GetType() const { return type; };
};


// Bond insertions don't have more depth to them than a simple bond construction.
// Just for convenience we create an alias.
typedef BondConstruction BondInsertion;


class BondDeletion : public TopologicalPerturbation {
  static const Type type = Type::BondDeletion_t;

public:
  BondDeletion(
    AtomIdx begin_atom_idx,
    AtomIdx end_atom_idx) {
    bond_destructions.emplace(begin_atom_idx, end_atom_idx);
  };

  BondDeletion(
    const RDKit::ROMol& molecule,
    BondIdx bond_idx) {
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    bond_destructions.emplace(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
  };

  BondDeletion(
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx,
    AtomIdx end_atom_idx,
    AtomIdx reroute_begin_atom_idx,
    AtomIdx reroute_end_atom_idx,
    bool preserve_reroute_bond_type = true,
    RDKit::Bond::BondType default_bond_type = RDKit::Bond::SINGLE) {
    const RDKit::Bond* bond =
      molecule.getBondBetweenAtoms(begin_atom_idx, end_atom_idx);
    RDKit::Bond::BondType bond_type = preserve_reroute_bond_type ?
      bond->getBondType() : default_bond_type;
    if (!molecule.getBondBetweenAtoms(
      reroute_begin_atom_idx, reroute_end_atom_idx)) {
      bond_constructions.emplace(
        molecule, reroute_begin_atom_idx, reroute_end_atom_idx, bond_type);
    };
    bond_destructions.emplace(begin_atom_idx, end_atom_idx);
  };

  Type GetType() const { return type; };
};


boost::dynamic_bitset<> IndexVectorToBitset(
  const std::vector<std::size_t>& indices,
  std::size_t size) {
  boost::dynamic_bitset<> bits (size);
  for (std::size_t idx : indices) {
    bits.set(idx);
  };
  return bits;
};


class SubgraphConstruction : public TopologicalPerturbation {
  static const Type type = Type::SubgraphConstruction_t;
  std::unordered_map<AtomIdx, AtomIdx> atom_mapping;

public:
  SubgraphConstruction(
    const RDKit::ROMol& target,
    const RDKit::ROMol& source,
    const boost::dynamic_bitset<>& subgraph_atom_mask) {
    // Adds an induced subgraph of the source molecule to the target molecule.
    // The subgraph is specified by a set of atoms, and bonds between said
    // atoms are considered part of the subgraph.
    std::size_t n_atoms = target.getNumAtoms();
    std::size_t max_atom_idx = n_atoms ? n_atoms - 1 : 0;
    Tag max_atom_tag = GetMaxAtomTag(target);
    Tag max_bond_tag = GetMaxBondTag(target);
    boost::dynamic_bitset<> bonds_mask (source.getNumBonds());
    for (AtomIdx atom_idx = subgraph_atom_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = subgraph_atom_mask.find_next(atom_idx)) {
      const RDKit::Atom* atom = source.getAtomWithIdx(atom_idx);
      atom_mapping[atom_idx] = ++max_atom_idx;
      atom_constructions.emplace(++max_atom_tag,
        atom->getAtomicNum(),
        atom->getFormalCharge(),
        atom->getNumExplicitHs());
      for (const RDKit::Bond* bond : source.atomBonds(atom)) {
        if (subgraph_atom_mask[bond->getBeginAtomIdx()] &&
          subgraph_atom_mask[bond->getEndAtomIdx()]) {
          bonds_mask.set(bond->getIdx());
        };
      };
    };
    for (BondIdx bond_idx = bonds_mask.find_first();
      bond_idx != boost::dynamic_bitset<>::npos;
      bond_idx = bonds_mask.find_next(bond_idx)) {
      const RDKit::Bond* bond = source.getBondWithIdx(bond_idx);
      bond_constructions.emplace(++max_bond_tag,
        atom_mapping[bond->getBeginAtomIdx()],
        atom_mapping[bond->getEndAtomIdx()],
        bond->getBondType());
    };
  };

  SubgraphConstruction(
    const RDKit::ROMol& target,
    const RDKit::ROMol& source,
    const std::vector<AtomIdx>& subgraph_atom_indices) :
    SubgraphConstruction(target, source, IndexVectorToBitset(
      subgraph_atom_indices, source.getNumAtoms())) {};
  
  const std::unordered_map<AtomIdx, AtomIdx>& GetAtomMapping() const {
    return atom_mapping;
  };

  Type GetType() const { return type; };
};


class SubgraphDestruction : public TopologicalPerturbation {
  static const Type type = Type::SubgraphDestruction_t;

public:
  SubgraphDestruction(
    const boost::dynamic_bitset<>& subgraph_atom_mask) {
    for (AtomIdx atom_idx = subgraph_atom_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = subgraph_atom_mask.find_next(atom_idx)) {
      atom_destructions.emplace(atom_idx);
    };
  };

  SubgraphDestruction(
    const RDKit::ROMol& molecule,
    const std::vector<AtomIdx>& subgraph_atom_indices) :
    SubgraphDestruction(IndexVectorToBitset(
      subgraph_atom_indices, molecule.getNumAtoms())) {};
  
  Type GetType() const { return type; };
};


class SubgraphPerturbation : public TopologicalPerturbation {
  static const Type type = Type::SubgraphPerturbation_t;
  SubgraphConstruction subgraph_construction;
  SubgraphDestruction subgraph_destruction;

public:
  SubgraphPerturbation(
    const RDKit::ROMol& target,
    const RDKit::ROMol& source,
    const boost::dynamic_bitset<>& constructed_subgraph,
    const boost::dynamic_bitset<>& destroyed_subgraph) :
    subgraph_construction(target, source, constructed_subgraph),
    subgraph_destruction(destroyed_subgraph) {};

  SubgraphPerturbation(
    const RDKit::ROMol& target,
    const RDKit::ROMol& source,
    const std::vector<AtomIdx>& constructed_subgraph_atom_indices,
    const std::vector<AtomIdx>& destroyed_subgraph_atom_indices) :
    subgraph_construction(target, source, constructed_subgraph_atom_indices),
    subgraph_destruction(target, destroyed_subgraph_atom_indices) {};

  void operator()(RDKit::RWMol& molecule) const {
    subgraph_construction(molecule);
    TopologicalPerturbation::operator()(molecule);
    subgraph_destruction(molecule);
  };

  void ProjectMolecularGraph(MolecularGraphProjection& projection) const {
    subgraph_construction.ProjectMolecularGraph(projection);
    TopologicalPerturbation::ProjectMolecularGraph(projection);
    subgraph_destruction.ProjectMolecularGraph(projection);
  };

  std::size_t ID() const {
    std::size_t id = TopologicalPerturbation::ID();
    boost::hash_combine(id, subgraph_construction.ID());
    boost::hash_combine(id, subgraph_destruction.ID());
    return id;
  };

  const SubgraphConstruction& GetSubgraphConstruction() const {
    return subgraph_construction;
  };

  const SubgraphDestruction& GetSubgraphDestruction() const {
    return subgraph_destruction;
  };

  Type GetType() const { return type; };
};

#endif // !_MOLECULAR_PERTURBATIONS_HPP_
