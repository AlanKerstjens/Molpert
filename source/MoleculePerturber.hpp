#pragma once
#ifndef _MOLECULE_PERTURBER_HPP_
#define _MOLECULE_PERTURBER_HPP_

#include "Combinatorics.hpp"
#include "RandomSampling.hpp"
#include "MolecularGraphSearch.hpp"
#include "MolecularPerturbations.hpp"
#include "MolecularConstraints.hpp"
#include <boost/core/null_deleter.hpp>

class MolecularPerturbationQueue {
  std::deque<std::shared_ptr<MolecularPerturbation>> queue;
  std::unordered_set<std::size_t> queued_perturbation_ids;

public:
  bool push(const std::shared_ptr<MolecularPerturbation>& perturbation) {
    std::size_t id = perturbation->ID();
    // If a perturbation was ever queued before it can't be queued again.
    if (queued_perturbation_ids.contains(id)) {
      return false;
    };
    queue.push_back(perturbation);
    return true;
  };

  void pop() {
    queue.pop_front();
  };

  void shuffle(std::mt19937& prng) {
    std::shuffle(queue.begin(), queue.end(), prng);
  };

  std::shared_ptr<MolecularPerturbation> front() const {
    return queue.front();
  };

  bool empty() const {
    return queue.empty();
  };

  std::size_t size() const {
    return queue.size();
  };
};


class MoleculePerturber {
public:
  MolecularPerturbation::TypeMask perturbation_types = 0xFFFF;
  std::array<double, MolecularPerturbation::n_types> perturbation_types_weights
    {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};

  // Allowed values for atom / bond properties.
  std::vector<std::uint8_t> atomic_numbers {6, 7, 8, 9, 15, 16, 17, 35, 53};
  std::vector<std::uint8_t> ring_atomic_numbers {6, 7, 8, 16};
  std::vector<std::int8_t> formal_charges {0, 1, -1};
  std::vector<std::uint8_t> n_explicit_hydrogens {0, 1};
  std::vector<RDKit::Bond::BondType> bond_types
    {RDKit::Bond::SINGLE, RDKit::Bond::DOUBLE, RDKit::Bond::TRIPLE};
  std::vector<RDKit::Bond::BondType> ring_bond_types
    {RDKit::Bond::SINGLE, RDKit::Bond::DOUBLE};

  // Sampling weights for each of the allowed atom / bond properties.
  std::vector<double> atomic_numbers_weights {1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> ring_atomic_numbers_weights {1, 1, 1, 1};
  std::vector<double> formal_charges_weights {1, 1, 1};
  std::vector<double> n_explicit_hydrogens_weights {1, 1};
  std::vector<double> bond_types_weights {1, 1, 1};
  std::vector<double> ring_bond_types_weights {1, 1};

  // Default values used for certain MolecularPerturbations. The purpose of
  // having these default values is to reduce the number of outputs of the
  // deterministic overloads.
  std::uint8_t default_atomic_number = 6;
  std::int8_t default_formal_charge = 0;
  std::uint8_t default_n_explicit_hydrogens = 0;
  RDKit::Bond::BondType default_bond_type = RDKit::Bond::SINGLE;

  // General settings
  static const unsigned max_unsigned = std::numeric_limits<unsigned>::max();
  bool allow_disconnections = false;
  bool assess_connectivity_with_sssr = false;
  bool assess_distances_with_distance_matrix = false;

  // Atomic number change settings
  bool cyclicity_based_atomic_numbers = false;
  // Bond type change settings
  bool cyclicity_based_bond_types = false;
  // These settings don't apply to atom/bond insertions because it is difficult
  // to predict a priori if the newly added atom or bond will be part of a cycle
  // or not.

  // Atom insertion settings
  unsigned atom_insertion_max_n_neighbors = 3;
  bool atom_insertion_drop_an_atom = true;
  bool atom_insertion_dropping_is_optional = false; // Only when above is true
  // Atom insertion settings (central atom overloads)
  bool atom_insertion_only_pericentral_atoms_as_neighbor_candidates = true;
  unsigned atom_insertion_min_distance_neighbor = 1;
  unsigned atom_insertion_max_distance_neighbor = max_unsigned;
  // Atom insertion settings (stochastic overloads)
  bool atom_insertion_randomize_atomic_number = true;
  bool atom_insertion_randomize_formal_charge = false;
  bool atom_insertion_randomize_n_explicit_hydrogens = false;
  bool atom_insertion_randomize_bond_types = true;
  // Atom insertion settings (deterministic overloads)
  bool atom_insertion_iterate_atomic_numbers = false;
  bool atom_insertion_iterate_formal_charges = false;
  bool atom_insertion_iterate_n_explicit_hydrogens = false;
  bool atom_insertion_iterate_bond_types = false;

  // Atom deletion settings
  // Deleting without reconnections (when possible) has priority, even if
  // reconnections are allowed.
  bool atom_deletion_allow_reconnections = true;
  bool atom_deletion_only_consider_neighbors_for_reconnection = true;
  bool atom_deletion_preserve_bond_types_during_reconnection = true;
  unsigned atom_deletion_min_distance_reconnection_atom = 1;
  unsigned atom_deletion_max_distance_reconnection_atom = max_unsigned;

  // Bond insertion settings
  unsigned bond_insertion_min_distance_partner = 2;
  unsigned bond_insertion_max_distance_partner = max_unsigned;
  unsigned bond_insertion_max_atom_n_rings_membership = max_unsigned;
  // Bond insertion settings (stochastic overloads)
  bool bond_insertion_randomize_bond_type = true;
  // Bond insertion settings (deterministic overloads)
  bool bond_insertion_iterate_bond_types = false;

  // Bond deletion settings
  bool bond_deletion_allow_reroutes = true;
  bool bond_deletion_preserve_bond_types_during_reroute = true;
  unsigned bond_deletion_min_distance_reroute = 2;
  unsigned bond_deletion_max_distance_reroute = 2;

private:
  std::vector<std::size_t> ShuffledAtomIndices(
    const RDKit::ROMol& molecule, std::mt19937& prng) const {
    std::vector<std::size_t> atom_indices (molecule.getNumAtoms());
    std::iota(atom_indices.begin(), atom_indices.end(), 0);
    std::shuffle(atom_indices.begin(), atom_indices.end(), prng);
    return atom_indices;
  };

  std::vector<std::size_t> ShuffledBondIndices(
    const RDKit::ROMol& molecule, std::mt19937& prng) const {
    std::vector<std::size_t> bond_indices (molecule.getNumBonds());
    std::iota(bond_indices.begin(), bond_indices.end(), 0);
    std::shuffle(bond_indices.begin(), bond_indices.end(), prng);
    return bond_indices;
  };

  bool CanDeleteAtomWithoutReconnection(
    const RDKit::ROMol& molecule, AtomIdx atom_idx) const {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    return allow_disconnections || atom->getDegree() <= 1 ||
      ConnectedAfterAtomDeletion(molecule, atom_idx);
  };

  bool CanDeleteBondWithoutReconnection(
    const RDKit::ROMol& molecule, BondIdx bond_idx) const {
    return allow_disconnections || ConnectedAfterBondDeletion(
      molecule, bond_idx, assess_connectivity_with_sssr);
  };

  std::vector<std::size_t> ExtractBits(
    const boost::dynamic_bitset<>& bits) const {
    std::vector<std::size_t> extracted_bit_indices;
    extracted_bit_indices.reserve(bits.count());
    for (std::size_t bit_idx = bits.find_first();
      bit_idx != boost::dynamic_bitset<>::npos;
      bit_idx = bits.find_next(bit_idx)) {
      extracted_bit_indices.push_back(bit_idx);
    };
    return extracted_bit_indices;
  };

  std::vector<std::size_t> ExtractBits(
    const boost::dynamic_bitset<>& bits,
    const Mask& mask) const {
    std::size_t n_mask_bits = mask.size();
    assert(bits.count() >= n_mask_bits);
    std::vector<std::size_t> extracted_bit_indices;
    extracted_bit_indices.reserve(n_mask_bits);
    std::size_t bit_idx = bits.find_first();
    for (std::size_t mask_idx = 0; mask_idx < n_mask_bits; ++mask_idx) {
      if (mask[mask_idx]) {
        extracted_bit_indices.push_back(bit_idx);
      };
      bit_idx = bits.find_next(bit_idx);
    };
    return extracted_bit_indices;
  };

  boost::dynamic_bitset<> AllAtomsMask(const RDKit::ROMol& molecule) const {
    boost::dynamic_bitset<> atom_mask (molecule.getNumAtoms());
    atom_mask.set();
    return atom_mask;
  };

  boost::dynamic_bitset<> AtomPeripheryMask(
    const RDKit::ROMol& molecule,
    AtomIdx central_atom_idx,
    bool include_central_atom = true) const {
    boost::dynamic_bitset<> atom_mask (molecule.getNumAtoms());
    if (include_central_atom) {
      atom_mask.set(central_atom_idx);
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(central_atom_idx);
    for (const RDKit::Atom* neighbor : molecule.atomNeighbors(atom)) {
      atom_mask.set(neighbor->getIdx());
    };
    return atom_mask;
  };

  boost::dynamic_bitset<> AtomInsertionNeighborMask(
    const RDKit::ROMol& molecule,
    AtomIdx central_atom_idx) const {
    if (atom_insertion_only_pericentral_atoms_as_neighbor_candidates) {
      return AtomPeripheryMask(molecule, central_atom_idx, false);
    };
    boost::dynamic_bitset<> atom_mask (molecule.getNumAtoms());
    if (atom_insertion_min_distance_neighbor > 1 ||
      atom_insertion_max_distance_neighbor < max_unsigned) {
      std::vector<std::size_t> candidate_neighbor_indices = AtomsWithinDistance(
        molecule, central_atom_idx,
        atom_insertion_min_distance_neighbor,
        atom_insertion_max_distance_neighbor,
        assess_distances_with_distance_matrix);
      for (std::size_t neighbor_idx : candidate_neighbor_indices) {
        atom_mask.set(neighbor_idx);
      };
    } else {
      atom_mask.set();
    };
    atom_mask.reset(central_atom_idx);
    return atom_mask;
  };

  boost::dynamic_bitset<> AtomDeletionReconnectionMask(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx) const {
    if (!atom_deletion_allow_reconnections) {
      return boost::dynamic_bitset<>();
    };
    if (atom_deletion_only_consider_neighbors_for_reconnection) {
      return AtomPeripheryMask(molecule, atom_idx, false);
    };
    boost::dynamic_bitset<> atom_mask (molecule.getNumAtoms());
    if (atom_deletion_min_distance_reconnection_atom > 1 ||
      atom_deletion_max_distance_reconnection_atom < max_unsigned) {
      std::vector<std::size_t> candidate_reconnection_indices =
        AtomsWithinDistance(molecule, atom_idx,
          atom_deletion_min_distance_reconnection_atom,
          atom_deletion_max_distance_reconnection_atom,
          assess_distances_with_distance_matrix);
      for (std::size_t reconnection_idx : candidate_reconnection_indices) {
        atom_mask.set(reconnection_idx);
      };
    } else {
      atom_mask.set();
    };
    atom_mask.reset(atom_idx);
    return atom_mask;
  };

  boost::dynamic_bitset<> RingMembershipMask(
    const RDKit::ROMol& molecule,
    unsigned min_n_rings = 0,
    unsigned max_n_rings = 1) const {
    const RDKit::RingInfo* ring_info = molecule.getRingInfo();
    if (!ring_info->isInitialized()) {
      RDKit::MolOps::findSSSR(molecule);
    };
    std::size_t n_atoms = molecule.getNumAtoms();
    boost::dynamic_bitset<> atom_mask (n_atoms);
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      unsigned n_rings = ring_info->numAtomRings(atom_idx);
      if (n_rings >= min_n_rings && n_rings <= max_n_rings) {
        atom_mask.set(atom_idx);
      };
    };
    return atom_mask;
  };

  boost::dynamic_bitset<> BondingPartnerMask(
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx) const {
    boost::dynamic_bitset<> atom_mask (molecule.getNumAtoms());
    if (bond_insertion_min_distance_partner > 2 ||
      bond_insertion_max_distance_partner < max_unsigned) {
      std::vector<std::size_t> candidate_partner_indices = AtomsWithinDistance(
        molecule, begin_atom_idx,
        bond_insertion_min_distance_partner,
        bond_insertion_max_distance_partner,
        assess_distances_with_distance_matrix);
      for (std::size_t partner_idx : candidate_partner_indices) {
        atom_mask.set(partner_idx);
      };
    } else {
      atom_mask.set();
    };
    if (bond_insertion_max_atom_n_rings_membership < max_unsigned) {
      atom_mask &= RingMembershipMask(
        molecule, 0, bond_insertion_max_atom_n_rings_membership - 1);
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(begin_atom_idx);
    for (const RDKit::Atom* neighbor : molecule.atomNeighbors(atom)) {
      atom_mask.reset(neighbor->getIdx());
    };
    atom_mask.reset(begin_atom_idx);
    return atom_mask;
  };

  MolecularPerturbation::TypeMask EligiblePerturbationTypes(
    const RDKit::ROMol& molecule) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    std::size_t n_bonds = molecule.getNumBonds();
    MolecularPerturbation::TypeMask eligible_types;
    // Atom insertions are always possible.
    if (perturbation_types[MolecularPerturbation::AtomInsertion_t]) {
      eligible_types.set(MolecularPerturbation::AtomInsertion_t);
    };
    if (n_atoms == 0) {
      return eligible_types;
    };
    // If we have at least one atom we can change its properties.
    if (n_atoms > 0) {
      if (perturbation_types[MolecularPerturbation::AtomicNumberChange_t]) {
        eligible_types.set(MolecularPerturbation::AtomicNumberChange_t);
      };
      if (perturbation_types[MolecularPerturbation::FormalChargeChange_t]) {
        eligible_types.set(MolecularPerturbation::FormalChargeChange_t);
      };
      if (perturbation_types[MolecularPerturbation::ExplicitHydrogensChange_t]) {
        eligible_types.set(MolecularPerturbation::ExplicitHydrogensChange_t);
      };
    };
    // If we have at least 2 atoms we can delete one or bond them together.
    if (n_atoms > 1) {
      if (perturbation_types[MolecularPerturbation::AtomDeletion_t]) {
        eligible_types.set(MolecularPerturbation::AtomDeletion_t);
      };
      if (perturbation_types[MolecularPerturbation::BondInsertion_t]) {
        eligible_types.set(MolecularPerturbation::BondInsertion_t);
      };
    };
    // If we have at least 1 bond we can change its properties or delete it.
    if (n_bonds > 0) {
      if (perturbation_types[MolecularPerturbation::BondTypeChange_t]) {
        eligible_types.set(MolecularPerturbation::BondTypeChange_t);
      };
      // If we wish to preserve a connected molecule a bond can only be deleted
      // when it's part of a cycle. Molecules with less than 3 atoms are
      // guaranteed to not have cycles.
      if (n_atoms > 2 || allow_disconnections) {
        if (perturbation_types[MolecularPerturbation::BondDeletion_t]) {
          eligible_types.set(MolecularPerturbation::BondDeletion_t);
        };
      };
    };
    return eligible_types;
  };

  std::vector<MolecularPerturbation::Type> ShuffledEligiblePerturbationTypes(
    const RDKit::ROMol& molecule, std::mt19937& prng) const {
    MolecularPerturbation::TypeMask eligible_perturbation_types =
      EligiblePerturbationTypes(molecule);
    std::vector<MolecularPerturbation::Type> types;
    std::vector<double> weights;
    types.reserve(MolecularPerturbation::n_types);
    weights.reserve(MolecularPerturbation::n_types);
    for (std::size_t type = 0; type < MolecularPerturbation::n_types; ++type) {
      if (eligible_perturbation_types[type]) {
        types.push_back(MolecularPerturbation::Type(type));
        weights.push_back(perturbation_types_weights[type]);
      };
    };
    WeightedShuffleInPlace(types, weights, prng);
    return types;
  };

public:
  MoleculePerturber(
    bool use_chembl_distribution = true,
    bool use_aromatic_bonds = false,
    bool acyclic_can_be_aromatic = false,
    bool cyclicity_based_atomic_numbers = false,
    bool cyclicity_based_bond_types = false) :
    cyclicity_based_atomic_numbers(cyclicity_based_atomic_numbers), 
    cyclicity_based_bond_types(cyclicity_based_bond_types) {
    if (use_chembl_distribution) {
      SetChEMBLProperties();
    };
    if (use_aromatic_bonds) {
      SetAromaticProperties(acyclic_can_be_aromatic, use_chembl_distribution);
    };
  };

  void SetChEMBLProperties() {
    // The default value for atom and bond properties is the most frequent value.
    // The probability of changing a property is proportional to how often said
    // property deviates from its default value. We multiply the latter value by
    // the number of (property) perturbations to sample from, that is, 4.
    default_atomic_number = 6;
    default_formal_charge = 0;
    default_n_explicit_hydrogens = 0;
    default_bond_type = RDKit::Bond::SINGLE;
    perturbation_types_weights[
      MolecularPerturbation::Type::AtomicNumberChange_t] = 0.2620 * 4;
    perturbation_types_weights[
      MolecularPerturbation::Type::FormalChargeChange_t] = 0.0042 * 4;
    perturbation_types_weights[
      MolecularPerturbation::Type::ExplicitHydrogensChange_t] = 0.0312 * 4;
    perturbation_types_weights[
      MolecularPerturbation::Type::BondTypeChange_t] = 0.5103 * 4;
    // The values listed below are allowed. Other values occur with negligible
    // (< 0.01%) frequencies in drug-like molecules (ChEMBL) and are ignored.
    // For reference, B and Si have about 0.006% and As has about 0.004%.
    atomic_numbers = {6, 8, 7, 9, 16, 17, 35, 15, 53};
    ring_atomic_numbers = {6, 7, 8, 16};
    formal_charges = {0, 1, -1};
    n_explicit_hydrogens = {0, 1};
    bond_types = {RDKit::Bond::SINGLE, RDKit::Bond::DOUBLE, RDKit::Bond::TRIPLE};
    ring_bond_types = {RDKit::Bond::SINGLE, RDKit::Bond::DOUBLE};
    if (cyclicity_based_atomic_numbers) {
      atomic_numbers_weights =
        {10010308, 5369046, 2200788, 719854, 414403, 413065, 85966, 36870, 11771};
    } else {
      atomic_numbers_weights =
        {39042049, 6055511, 5821996, 719854, 705007, 413065, 85966, 39562, 11815};
    };
    ring_atomic_numbers_weights = {29031741, 3621208, 686465, 290604};
    formal_charges_weights = {52682232, 121749, 100139};
    n_explicit_hydrogens_weights = {51251210, 1652859};
    if (cyclicity_based_bond_types) {
      bond_types_weights = {18642643, 3488616, 120889};
    } else {
      bond_types_weights = {40892139, 16473946, 121206};
    };
    ring_bond_types_weights = {22249495, 12985329};
  };

  void SetAromaticProperties(
    bool acyclic_can_be_aromatic = false,
    bool use_chembl_distribution = true) {
    if (acyclic_can_be_aromatic) {
      bond_types = {
        RDKit::Bond::SINGLE, RDKit::Bond::DOUBLE, 
        RDKit::Bond::TRIPLE, RDKit::Bond::AROMATIC};
      if (use_chembl_distribution) {
        bond_types_weights = {28149153, 3730960, 121206, 25485971};
      } else {
        bond_types_weights = {1, 1, 1, 1};
      };
    };
    ring_bond_types = {
        RDKit::Bond::SINGLE, RDKit::Bond::DOUBLE, RDKit::Bond::AROMATIC};
    if (use_chembl_distribution) {
      ring_bond_types_weights = {9506510, 242344, 25485971};
    } else {
      ring_bond_types_weights = {1, 1, 1};
    };
  };

  std::shared_ptr<AtomicNumberChange> ChangeAtomicNumber(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    const std::vector<std::uint8_t>* values = &atomic_numbers;
    const std::vector<double>* weights = &atomic_numbers_weights;
    if (cyclicity_based_atomic_numbers && AtomIsInRing(molecule, atom_idx)) {
      values = &ring_atomic_numbers;
      weights = &ring_atomic_numbers_weights;
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    std::uint8_t current_atomic_number = atom->getAtomicNum();
    for (std::uint8_t atomic_number : WeightedShuffle(*values, *weights, prng)) {
      if (atomic_number == current_atomic_number) {
        continue;
      };
      auto perturbation = std::make_shared<AtomicNumberChange>(
        atom_idx, atomic_number);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<AtomicNumberChange> ChangeAtomicNumber(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    for (std::size_t atom_idx : ShuffledAtomIndices(molecule, prng)) {
      auto perturbation = ChangeAtomicNumber(
        molecule, atom_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void AtomicNumberChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<std::uint8_t>* allowed_values = nullptr) const {
    const std::vector<std::uint8_t>* values = &atomic_numbers;
    if (allowed_values) {
      values = allowed_values;
    } else if (
      cyclicity_based_atomic_numbers && AtomIsInRing(molecule, atom_idx)) {
      values = &ring_atomic_numbers;
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    std::uint8_t current_atomic_number = atom->getAtomicNum();
    for (std::uint8_t atomic_number : *values) {
      if (atomic_number == current_atomic_number) {
        continue;
      };
      auto perturbation = std::make_shared<AtomicNumberChange>(
        atom_idx, atomic_number);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void AtomicNumberChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<std::uint8_t>* allowed_values = nullptr) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      AtomicNumberChanges(
        queue, molecule, atom_idx, constraints, allowed_values);
    };
  };


  std::shared_ptr<FormalChargeChange> ChangeFormalCharge(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    std::int8_t current_formal_charge = atom->getFormalCharge();
    for (std::int8_t formal_charge :
      WeightedShuffle(formal_charges, formal_charges_weights, prng)) {
      if (formal_charge == current_formal_charge) {
        continue;
      };
      auto perturbation = std::make_shared<FormalChargeChange>(
        atom_idx, formal_charge);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<FormalChargeChange> ChangeFormalCharge(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    for (std::size_t atom_idx : ShuffledAtomIndices(molecule, prng)) {
      auto perturbation = ChangeFormalCharge(
        molecule, atom_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void FormalChargeChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<std::int8_t>* allowed_values = nullptr) const {
    const std::vector<std::int8_t>* values = &formal_charges;
    if (allowed_values) {
      values = allowed_values;
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    std::int8_t current_formal_charge = atom->getFormalCharge();
    for (std::int8_t formal_charge : *values) {
      if (formal_charge == current_formal_charge) {
        continue;
      };
      auto perturbation = std::make_shared<FormalChargeChange>(
        atom_idx, formal_charge);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void FormalChargeChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<std::int8_t>* allowed_values = nullptr) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      FormalChargeChanges(
        queue, molecule, atom_idx, constraints, allowed_values);
    };
  };


  std::shared_ptr<ExplicitHydrogensChange> ChangeExplicitHydrogens(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    std::uint8_t current_explicit_hydrogens = atom->getNumExplicitHs();
    for (std::uint8_t explicit_hydrogens :
      WeightedShuffle(n_explicit_hydrogens, n_explicit_hydrogens_weights, prng)) {
      if (explicit_hydrogens == current_explicit_hydrogens) {
        continue;
      };
      auto perturbation = std::make_shared<ExplicitHydrogensChange>(
        atom_idx, explicit_hydrogens);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<ExplicitHydrogensChange> ChangeExplicitHydrogens(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    for (std::size_t atom_idx : ShuffledAtomIndices(molecule, prng)) {
      auto perturbation = ChangeExplicitHydrogens(
        molecule, atom_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void ExplicitHydrogenChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<std::uint8_t>* allowed_values = nullptr) const {
    const std::vector<std::uint8_t>* values = &n_explicit_hydrogens;
    if (allowed_values) {
      values = allowed_values;
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    std::uint8_t current_explicit_hydrogens = atom->getNumExplicitHs();
    for (std::uint8_t explicit_hydrogens : *values) {
      if (explicit_hydrogens == current_explicit_hydrogens) {
        continue;
      };
      auto perturbation = std::make_shared<ExplicitHydrogensChange>(
        atom_idx, explicit_hydrogens);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void ExplicitHydrogenChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<std::uint8_t>* allowed_values = nullptr) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      ExplicitHydrogenChanges(
        queue, molecule, atom_idx, constraints, allowed_values);
    };
  };


  std::shared_ptr<BondTypeChange> ChangeBondType(
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    const std::vector<RDKit::Bond::BondType>* values = &bond_types;
    const std::vector<double>* weights = &bond_types_weights;
    if (cyclicity_based_bond_types && BondIsInRing(molecule, bond_idx)) {
      values = &ring_bond_types;
      weights = &ring_bond_types_weights;
    };
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    RDKit::Bond::BondType current_bond_type = bond->getBondType();
    for (RDKit::Bond::BondType bond_type :
      WeightedShuffle(*values, *weights, prng)) {
      if (bond_type == current_bond_type) {
        continue;
      };
      auto perturbation = std::make_shared<BondTypeChange>(
        molecule, bond_idx, bond_type);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<BondTypeChange> ChangeBondType(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    for (std::size_t bond_idx : ShuffledBondIndices(molecule, prng)) {
      auto perturbation = ChangeBondType(molecule, bond_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void BondTypeChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<RDKit::Bond::BondType>* allowed_values = nullptr) const {
    const std::vector<RDKit::Bond::BondType>* values = &bond_types;
    if (allowed_values) {
      values = allowed_values;
    } else if (cyclicity_based_bond_types && BondIsInRing(molecule, bond_idx)) {
      values = &ring_bond_types;
    };
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    RDKit::Bond::BondType current_bond_type = bond->getBondType();
    for (RDKit::Bond::BondType bond_type : *values) {
      if (bond_type == current_bond_type) {
        continue;
      };
      auto perturbation = std::make_shared<BondTypeChange>(
        molecule, bond_idx, bond_type);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void BondTypeChanges(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr,
    const std::vector<RDKit::Bond::BondType>* allowed_values = nullptr) const {
    std::size_t n_bonds = molecule.getNumBonds();
    for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
      BondTypeChanges(queue, molecule, bond_idx, constraints, allowed_values);
    };
  };


  std::shared_ptr<AtomInsertion> InsertAtom(
    const RDKit::ROMol& molecule,
    const std::vector<AtomIdx>& neighbor_atom_indices,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    // Determine possible atomic numbers for the inserted atom.
    std::vector<std::uint8_t> candidate_atomic_numbers {default_atomic_number};
    if (atom_insertion_randomize_atomic_number) {
      candidate_atomic_numbers = WeightedShuffle(
        atomic_numbers, atomic_numbers_weights, prng);
    };
    // Determine possible formal charges for the inserted atom.
    std::vector<std::int8_t> candidate_formal_charges {default_formal_charge};
    if (atom_insertion_randomize_formal_charge) {
      candidate_formal_charges = WeightedShuffle(
        formal_charges, formal_charges_weights, prng);
    };
    // Determine possible numbers of explicit hydrogens for the inserted atom.
    std::vector<std::uint8_t> candidate_n_explicit_hydrogens {
      default_n_explicit_hydrogens};
    if (atom_insertion_randomize_n_explicit_hydrogens) {
      candidate_n_explicit_hydrogens = WeightedShuffle(
        n_explicit_hydrogens, n_explicit_hydrogens_weights, prng);
    };
    // Determine the candidate dropped neighbor atoms.
    std::vector<int> candidate_dropped_atom_indices {-1};
    if (atom_insertion_drop_an_atom && !neighbor_atom_indices.empty()) {
      if (!atom_insertion_dropping_is_optional) {
        candidate_dropped_atom_indices.clear();
      };
      candidate_dropped_atom_indices.insert(
        candidate_dropped_atom_indices.end(),
        neighbor_atom_indices.cbegin(), neighbor_atom_indices.cend());
      std::shuffle(candidate_dropped_atom_indices.begin(),
        candidate_dropped_atom_indices.end(), prng);
    };
    // Determine possible bond types for the formed bonds.
    std::size_t n_neighbors = neighbor_atom_indices.size();
    std::vector<std::vector<RDKit::Bond::BondType>>
      candidate_bond_types (n_neighbors, {default_bond_type});
    if (atom_insertion_randomize_bond_types) {
      for (std::size_t i = 0; i < n_neighbors; ++i) {
        candidate_bond_types[i] = WeightedShuffle(
          bond_types, bond_types_weights, prng);
      };
    };
    // Iterate over all combinations of dropped neighbor atoms, atom
    // properties and bond properties.
    // The order of the wheels in the odometer isn't capricious. The first
    // wheels are the most dynamic while the last wheels are the most static.
    // Properties one wishes to preserve the most should be on static wheels.
    // In our case we prioritize iterating through bond types, as this is most
    // likely to accomodate constraints with a small number of iterations.
    Odometer odometer;
    if (n_neighbors) {
      odometer.AddWheels(candidate_bond_types.front().size() - 1, n_neighbors);
    };
    odometer.AddWheels(candidate_atomic_numbers.size() - 1);
    odometer.AddWheels(candidate_formal_charges.size() - 1);
    odometer.AddWheels(candidate_n_explicit_hydrogens.size() - 1);
    odometer.AddWheels(candidate_dropped_atom_indices.size() - 1);
    while (!odometer.revs()) {
      std::vector<RDKit::Bond::BondType> neighbor_bond_types (n_neighbors);
      for (std::size_t i = 0; i < n_neighbors; ++i) {
        neighbor_bond_types[i] = candidate_bond_types[i][odometer[i]];
      };
      auto perturbation = std::make_shared<AtomInsertion>(
        molecule,
        neighbor_atom_indices,
        neighbor_bond_types,
        candidate_atomic_numbers[odometer[n_neighbors]],
        candidate_formal_charges[odometer[n_neighbors + 1]],
        candidate_n_explicit_hydrogens[odometer[n_neighbors + 2]],
        candidate_dropped_atom_indices[odometer[n_neighbors + 3]]);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        ++odometer;
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<AtomInsertion> InsertAtom(
    const RDKit::ROMol& molecule,
    const boost::dynamic_bitset<>& candidate_neighbor_atoms_mask,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    // Iterate over all combinations of candidate neighbors in a random order.
    RandomCombinator combinator (candidate_neighbor_atoms_mask.count(),
      allow_disconnections ? 0 : 1, atom_insertion_max_n_neighbors, prng);
    while (!combinator.revs) {
      auto perturbation = InsertAtom(molecule,
        ExtractBits(candidate_neighbor_atoms_mask, combinator.mask),
        prng, constraints);
      if (perturbation) {
        return perturbation;
      };
      ++combinator;
    };
    return nullptr;
  };

  std::shared_ptr<AtomInsertion> InsertAtom(
    const RDKit::ROMol& molecule,
    AtomIdx central_atom_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    boost::dynamic_bitset<> candidate_neighbor_atoms_mask =
      AtomInsertionNeighborMask(molecule, central_atom_idx);
    RandomCombinator combinator (candidate_neighbor_atoms_mask.count(),
      0, atom_insertion_max_n_neighbors, prng);
    while (!combinator.revs) {
      std::vector<std::size_t> neighbor_atom_indices = ExtractBits(
        candidate_neighbor_atoms_mask, combinator.mask);
      neighbor_atom_indices.push_back(central_atom_idx);
      auto perturbation = InsertAtom(
        molecule, neighbor_atom_indices, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
      ++combinator;
    };
    return nullptr;
  };

  std::shared_ptr<AtomInsertion> InsertAtom(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    if (!molecule.getNumAtoms()) {
      auto perturbation = InsertAtom(
        molecule, std::vector<AtomIdx>(), prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    for (std::size_t central_atom_idx : ShuffledAtomIndices(molecule, prng)) {
      auto perturbation = InsertAtom(
        molecule, central_atom_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void AtomInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const std::vector<AtomIdx>& neighbor_atom_indices,
    const MolecularConstraints* constraints = nullptr) const {
    std::shared_ptr<const std::vector<std::uint8_t>> candidate_atomic_numbers;
    if (atom_insertion_iterate_atomic_numbers) {
      candidate_atomic_numbers.reset(&atomic_numbers, boost::null_deleter());
    } else {
      candidate_atomic_numbers.reset(
        new std::vector<std::uint8_t> {default_atomic_number});
    };
    std::shared_ptr<const std::vector<std::int8_t>> candidate_formal_charges;
    if (atom_insertion_iterate_formal_charges) {
      candidate_formal_charges.reset(&formal_charges, boost::null_deleter());
    } else {
      candidate_formal_charges.reset(
        new std::vector<std::int8_t> {default_formal_charge});
    };
    std::shared_ptr<const std::vector<std::uint8_t>>
      candidate_n_explicit_hydrogens;
    if (atom_insertion_iterate_n_explicit_hydrogens) {
      candidate_n_explicit_hydrogens.reset(
        &n_explicit_hydrogens, boost::null_deleter());
    } else {
      candidate_n_explicit_hydrogens.reset(
        new std::vector<std::uint8_t> {default_n_explicit_hydrogens});
    };
    std::shared_ptr<const std::vector<RDKit::Bond::BondType>> candidate_bond_types;
    if (atom_insertion_iterate_bond_types) {
      candidate_bond_types.reset(&bond_types, boost::null_deleter());
    } else {
      candidate_bond_types.reset(
        new std::vector<RDKit::Bond::BondType> {default_bond_type});
    };
    std::vector<int> candidate_dropped_atom_indices {-1};
    if (atom_insertion_drop_an_atom && !neighbor_atom_indices.empty()) {
      if (!atom_insertion_dropping_is_optional) {
        candidate_dropped_atom_indices.clear();
      };
      candidate_dropped_atom_indices.insert(
        candidate_dropped_atom_indices.end(),
        neighbor_atom_indices.cbegin(), neighbor_atom_indices.cend());
    };
    std::size_t n_neighbors = neighbor_atom_indices.size();
    Odometer odometer (candidate_bond_types->size() - 1, n_neighbors);
    odometer.AddWheels(candidate_atomic_numbers->size() - 1);
    odometer.AddWheels(candidate_formal_charges->size() - 1);
    odometer.AddWheels(candidate_n_explicit_hydrogens->size() - 1);
    odometer.AddWheels(candidate_dropped_atom_indices.size() - 1);
    while (!odometer.revs()) {
      std::vector<RDKit::Bond::BondType> neighbor_bond_types (n_neighbors);
      for (std::size_t i = 0; i < n_neighbors; ++i) {
        neighbor_bond_types[i] = (*candidate_bond_types)[odometer[i]];
      };
      auto perturbation = std::make_shared<AtomInsertion>(
        molecule,
        neighbor_atom_indices,
        neighbor_bond_types,
        (*candidate_atomic_numbers)[odometer[n_neighbors]],
        (*candidate_formal_charges)[odometer[n_neighbors + 1]],
        (*candidate_n_explicit_hydrogens)[odometer[n_neighbors + 2]],
        candidate_dropped_atom_indices[odometer[n_neighbors + 3]]);
      ++odometer;
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void AtomInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const boost::dynamic_bitset<>& candidate_neighbor_atoms_mask,
    const MolecularConstraints* constraints = nullptr) const {
    assert(molecule.getNumAtoms() == candidate_neighbor_atoms_mask.size());
    Combinator combinator (candidate_neighbor_atoms_mask.count(),
      allow_disconnections || !molecule.getNumAtoms() ? 0 : 1,
      atom_insertion_max_n_neighbors);
    while (!combinator.revs) {
      AtomInsertions(queue, molecule,
        ExtractBits(candidate_neighbor_atoms_mask, combinator.mask), constraints);
      ++combinator;
    };
  };

  void AtomInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx central_atom_idx,
    const MolecularConstraints* constraints = nullptr) const {
    boost::dynamic_bitset<> candidate_neighbor_atoms_mask =
      AtomInsertionNeighborMask(molecule, central_atom_idx);
    Combinator combinator (candidate_neighbor_atoms_mask.count(),
      0, atom_insertion_max_n_neighbors);
    while (!combinator.revs) {
      std::vector<std::size_t> neighbor_atom_indices = ExtractBits(
        candidate_neighbor_atoms_mask, combinator.mask);
      neighbor_atom_indices.push_back(central_atom_idx);
      AtomInsertions(
        queue, molecule, neighbor_atom_indices, constraints);
      ++combinator;
    };
  };

  void AtomInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    if (!n_atoms) {
      AtomInsertions(queue, molecule,
        boost::dynamic_bitset<>(n_atoms), constraints);
      return;
    };
    for (AtomIdx atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      AtomInsertions(queue, molecule, atom_idx, constraints);
    };
  };


  std::shared_ptr<AtomDeletion> DeleteAtom(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    const boost::dynamic_bitset<>& candidate_reconnection_atoms_mask,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    if (CanDeleteAtomWithoutReconnection(molecule, atom_idx)) {
      auto perturbation = std::make_shared<AtomDeletion>(atom_idx);
      if (!constraints || constraints->IsAllowed(*perturbation, molecule)) {
        return perturbation;
      };
    };
    if (!atom_deletion_allow_reconnections) {
      return nullptr;
    };
    std::vector<std::size_t> reconnection_atom_indices =
      ExtractBits(candidate_reconnection_atoms_mask);
    std::shuffle(
      reconnection_atom_indices.begin(), reconnection_atom_indices.end(), prng);
    for (std::size_t reconnection_atom_idx : reconnection_atom_indices) {
      auto perturbation = std::make_shared<AtomDeletion>(
        molecule, atom_idx, reconnection_atom_idx,
        atom_deletion_preserve_bond_types_during_reconnection,
        default_bond_type);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<AtomDeletion> DeleteAtom(
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    return DeleteAtom(molecule, atom_idx,
      AtomDeletionReconnectionMask(molecule, atom_idx), prng, constraints);
  };

  std::shared_ptr<AtomDeletion> DeleteAtom(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    for (std::size_t atom_idx : ShuffledAtomIndices(molecule, prng)) {
      auto perturbation = DeleteAtom(molecule, atom_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void AtomDeletions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    const boost::dynamic_bitset<>& candidate_reconnection_atoms_mask,
    const MolecularConstraints* constraints = nullptr) const {
    if (CanDeleteAtomWithoutReconnection(molecule, atom_idx)) {
      auto perturbation = std::make_shared<AtomDeletion>(atom_idx);
      if (!constraints || constraints->IsAllowed(*perturbation, molecule)) {
        queue.push(perturbation);
      };
    };
    if (!atom_deletion_allow_reconnections) {
      return;
    };
    for (std::size_t raix = candidate_reconnection_atoms_mask.find_first();
      raix != boost::dynamic_bitset<>::npos;
      raix = candidate_reconnection_atoms_mask.find_next(raix)) {
      auto perturbation = std::make_shared<AtomDeletion>(
        molecule, atom_idx, raix,
        atom_deletion_preserve_bond_types_during_reconnection,
        default_bond_type);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void AtomDeletions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx atom_idx,
    const MolecularConstraints* constraints = nullptr) const {
    AtomDeletions(queue, molecule, atom_idx,
      AtomDeletionReconnectionMask(molecule, atom_idx), constraints);
  };

  void AtomDeletions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      AtomDeletions(queue, molecule, atom_idx, constraints);
    };
  };


  std::shared_ptr<BondInsertion> InsertBond(
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx,
    const boost::dynamic_bitset<>& candidate_end_atom_mask,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    std::vector<std::size_t> end_atom_indices =
      ExtractBits(candidate_end_atom_mask);
    std::shuffle(end_atom_indices.begin(), end_atom_indices.end(), prng);
    for (std::size_t end_atom_idx : end_atom_indices) {
      std::vector<RDKit::Bond::BondType> candidate_bond_types {default_bond_type};
      if (bond_insertion_randomize_bond_type) {
        candidate_bond_types = WeightedShuffle(
          bond_types, bond_types_weights, prng);
      };
      for (RDKit::Bond::BondType bond_type : candidate_bond_types) {
        auto perturbation = std::make_shared<BondInsertion>(
          molecule, begin_atom_idx, end_atom_idx, bond_type);
        if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
          continue;
        };
        return perturbation;
      };
    };
    return nullptr;
  };

  std::shared_ptr<BondInsertion> InsertBond(
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    return InsertBond(molecule, begin_atom_idx,
      BondingPartnerMask(molecule, begin_atom_idx), prng, constraints);
  };

  std::shared_ptr<BondInsertion> InsertBond(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    boost::dynamic_bitset<> candidate_atom_mask (molecule.getNumAtoms());
    candidate_atom_mask.set();
    if (bond_insertion_max_atom_n_rings_membership < max_unsigned) {
      candidate_atom_mask = RingMembershipMask(
        molecule, 0, bond_insertion_max_atom_n_rings_membership - 1);
    };
    for (std::size_t begin_atom_idx : ShuffledAtomIndices(molecule, prng)) {
      if (!candidate_atom_mask[begin_atom_idx]) {
        continue;
      };
      auto perturbation = InsertBond(
        molecule, begin_atom_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void BondInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx,
    const boost::dynamic_bitset<>& candidate_end_atom_mask,
    const MolecularConstraints* constraints = nullptr) const {
    std::shared_ptr<const std::vector<RDKit::Bond::BondType>> candidate_bond_types;
    if (bond_insertion_iterate_bond_types) {
      candidate_bond_types.reset(&bond_types, boost::null_deleter());
    } else {
      candidate_bond_types.reset(
        new std::vector<RDKit::Bond::BondType> {default_bond_type});
    };
    for (std::size_t end_atom_idx = candidate_end_atom_mask.find_first();
      end_atom_idx != boost::dynamic_bitset<>::npos;
      end_atom_idx = candidate_end_atom_mask.find_next(end_atom_idx)) {
      for (RDKit::Bond::BondType bond_type : *candidate_bond_types) {
        auto perturbation = std::make_shared<BondInsertion>(
          molecule, begin_atom_idx, end_atom_idx, bond_type);
        if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
          continue;
        };
        queue.push(perturbation);
      };
    };
  };

  void BondInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    AtomIdx begin_atom_idx,
    const MolecularConstraints* constraints = nullptr) const {
    BondInsertions(queue, molecule, begin_atom_idx,
      BondingPartnerMask(molecule, begin_atom_idx), constraints);
  };

  void BondInsertions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr) const {
    std::size_t n_atoms = molecule.getNumAtoms();
    boost::dynamic_bitset<> candidate_atom_mask (n_atoms);
    candidate_atom_mask.set();
    if (bond_insertion_max_atom_n_rings_membership < max_unsigned) {
      candidate_atom_mask = RingMembershipMask(
        molecule, 0, bond_insertion_max_atom_n_rings_membership - 1);
    };
    for (std::size_t atom_idx = candidate_atom_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = candidate_atom_mask.find_next(atom_idx)) {
      BondInsertions(queue, molecule, atom_idx, constraints);
    };
  };


  std::shared_ptr<BondDeletion> DeleteBond(
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    const std::vector<std::pair<std::size_t, std::size_t>>& candidate_reroutes,
    const MolecularConstraints* constraints = nullptr) const {
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    AtomIdx begin_atom_idx = bond->getBeginAtomIdx();
    AtomIdx end_atom_idx = bond->getEndAtomIdx();
    if (CanDeleteBondWithoutReconnection(molecule, bond_idx)) {
      auto perturbation = std::make_shared<BondDeletion>(
        begin_atom_idx, end_atom_idx);
      if (!constraints || constraints->IsAllowed(*perturbation, molecule)) {
        return perturbation;
      };
    };
    if (!bond_deletion_allow_reroutes) {
      return nullptr;
    };
    for (const auto [reroute_begin_idx, reroute_end_idx] : candidate_reroutes) {
      auto perturbation = std::make_shared<BondDeletion>(
          molecule, begin_atom_idx, end_atom_idx,
          reroute_begin_idx, reroute_end_idx,
          bond_deletion_preserve_bond_types_during_reroute, default_bond_type);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      return perturbation;
    };
    return nullptr;
  };

  std::shared_ptr<BondDeletion> DeleteBond(
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    auto bond_reroutes = BondReroutes(molecule, bond_idx,
      bond_deletion_min_distance_reroute,
      bond_deletion_max_distance_reroute);
    std::shuffle(bond_reroutes.begin(), bond_reroutes.end(), prng);
    return DeleteBond(molecule, bond_idx, bond_reroutes, constraints);
  };

  std::shared_ptr<BondDeletion> DeleteBond(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    for (std::size_t bond_idx : ShuffledBondIndices(molecule, prng)) {
      auto perturbation = DeleteBond(molecule, bond_idx, prng, constraints);
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void BondDeletions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    const std::vector<std::pair<std::size_t, std::size_t>>& candidate_reroutes,
    const MolecularConstraints* constraints = nullptr) const {
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    AtomIdx begin_atom_idx = bond->getBeginAtomIdx();
    AtomIdx end_atom_idx = bond->getEndAtomIdx();
    if (CanDeleteBondWithoutReconnection(molecule, bond_idx)) {
      auto perturbation = std::make_shared<BondDeletion>(
        begin_atom_idx, end_atom_idx);
      if (!constraints || constraints->IsAllowed(*perturbation, molecule)) {
        queue.push(perturbation);
      };
    };
    if (!bond_deletion_allow_reroutes) {
      return;
    };
    for (const auto [reroute_begin_idx, reroute_end_idx] : candidate_reroutes) {
      auto perturbation = std::make_shared<BondDeletion>(
        molecule, begin_atom_idx, end_atom_idx,
        reroute_begin_idx, reroute_end_idx,
        bond_deletion_preserve_bond_types_during_reroute, default_bond_type);
      if (constraints && !constraints->IsAllowed(*perturbation, molecule)) {
        continue;
      };
      queue.push(perturbation);
    };
  };

  void BondDeletions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    BondIdx bond_idx,
    const MolecularConstraints* constraints = nullptr) const {
    BondDeletions(queue, molecule, bond_idx,
      BondReroutes(molecule, bond_idx,
        bond_deletion_min_distance_reroute,
        bond_deletion_min_distance_reroute), constraints);
  };

  void BondDeletions(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr) const {
    std::size_t n_bonds = molecule.getNumBonds();
    for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
      BondDeletions(queue, molecule, bond_idx, constraints);
    };
  };


  std::shared_ptr<MolecularPerturbation> operator()(
    const RDKit::ROMol& molecule,
    std::mt19937& prng,
    const MolecularConstraints* constraints = nullptr) const {
    std::vector<MolecularPerturbation::Type> eligible_perturbation_types =
      ShuffledEligiblePerturbationTypes(molecule, prng);
    if (eligible_perturbation_types.empty()) {
      return nullptr;
    };
    std::shared_ptr<MolecularPerturbation> perturbation;
    for (MolecularPerturbation::Type type : eligible_perturbation_types) {
      // std::cout << "Creating perturbation: " << type << std::endl;
      switch (type) {
        case MolecularPerturbation::AtomicNumberChange_t:
          perturbation = ChangeAtomicNumber(molecule, prng, constraints);
          break;
        case MolecularPerturbation::FormalChargeChange_t:
          perturbation = ChangeFormalCharge(molecule, prng, constraints);
          break;
        case MolecularPerturbation::ExplicitHydrogensChange_t:
          perturbation = ChangeExplicitHydrogens(molecule, prng, constraints);
          break;
        case MolecularPerturbation::BondTypeChange_t:
          perturbation = ChangeBondType(molecule, prng, constraints);
          break;
        case MolecularPerturbation::AtomInsertion_t:
          perturbation = InsertAtom(molecule, prng, constraints);
          break;
        case MolecularPerturbation::AtomDeletion_t:
          perturbation = DeleteAtom(molecule, prng, constraints);
          break;
        case MolecularPerturbation::BondInsertion_t:
          perturbation = InsertBond(molecule, prng, constraints);
          break;
        case MolecularPerturbation::BondDeletion_t:
          perturbation = DeleteBond(molecule, prng, constraints);
          break;
      };
      if (perturbation) {
        return perturbation;
      };
    };
    return nullptr;
  };

  void operator()(
    MolecularPerturbationQueue& queue,
    const RDKit::ROMol& molecule,
    const MolecularConstraints* constraints = nullptr) const {
    if (perturbation_types[MolecularPerturbation::AtomicNumberChange_t]) {
      AtomicNumberChanges(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::FormalChargeChange_t]) {
      FormalChargeChanges(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::ExplicitHydrogensChange_t]) {
      ExplicitHydrogenChanges(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::BondTypeChange_t]) {
      BondTypeChanges(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::AtomInsertion_t]) {
      AtomInsertions(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::AtomDeletion_t]) {
      AtomDeletions(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::BondInsertion_t]) {
      BondInsertions(queue, molecule, constraints);
    };
    if (perturbation_types[MolecularPerturbation::BondDeletion_t]) {
      BondDeletions(queue, molecule, constraints);
    };
  };
};

#endif // !_MOLECULE_PERTURBER_HPP_
