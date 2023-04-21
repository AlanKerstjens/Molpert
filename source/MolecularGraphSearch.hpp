#pragma once
#ifndef _MOLECULAR_GRAPH_SEARCH_HPP_
#define _MOLECULAR_GRAPH_SEARCH_HPP_

#include "RandomSampling.hpp"
#include <GraphMol/MolOps.h>
#include <queue>

// First = Atom/Bond index, Second = Discovery time (e.g. search depth)
typedef std::pair<std::size_t, std::size_t> MolecularGraphDiscovery;

class MolecularGraphDiscoveryJournal {
  std::vector<MolecularGraphDiscovery> atom_discoveries, bond_discoveries;
  boost::dynamic_bitset<> discovered_atoms_mask, discovered_bonds_mask;
  bool record_discovery_times = true;

public:
  MolecularGraphDiscoveryJournal(
    std::size_t n_atoms,
    std::size_t n_bonds,
    bool record_discovery_times = true) :
    record_discovery_times(record_discovery_times) {
    discovered_atoms_mask.resize(n_atoms);
    discovered_bonds_mask.resize(n_bonds);
  };

  void DiscoverAtom(std::size_t atom_idx, std::size_t time) {
    discovered_atoms_mask.set(atom_idx);
    if (record_discovery_times) {
      atom_discoveries.emplace_back(atom_idx, time);
    };
  };

  void DiscoverAtoms(
    const std::vector<std::size_t>& atom_indices, std::size_t time) {
    for (std::size_t atom_idx : atom_indices) {
      DiscoverAtom(atom_idx, time);
    };
  };

  void DiscoverAtoms(
    const boost::dynamic_bitset<>& atoms_mask, std::size_t time) {
    for (std::size_t atom_idx = atoms_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = atoms_mask.find_next(atom_idx)) {
      DiscoverAtom(atom_idx, time);
    };
  };

  void DiscoverBond(
    std::size_t bond_idx,
    std::size_t time) {
    discovered_bonds_mask.set(bond_idx);
    if (record_discovery_times) {
      bond_discoveries.emplace_back(bond_idx, time);
    };
  };

  void DiscoverBonds(
    const std::vector<std::size_t>& bond_indices, std::size_t time) {
    for (std::size_t bond_idx : bond_indices) {
      DiscoverBond(bond_idx, time);
    };
  };

  void DiscoverBonds(
    const boost::dynamic_bitset<>& bonds_mask, std::size_t time) {
    for (std::size_t bond_idx = bonds_mask.find_first();
      bond_idx != boost::dynamic_bitset<>::npos;
      bond_idx = bonds_mask.find_next(bond_idx)) {
      DiscoverBond(bond_idx, time);
    };
  };

  std::pair<std::size_t, bool> GetAtomDiscoveryTime(std::size_t atom_idx) const {
    if (!discovered_atoms_mask[atom_idx]) {
      return {0, false};
    };
    for (auto [discovered_atom_idx, discovery_time] : atom_discoveries) {
      if (discovered_atom_idx == atom_idx) {
        return {discovery_time, true};
      };
    };
    return {0, false};
  };

  std::pair<std::size_t, bool> GetBondDiscoveryTime(std::size_t bond_idx) const {
    if (!discovered_bonds_mask[bond_idx]) {
      return {0, false};
    };
    for (auto [discovered_bond_idx, discovery_time] : bond_discoveries) {
      if (discovered_bond_idx == bond_idx) {
        return {discovery_time, true};
      };
    };
    return {0, false};
  };

  bool AtomHasBeenDiscovered(std::size_t atom_idx) const {
    return discovered_atoms_mask[atom_idx];
  };

  bool BondHasBeenDiscovered(std::size_t bond_idx) const {
    return discovered_bonds_mask[bond_idx];
  };

  std::size_t NDiscoveredAtoms() const {
    return discovered_atoms_mask.count();
  };

  std::size_t NDiscoveredBonds() const {
    return discovered_bonds_mask.count();
  };

  const MolecularGraphDiscovery& GetLastAtomDiscovery() const {
    return atom_discoveries.back();
  };

  const MolecularGraphDiscovery& GetLastBondDiscovery() const {
    return bond_discoveries.back();
  };

  const std::vector<MolecularGraphDiscovery>& GetAtomDiscoveries() const {
    return atom_discoveries;
  };

  const std::vector<MolecularGraphDiscovery>& GetBondDiscoveries() const {
    return bond_discoveries;
  };

  const boost::dynamic_bitset<>& GetDiscoveredAtomsMask() const {
    return discovered_atoms_mask;
  };

  const boost::dynamic_bitset<>& GetDiscoveredBondsMask() const {
    return discovered_bonds_mask;
  };
};

MolecularGraphDiscoveryJournal MoleculeBFS(
  const RDKit::ROMol& molecule,
  std::size_t start_atom_idx,
  std::size_t max_search_depth,
  const boost::dynamic_bitset<>& untraversable_atoms_mask,
  const boost::dynamic_bitset<>& untraversable_bonds_mask,
  const std::function<
    bool(const MolecularGraphDiscoveryJournal&)>& termination_criterion,
  bool record_discovery_times = true) {
  MolecularGraphDiscoveryJournal journal (
    molecule.getNumAtoms(), molecule.getNumBonds(), record_discovery_times);
  std::queue<MolecularGraphDiscovery> atom_queue;
  journal.DiscoverAtom(start_atom_idx, 0);
  if (termination_criterion(journal)) {
    return journal;
  };
  atom_queue.emplace(start_atom_idx, 0);
  while (!atom_queue.empty()) {
    const auto [atom_idx, search_depth] = atom_queue.front();
    if (search_depth >= max_search_depth) {
      atom_queue.pop();
      continue;
    };
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
      std::size_t bond_idx = bond->getIdx();
      if (!journal.BondHasBeenDiscovered(bond_idx) &&
          !untraversable_bonds_mask[bond_idx]) {
        journal.DiscoverBond(bond_idx, search_depth);
        if (termination_criterion(journal)) {
          return journal;
        };
        std::size_t neighbor_idx = bond->getOtherAtomIdx(atom_idx);
        if (!journal.AtomHasBeenDiscovered(neighbor_idx) &&
            !untraversable_atoms_mask[neighbor_idx]) {
          journal.DiscoverAtom(neighbor_idx, search_depth + 1);
          if (termination_criterion(journal)) {
            return journal;
          };
          atom_queue.emplace(neighbor_idx, search_depth + 1);
        };
      };
    };
    atom_queue.pop();
  };
  return journal;
};

std::vector<boost::dynamic_bitset<>> RingSystems(
  const RDKit::ROMol& molecule) {
  // Initialize the ring systems as the rings of the SSSR.
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  std::size_t n_ring_systems = ring_info->numRings();
  std::vector<boost::dynamic_bitset<>> ring_systems;
  ring_systems.reserve(n_ring_systems);
  for (const std::vector<int>& ring : ring_info->atomRings()) {
    ring_systems.emplace_back(molecule.getNumAtoms());
    boost::dynamic_bitset<>& ring_system = ring_systems.back();
    for (int atom_idx : ring) {
      ring_system.set(atom_idx);
    };
  };
  // Merge ring systems if they share an atom.
  for (std::size_t i = 0; i < n_ring_systems; ++i) {
    for (std::size_t j = i + 1; j < n_ring_systems; ++j) {
      if (ring_systems[i].intersects(ring_systems[j])) {
        ring_systems[i] |= ring_systems[j];
        std::swap(ring_systems[j], ring_systems[n_ring_systems - 1]);
        --j;
        --n_ring_systems;
      };
    };
  };
  ring_systems.erase(ring_systems.begin() + n_ring_systems, ring_systems.end());
  return ring_systems;
};

boost::dynamic_bitset<> BondsMask(
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask) {
  boost::dynamic_bitset<> bonds_mask (molecule.getNumBonds());
  for (std::size_t atom_idx = atoms_mask.find_first();
    atom_idx != boost::dynamic_bitset<>::npos;
    atom_idx = atoms_mask.find_next(atom_idx)) {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
      std::size_t neighbor_idx = bond->getOtherAtomIdx(atom_idx);
      if (atoms_mask[neighbor_idx]) {
        bonds_mask.set(bond->getIdx());
      };
    };
  };
  return bonds_mask;
};

MolecularGraphDiscoveryJournal MoleculeRandomWalk(
  const RDKit::ROMol& molecule,
  std::size_t start_atom_idx,
  std::size_t max_walk_length,
  std::mt19937& prng,
  const boost::dynamic_bitset<>& untraversable_atoms_mask,
  const boost::dynamic_bitset<>& untraversable_bonds_mask,
  bool walk_cycle_upon_discovery = false) {
  std::optional<std::vector<boost::dynamic_bitset<>>> ring_systems;
  if (walk_cycle_upon_discovery) {
    ring_systems = RingSystems(molecule);
  };
  MolecularGraphDiscoveryJournal journal (
    molecule.getNumAtoms(), molecule.getNumBonds(), false);
  boost::dynamic_bitset<> expandable_atoms_mask (molecule.getNumAtoms());
  journal.DiscoverAtom(start_atom_idx, 0);
  expandable_atoms_mask.set(start_atom_idx);
  while (journal.NDiscoveredAtoms() < max_walk_length &&
    expandable_atoms_mask.count()) {
    // Sample a random expandable atom.
    std::size_t atom_idx = Sample(expandable_atoms_mask, prng);
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    // Iterate over its bonds in a random order.
    std::vector<const RDKit::Bond*> bonds;
    bonds.reserve(atom->getDegree());
    for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
      bonds.push_back(bond);
    };
    std::shuffle(bonds.begin(), bonds.end(), prng);
    bool expanded_atom = false;
    for (const RDKit::Bond* bond : bonds) {
      std::size_t bond_idx = bond->getIdx();
      // If the bond or adjacent atom have been discovered or are untraversable
      // skip the bond.
      if (journal.BondHasBeenDiscovered(bond_idx) ||
        untraversable_bonds_mask[bond_idx]) {
        continue;
      };
      journal.DiscoverBond(bond_idx, 0);
      std::size_t neighbor_idx = bond->getOtherAtomIdx(atom_idx);
      if (journal.AtomHasBeenDiscovered(neighbor_idx) ||
        untraversable_atoms_mask[neighbor_idx]) {
        continue;
      };
      journal.DiscoverAtom(neighbor_idx, 0);
      expandable_atoms_mask.set(neighbor_idx);
      // If the adjacent atom is part of a cycle we may instantly discover the
      // whole ring system.
      if (walk_cycle_upon_discovery) {
        for (const boost::dynamic_bitset<>& ring_system : *ring_systems) {
          if (ring_system[neighbor_idx]) {
            journal.DiscoverAtoms(ring_system, 0);
            journal.DiscoverBonds(BondsMask(molecule, ring_system), 0);
            expandable_atoms_mask |= ring_system;
            break;
          };
        };
      };
      expanded_atom = true;
      break;
    };
    // If the sampled atom wasn't expanded it isn't expandable. Flag it as such.
    if (!expanded_atom) {
      expandable_atoms_mask.reset(atom_idx);
    };
  };
  return journal;
};

bool NeverMetCriterion(const MolecularGraphDiscoveryJournal&) {
  return false;
};

class DiscoveryCriterion {
  boost::dynamic_bitset<> atoms_to_discover_mask, bonds_to_discover_mask;

public:
  DiscoveryCriterion(
    const boost::dynamic_bitset<>& atoms_to_discover_mask,
    const boost::dynamic_bitset<>& bonds_to_discover_mask) :
    atoms_to_discover_mask(atoms_to_discover_mask),
    bonds_to_discover_mask(bonds_to_discover_mask) {};

  bool operator()(const MolecularGraphDiscoveryJournal& journal) const {
    if (!atoms_to_discover_mask.is_subset_of(journal.GetDiscoveredAtomsMask())) {
      return false;
    };
    return bonds_to_discover_mask.is_subset_of(journal.GetDiscoveredBondsMask());
  };
};

bool AtomIsInRing(const RDKit::ROMol& molecule, std::size_t atom_idx) {
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  return ring_info->numAtomRings(atom_idx);
};

bool BondIsInRing(const RDKit::ROMol& molecule, std::size_t bond_idx) {
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  return ring_info->numBondRings(bond_idx);
};

bool ConnectedAfterAtomDeletion(
  const RDKit::ROMol& molecule,
  std::size_t atom_idx) {
  // The graph can't be disconnected if the degree of the vertex is 0 or 1.
  const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
  if (atom->getDegree() <= 1) {
    return true;
  };
  // Otherwise, the graph remains connected after a deletion if all of the atom's
  // neighbors are reachable from one another without traversing the deleted atom.
  std::size_t n_atoms = molecule.getNumAtoms(), n_bonds = molecule.getNumBonds();
  boost::dynamic_bitset<> atoms_to_discover_mask (n_atoms);
  boost::dynamic_bitset<> untraversable_atoms_mask (n_atoms);
  std::size_t start_neighbor_idx;
  for (const RDKit::Atom* neighbor : molecule.atomNeighbors(atom)) {
    std::size_t neighbor_idx = neighbor->getIdx();
    start_neighbor_idx = neighbor_idx;
    atoms_to_discover_mask.set(neighbor_idx);
  };
  untraversable_atoms_mask.set(atom_idx);
  DiscoveryCriterion criterion (
    atoms_to_discover_mask,
    boost::dynamic_bitset<>(n_bonds));
  // Perform a BFS to check if this is the case or not.
  MolecularGraphDiscoveryJournal journal = MoleculeBFS(
    molecule,
    start_neighbor_idx,
    std::numeric_limits<std::size_t>::max(),
    untraversable_atoms_mask,
    boost::dynamic_bitset<>(n_bonds),
    criterion,
    false);
  return criterion(journal);
};

bool ConnectedAfterBondDeletion(
  const RDKit::ROMol& molecule,
  std::size_t bond_idx,
  bool calculate_sssr = false) {
  // If we have SSSR data (or we are willing to compute it) we can check if the
  // molecule remains connected by checking if the atom is part of a ring system.
  if (calculate_sssr) {
    // This call does nothing if RingInfo is already populated.
    RDKit::MolOps::findSSSR(molecule);
  };
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  if (ring_info->isInitialized()) {
    return ring_info->numBondRings(bond_idx) > 0;
  };
  std::size_t n_atoms = molecule.getNumAtoms(), n_bonds = molecule.getNumBonds();
  boost::dynamic_bitset<> atoms_to_discover_mask (n_atoms);
  boost::dynamic_bitset<> untraversable_bonds_mask (n_bonds);
  const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
  atoms_to_discover_mask.set(bond->getEndAtomIdx());
  untraversable_bonds_mask.set(bond_idx);
  DiscoveryCriterion criterion (
    atoms_to_discover_mask,
    boost::dynamic_bitset<>(n_bonds));
  MolecularGraphDiscoveryJournal journal = MoleculeBFS(
    molecule,
    bond->getBeginAtomIdx(),
    std::numeric_limits<std::size_t>::max(),
    boost::dynamic_bitset<>(n_atoms),
    untraversable_bonds_mask,
    criterion,
    false);
  return criterion(journal);
};

std::size_t DistanceBetweenAtoms(
  const RDKit::ROMol& molecule,
  std::size_t atom_idx1,
  std::size_t atom_idx2,
  bool use_topological_distance_matrix = false) {
  if (use_topological_distance_matrix) {
    // If the matrix was already calculated it is cached and we can fetch it.
    // Otherwise we compute it now. Note that the getDistanceMat() call doesn't
    // force the calculation of the matrix, meaning that if the cached matrix is
    // no longer valid due to molecule edition this function won't behave as
    // expected. The user is responsible for keeping the matrix updated.
    double* distance_matrix = RDKit::MolOps::getDistanceMat(molecule);
    return distance_matrix[atom_idx1 * molecule.getNumAtoms() + atom_idx2];
  };
  std::list<int> shortest_path =
    RDKit::MolOps::getShortestPath(molecule, atom_idx1, atom_idx2);
  if (shortest_path.empty()) {
    return std::numeric_limits<std::size_t>::max();
  };
  return shortest_path.size() - 1; // -1 because the input atoms are included
};

std::pair<std::size_t, bool> DistanceBetweenAtoms(
  const RDKit::ROMol& molecule,
  std::size_t atom_idx1,
  std::size_t atom_idx2,
  const boost::dynamic_bitset<>& untraversable_atoms_mask,
  const boost::dynamic_bitset<>& untraversable_bonds_mask) {
  if (atom_idx1 == atom_idx2) {
    return {0, true};
  };
  boost::dynamic_bitset<> atoms_to_discover_mask (molecule.getNumAtoms());
  atoms_to_discover_mask.set(atom_idx2);
  DiscoveryCriterion criterion (
    atoms_to_discover_mask,
    boost::dynamic_bitset<>(molecule.getNumBonds()));
  MolecularGraphDiscoveryJournal journal = MoleculeBFS(
    molecule,
    atom_idx1,
    std::numeric_limits<std::size_t>::max(),
    untraversable_atoms_mask,
    untraversable_bonds_mask,
    criterion,
    true);
  return journal.GetAtomDiscoveryTime(atom_idx2);
};

std::vector<std::size_t> AtomsWithinDistance(
  const RDKit::ROMol& molecule,
  std::size_t start_atom_idx,
  std::size_t min_distance,
  std::size_t max_distance,
  bool use_topological_distance_matrix = false) {
  std::size_t n_atoms = molecule.getNumAtoms();
  std::vector<std::size_t> atom_indices;
  // If the user requests it we can use the topological distance matrix to find
  // the relevant atoms.
  if (use_topological_distance_matrix) {
    // If the matrix was already calculated it is cached and we can fetch it.
    // Otherwise we compute it now. Note that the getDistanceMat() call doesn't
    // force the calculation of the matrix, meaning that if the cached matrix is
    // no longer valid due to molecule edition this function won't behave as
    // expected. The user is responsible for keeping the matrix updated.
    double* distance_matrix = RDKit::MolOps::getDistanceMat(molecule);
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      double distance = distance_matrix[start_atom_idx * n_atoms + atom_idx];
      if (distance >= min_distance && distance <= max_distance) {
        atom_indices.push_back(atom_idx);
      };
    };
    return atom_indices;
  };
  MolecularGraphDiscoveryJournal journal = MoleculeBFS(
    molecule,
    start_atom_idx,
    max_distance,
    boost::dynamic_bitset<>(n_atoms),
    boost::dynamic_bitset<>(molecule.getNumBonds()),
    NeverMetCriterion,
    true);
  const std::vector<MolecularGraphDiscovery>& atom_discoveries =
    journal.GetAtomDiscoveries();
  atom_indices.reserve(atom_discoveries.size());
  for (const auto& [atom_idx, distance] : atom_discoveries) {
    if (distance >= min_distance && distance <= max_distance) {
      atom_indices.push_back(atom_idx);
    };
  };
  return atom_indices;
};

std::vector<std::pair<std::size_t, std::size_t>> BondReroutes(
  const RDKit::ROMol& molecule,
  std::size_t bond_idx,
  std::size_t min_distance,
  std::size_t max_distance) {
  const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
  std::size_t begin_atom_idx = bond->getBeginAtomIdx();
  std::size_t end_atom_idx = bond->getEndAtomIdx();
  std::size_t bfs_depth = max_distance / 2;
  boost::dynamic_bitset<> untraversable_atoms_mask (molecule.getNumAtoms());
  untraversable_atoms_mask.set(end_atom_idx);
  MolecularGraphDiscoveryJournal begin_journal = MoleculeBFS(
    molecule,
    begin_atom_idx,
    bfs_depth,
    untraversable_atoms_mask,
    boost::dynamic_bitset<>(molecule.getNumBonds()),
    NeverMetCriterion,
    true);
  const std::vector<MolecularGraphDiscovery>& begin_discoveries =
    begin_journal.GetAtomDiscoveries();
  untraversable_atoms_mask.reset();
  untraversable_atoms_mask.set(begin_atom_idx);
  MolecularGraphDiscoveryJournal end_journal = MoleculeBFS(
    molecule,
    end_atom_idx,
    bfs_depth,
    untraversable_atoms_mask,
    boost::dynamic_bitset<>(molecule.getNumBonds()),
    NeverMetCriterion,
    true);
  const std::vector<MolecularGraphDiscovery>& end_discoveries =
    end_journal.GetAtomDiscoveries();
  std::vector<std::pair<std::size_t, std::size_t>> bond_reroutes;
  bond_reroutes.reserve(begin_discoveries.size() * end_discoveries.size());
  for (const auto& [begin_atom_idx, begin_distance] : begin_discoveries) {
    for (const auto& [end_atom_idx, end_distance] : end_discoveries) {
      // +1 because the distances are expressed starting from the begin and
      // end atoms respectively, which are separated by one bond.
      std::size_t distance = begin_distance + end_distance + 1;
      if (distance >= min_distance && distance <= max_distance) {
        bond_reroutes.emplace_back(begin_atom_idx, end_atom_idx);
      };
    };
  };
  return bond_reroutes;
};

#endif // !_MOLECULAR_GRAPH_SEARCH_HPP_
