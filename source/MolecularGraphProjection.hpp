#pragma once
#ifndef _MOLECULAR_GRAPH_PROJECTION_HPP_
#define _MOLECULAR_GRAPH_PROJECTION_HPP_

#include "MolecularKeys.hpp"
#include "MolecularTags.hpp"
#include "CircularAtomicEnvironment.hpp"
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

// Some operations may reference atoms/bonds that aren't part of the molecular
// graph, either because they were removed or because they weren't added yet,
// despite being projected to be a part in the future. This makes iterating over
// bonds associated with atoms and viceversa difficult. This class specifies
// the expected topology. Despite it involving a copy of the topology it is much
// more lightweight than a full copy of a RDKit::RWMol.
class MolecularGraphProjection {
public:

  struct Atom {
  private:
    Tag atom_tag;
    bool is_null = true;

  public:
    std::uint8_t atomic_number;
    std::int8_t formal_charge;
    std::uint8_t n_explicit_hydrogens;

    Atom() = default;
    Atom(
      Tag atom_tag,
      std::uint8_t atomic_number,
      std::int8_t formal_charge,
      std::uint8_t n_explicit_hydrogens) :
      atom_tag(atom_tag),
      atomic_number(atomic_number),
      formal_charge(formal_charge),
      n_explicit_hydrogens(n_explicit_hydrogens),
      is_null(false) {};
    Atom(const RDKit::Atom* atom) {
      Set(atom);
    };

    void Set(const RDKit::Atom* atom) {
      atom_tag = GetTag(atom);
      atomic_number = atom->getAtomicNum();
      formal_charge = atom->getFormalCharge();
      n_explicit_hydrogens = atom->getNumExplicitHs();
      is_null = false;
    };

    void Null() {
      is_null = true;
    };

    Tag GetAtomTag() const {
      return atom_tag;
    };

    bool IsNull() const {
      return is_null;
    };
  };


  struct Bond {
  private:
    Tag bond_tag;
    bool is_null = true;

  public:
    RDKit::Bond::BondType bond_type;

    Bond() = default;
    Bond(Tag bond_tag, RDKit::Bond::BondType bond_type) :
      bond_tag(bond_tag), bond_type(bond_type), is_null(false) {};
    Bond(const RDKit::Bond* bond) {
      Set(bond);
    };

    void Set(const RDKit::Bond* bond) {
      bond_tag = GetTag(bond);
      bond_type = bond->getBondType();
      is_null = false;
    };

    void Null() {
      is_null = true;
    };

    Tag GetBondTag() const {
      return bond_tag;
    };

    bool IsNull() const {
      return is_null;
    };
  };

private:
  typedef boost::property<boost::edge_weight_t, std::size_t, Bond> EdgeProperty;
  typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS, Atom, EdgeProperty> Topology;

  typedef boost::exterior_vertex_property<
    MolecularGraphProjection::Topology, std::size_t> DistanceProperty;
  typedef DistanceProperty::matrix_type DistanceMatrix;

  static constexpr std::size_t inf = std::numeric_limits<std::size_t>::max();

  Topology topology;
  boost::dynamic_bitset<> affected_atoms;

private:
  AtomKey GetAtomKey(std::size_t atom_idx) const {
    const Atom& atom = topology[atom_idx];
    if (atom.IsNull()) {
      return AtomKey();
    };
    double valence = 0.0;
    for (std::size_t neighbor_idx : boost::make_iterator_range(
      boost::adjacent_vertices(atom_idx, topology))) {
      auto [edge, exists] = boost::edge(atom_idx, neighbor_idx, topology);
      const Bond& bond = topology[edge];
      valence += BondTypeValenceContribution(bond.bond_type);
    };
    return AtomKey(
      boost::out_degree(atom_idx, topology),
      valence,
      atom.atomic_number,
      atom.formal_charge,
      atom.n_explicit_hydrogens);
  };

  DistanceMatrix AllPairsShortestPaths() const {
    DistanceMatrix distance_matrix (boost::num_vertices(topology));
    boost::johnson_all_pairs_shortest_paths(topology, distance_matrix,
      boost::distance_inf(inf));
    return distance_matrix;
  };

  boost::dynamic_bitset<> ShortestPath(
    std::size_t i, std::size_t j,
    const DistanceMatrix& distance_matrix) const {
    boost::dynamic_bitset<> path (boost::num_vertices(topology));
    if (distance_matrix[i][j] >= inf) {
      return path;
    };
    path.set(i);
    path.set(j);
    while (i != j) {
      for (std::size_t adj : boost::make_iterator_range(
        boost::adjacent_vertices(i, topology))) {
        if (distance_matrix[adj][j] < distance_matrix[i][j]) {
          path.set(adj);
          i = adj;
          break;
        };
      };
    };
    return path;
  };

  std::vector<boost::dynamic_bitset<>> HortonCycles(
    const DistanceMatrix& distance_matrix) const {
    std::size_t n = boost::num_vertices(topology);
    std::vector<boost::dynamic_bitset<>> cycles;
    for (Topology::edge_descriptor e :
      boost::make_iterator_range(boost::edges(topology))) {
      Topology::vertex_descriptor s = boost::source(e, topology);
      Topology::vertex_descriptor t = boost::target(e, topology);
      for (Topology::vertex_descriptor v :
        boost::make_iterator_range(boost::vertices(topology))) {
        if (v == s || v == t) {
          continue;
        };
        boost::dynamic_bitset<> path_to_s = ShortestPath(v, s, distance_matrix);
        boost::dynamic_bitset<> path_to_t = ShortestPath(v, t, distance_matrix);
        boost::dynamic_bitset<> cycle  = path_to_s | path_to_t;
        if (cycle.count() + 1 == path_to_s.count() + path_to_t.count()) {
          cycles.push_back(std::move(cycle));
        };
      };
    };
    return cycles;
  };

  boost::dynamic_bitset<> CircularEnvironment(
    Topology::vertex_descriptor vertex,
    std::uint8_t radius) const {
    std::queue<Topology::vertex_descriptor> queue;
    boost::dynamic_bitset<> queued (boost::num_vertices(topology));
    boost::dynamic_bitset<> environment (boost::num_vertices(topology));
    queue.push(vertex);
    queued.set(vertex);
    environment.set(vertex);
    for (std::uint8_t r = 0; r <= radius; ++r) {
      std::size_t n = queue.size();
      if (!n) {
        break;
      };
      for (std::size_t i = 0; i < n; ++i) {
        Topology::vertex_descriptor v = queue.front();
        environment.set(v);
        for (Topology::vertex_descriptor nbr : boost::make_iterator_range(
          boost::adjacent_vertices(v, topology))) {
          if (!queued[nbr]) {
            queue.push(nbr);
            queued.set(nbr);
          };
        };
        queue.pop();
      };
    };
    return environment;
  };

  std::vector<std::uint64_t> AtomHashes() const {
    std::vector<std::uint64_t> atom_hashes;
    atom_hashes.reserve(boost::num_vertices(topology));
    for (auto vertex : boost::make_iterator_range(boost::vertices(topology))) {
      AtomKey atom_key = GetAtomKey(vertex);
      atom_hashes.push_back(hash_value(atom_key));
    };
    return atom_hashes;
  };

  static std::uint64_t BondTypeAsHash(const Bond& bond) {
    return bond.bond_type;
  };

  EnvironmentKey CircularEnvironmentKey(
    std::size_t atom_idx,
    const std::vector<std::uint64_t>& atom_hashes,
    const boost::dynamic_bitset<>& environment_mask,
    std::uint8_t radius) const {
    if (topology[atom_idx].IsNull()) {
      return 0;
    };
    return MorganHash<Topology, Bond>(
      atom_hashes,
      environment_mask,
      topology,
      (radius + 1) / 2,
      BondTypeAsHash);
  };

public:
  MolecularGraphProjection(const RDKit::ROMol& molecule) {
    std::size_t n = molecule.getNumAtoms();
    topology = Topology(n);
    for (const RDKit::Atom* atom : molecule.atoms()) {
      topology[atom->getIdx()].Set(atom);
    };
    for (const RDKit::Bond* bond : molecule.bonds()) {
      boost::add_edge(bond->getBeginAtomIdx(), bond->getEndAtomIdx(),
        {1, Bond(bond)}, topology);
    };
    affected_atoms.resize(n);
  };

  std::size_t AddAtom(
    Tag atom_tag,
    std::uint8_t atomic_number,
    std::int8_t formal_charge,
    std::uint8_t n_explicit_hydrogens) {
    std::size_t atom_idx = boost::add_vertex(
      Atom(atom_tag, atomic_number, formal_charge, n_explicit_hydrogens),
        topology);
    affected_atoms.resize(atom_idx + 1, true);
    return atom_idx;
  };

  std::size_t AddBond(
    Tag bond_tag,
    std::size_t begin_atom_idx,
    std::size_t end_atom_idx,
    RDKit::Bond::BondType bond_type) {
    boost::add_edge(begin_atom_idx, end_atom_idx,
      {1, Bond(bond_tag, bond_type)}, topology);
    affected_atoms.set(begin_atom_idx);
    affected_atoms.set(end_atom_idx);
    return boost::num_edges(topology) - 1;
  };

  void RemoveAtom(std::size_t atom_idx) {
    // Be careful! This won't actually remove the atom from the graph but
    // rather null it and remove all of its associated bonds. This is done on
    // purpose to avoid vertex descriptor invalidation.
    Atom& atom = topology[atom_idx];
    atom.Null();
    for (auto neighbor_idx : boost::make_iterator_range(
      boost::adjacent_vertices(atom_idx, topology))) {
      affected_atoms.set(neighbor_idx);
    };
    boost::clear_vertex(atom_idx, topology);
    affected_atoms.set(atom_idx);
  };

  void RemoveBond(std::size_t begin_atom_idx, std::size_t end_atom_idx) {
    boost::remove_edge(begin_atom_idx, end_atom_idx, topology);
    affected_atoms.set(begin_atom_idx);
    affected_atoms.set(end_atom_idx);
  };

  Atom& EditAtom(std::size_t atom_idx) {
    affected_atoms.set(atom_idx);
    return topology[atom_idx];
  };

  Bond& EditBond(std::size_t begin_atom_idx, std::size_t end_atom_idx) {
    auto [edge, exists] = boost::edge(begin_atom_idx, end_atom_idx, topology);
    if (!exists) {
      throw std::out_of_range("No Bond exists between input Atoms.");
    };
    affected_atoms.set(begin_atom_idx);
    affected_atoms.set(end_atom_idx);
    return topology[edge];
  };

  std::size_t Size() const {
    return boost::num_vertices(topology);
  };

  std::vector<boost::dynamic_bitset<>> MinimumCycleBasis() const {
    // Finds the Minimum Cycle Basis of the graph using the Horton algorithm.
    // Horton, J.D. A polynomial-time algorithm to find the shortest cycle
    // basis of a graph. SIAM Journal on Computing 16(2), 358-366 (1987)
    // Note that this isn't the algorithm with the best complexity to solve the
    // problem, but it's easy to implement and understand, hence its use here.
    boost::dynamic_bitset<> assigned_vertices (boost::num_vertices(topology));
    std::vector<boost::dynamic_bitset<>> minimum_cycle_basis;
    std::vector<boost::dynamic_bitset<>> cycles =
      HortonCycles(AllPairsShortestPaths());
    std::sort(cycles.begin(), cycles.end(),
      [](const auto& c1, const auto& c2) {
        return c1.count() < c2.count();
      });
    for (const auto& cycle : cycles) {
      if (!cycle.is_subset_of(assigned_vertices)) {
        assigned_vertices |= cycle;
        minimum_cycle_basis.push_back(cycle);
      };
    };
    return minimum_cycle_basis;
  };

  std::pair<std::map<Tag, AtomKeyChange>, std::map<Tag, BondKeyChange>>
  MolecularKeyChanges(
    const RDKit::ROMol& molecule,
    bool calculate_bond_key_changes = true) {
    std::map<Tag, AtomKeyChange> atom_key_changes;
    std::map<Tag, BondKeyChange> bond_key_changes;
    if (affected_atoms.none()) {
      return {std::move(atom_key_changes), std::move(bond_key_changes)};
    };
    std::size_t n_atoms = molecule.getNumAtoms();
    for (std::size_t atom_idx = affected_atoms.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = affected_atoms.find_next(atom_idx)) {
      // Calculate AtomKey change.
      AtomKey prior_atom_key;
      if (atom_idx < n_atoms) {
        prior_atom_key = AtomKey(molecule.getAtomWithIdx(atom_idx));
      };
      AtomKey posterior_atom_key = GetAtomKey(atom_idx);
      if (prior_atom_key != posterior_atom_key) {
        atom_key_changes.try_emplace(
          topology[atom_idx].GetAtomTag(), prior_atom_key, posterior_atom_key);
      };
      if (!calculate_bond_key_changes) {
        continue;
      };
      // Calculate BondKey changes for bonds that exist in the original molecule.
      if (atom_idx < n_atoms) {
        const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
        for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
          Tag bond_tag = GetTag(bond);
          if (bond_key_changes.contains(bond_tag)) {
            continue;
          };
          const RDKit::Atom* neighbor = bond->getOtherAtom(atom);
          std::size_t neighbor_idx = neighbor->getIdx();
          BondKey prior_bond_key (
            prior_atom_key, AtomKey(neighbor), bond->getBondType());
          BondKey posterior_bond_key;
          auto [edge, exists] = boost::edge(atom_idx, neighbor_idx, topology);
          if (exists) {
            posterior_bond_key = BondKey(
              posterior_atom_key,
              GetAtomKey(neighbor_idx),
              topology[edge].bond_type);
          };
          if (prior_bond_key == posterior_bond_key) {
            continue;
          };
          bond_key_changes.try_emplace(
            bond_tag, prior_bond_key, posterior_bond_key);
        };
      };
      // Calculate BondKey changes for bonds that exist in the projection.
      for (std::size_t neighbor_idx : boost::make_iterator_range(
        boost::adjacent_vertices(atom_idx, topology))) {
        auto [edge, exists] = boost::edge(atom_idx, neighbor_idx, topology);
        const Bond& bond = topology[edge];
        if (bond_key_changes.contains(bond.GetBondTag())) {
          continue;
        };
        BondKey prior_bond_key;
        const RDKit::Bond* b = nullptr;
        if (atom_idx < n_atoms && neighbor_idx < n_atoms) {
          b = molecule.getBondBetweenAtoms(atom_idx, neighbor_idx);
        };
        if (b) {
          prior_bond_key = BondKey(
            prior_atom_key,
            AtomKey(molecule.getAtomWithIdx(neighbor_idx)),
            b->getBondType());
        };
        BondKey posterior_bond_key (
          posterior_atom_key, GetAtomKey(neighbor_idx), bond.bond_type);
        if (prior_bond_key == posterior_bond_key) {
          continue;
        };
        bond_key_changes.try_emplace(
          bond.GetBondTag(), prior_bond_key, posterior_bond_key);
      };
    };
    return {std::move(atom_key_changes), std::move(bond_key_changes)};
  };

  std::map<Tag, EnvironmentKeyChange> EnvironmentKeyChanges(
    const RDKit::ROMol& molecule,
    std::uint8_t radius) const {
    std::map<Tag, EnvironmentKeyChange> environment_key_changes;
    if (affected_atoms.none()) {
      return environment_key_changes;
    };
    std::vector<boost::dynamic_bitset<>> posterior_environments (
      boost::num_vertices(topology));
    boost::dynamic_bitset<> affected_environments (affected_atoms);
    for (std::size_t atom_idx = affected_atoms.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = affected_atoms.find_next(atom_idx)) {
      boost::dynamic_bitset<> environment = CircularEnvironment(atom_idx, radius);
      affected_environments |= environment;
      posterior_environments[atom_idx] = std::move(environment);
    };
    std::size_t n_atoms = molecule.getNumAtoms();
    std::vector<std::uint64_t> prior_atom_hashes = AtomKeyHashes(molecule);
    std::vector<std::uint64_t> posterior_atom_hashes = AtomHashes();
    for (std::size_t atom_idx = affected_environments.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = affected_environments.find_next(atom_idx)) {
      EnvironmentKey prior_environment_key = 0;
      if (atom_idx < n_atoms) {
        CircularAtomicEnvironment prior_environment (
          molecule.getAtomWithIdx(atom_idx), radius);
        prior_environment_key = prior_environment.Key(prior_atom_hashes);
      };
      boost::dynamic_bitset<> posterior_environment =
        posterior_environments[atom_idx];
      if (!posterior_environment.size()) {
        posterior_environment = CircularEnvironment(atom_idx, radius);
      };
      EnvironmentKey posterior_environment_key = CircularEnvironmentKey(
        atom_idx, posterior_atom_hashes, posterior_environment, radius);
      if (prior_environment_key != posterior_environment_key) {
        environment_key_changes.try_emplace(
          topology[atom_idx].GetAtomTag(),
          prior_environment_key,
          posterior_environment_key);
      };
    };
    return environment_key_changes;
  };
};

#endif // !_MOLECULAR_GRAPH_PROJECTION_HPP_
