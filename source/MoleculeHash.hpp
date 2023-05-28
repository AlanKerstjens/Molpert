#ifndef _MOLECULE_HASH_HPP_
#define _MOLECULE_HASH_HPP_

#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/container_hash/hash.hpp>
#include "MolecularKeys.hpp"

typedef std::function<
  std::vector<std::uint64_t>(const RDKit::ROMol&)> AtomsHasher;

std::uint64_t LoneAtomHash(const RDKit::Atom* atom) {
  std::uint64_t hash = atom->getAtomicNum();
  boost::hash_combine(hash, atom->getFormalCharge());
  // I used to include the total (implicit + explicit ) as opposed to the
  // explicit number of hydrogens, but in the RDKit the number of implicit
  // hydrogens isn't an atom member variable but rather inferred. This makes the
  // results rather inconsistent for aromatic environments, where different
  // kekule structures can lead to different distributions of implicit hydrogens.
  // Moreover, the number of implicit hydrogens depends on the number of non-H
  // neighboring atoms, infringing upon the "lone atom" concept.
  boost::hash_combine(hash, atom->getNumExplicitHs());
  return hash;
};

std::vector<std::uint64_t> LoneAtomHashes(const RDKit::ROMol& molecule) {
  std::size_t n = molecule.getNumAtoms();
  std::vector<std::uint64_t> atom_hashes (n);
  for (std::size_t atom_idx = 0; atom_idx < n; ++atom_idx) {
    atom_hashes[atom_idx] = LoneAtomHash(molecule.getAtomWithIdx(atom_idx));
  };
  return atom_hashes;
};

std::uint64_t RingAwareAtomHash(
  const RDKit::Atom* atom,
  bool is_in_ring) {
  std::uint64_t hash = LoneAtomHash(atom);
  boost::hash_combine(hash, is_in_ring);
  boost::hash_combine(hash, atom->getIsAromatic());
  return hash;
};

std::vector<std::uint64_t> RingAwareAtomHashes(const RDKit::ROMol& molecule) {
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  if (!ring_info->isInitialized()) {
    RDKit::MolOps::findSSSR(molecule);
  };
  std::size_t n = molecule.getNumAtoms();
  std::vector<std::uint64_t> atom_hashes (n);
  for (std::size_t atom_idx = 0; atom_idx < n; ++atom_idx) {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    atom_hashes[atom_idx] = RingAwareAtomHash(
      atom, ring_info->numAtomRings(atom_idx));
  };
  return atom_hashes;
};

std::uint64_t AtomKeyHash(const RDKit::Atom* atom) {
  AtomKey atom_key (atom);
  return hash_value(atom_key);
};

std::vector<std::uint64_t> AtomKeyHashes(const RDKit::ROMol& molecule) {
  std::size_t n = molecule.getNumAtoms();
  std::vector<std::uint64_t> atom_hashes (n);
  for (std::size_t atom_idx = 0; atom_idx < n; ++atom_idx) {
    atom_hashes[atom_idx] = AtomKeyHash(molecule.getAtomWithIdx(atom_idx));
  };
  return atom_hashes;
};


// The hash functions are templated to work with boost::graph::adjacency_lists.
// This allows us to re-use the code later on.
template <class Graph, class Edge>
void MorganHashUpdate(
  std::vector<std::uint64_t>& vertex_hashes,
  const boost::dynamic_bitset<>& vertex_mask,
  const Graph& graph,
  const std::function<
    std::uint64_t(const Edge&)>& edge_hasher = nullptr) {
  std::vector<std::uint64_t> updated_vertex_hashes (vertex_hashes.size());
  for (std::size_t vertex = vertex_mask.find_first();
    vertex != boost::dynamic_bitset<>::npos;
    vertex = vertex_mask.find_next(vertex)) {
    std::vector<std::uint64_t> neighbor_hashes;
    neighbor_hashes.reserve(boost::out_degree(vertex, graph));
    for (std::size_t neighbor : boost::make_iterator_range(
      boost::adjacent_vertices(vertex, graph))) {
      if (!vertex_mask[neighbor]) {
        continue;
      };
      std::uint64_t neighbor_hash = vertex_hashes[neighbor];
      if (edge_hasher) {
        auto [edge, _] = boost::edge(vertex, neighbor, graph);
        boost::hash_combine(neighbor_hash, edge_hasher(graph[edge]));
      };
      neighbor_hashes.push_back(neighbor_hash);
    };
    std::sort(neighbor_hashes.begin(), neighbor_hashes.end());
    std::uint64_t updated_hash = vertex_hashes[vertex];
    boost::hash_combine(updated_hash,
      boost::hash_range(neighbor_hashes.begin(), neighbor_hashes.end()));
    updated_vertex_hashes[vertex] = updated_hash;
  };
  vertex_hashes = updated_vertex_hashes;
};

template <class Graph, class Edge>
std::uint64_t MorganHash(
  const std::vector<std::uint64_t>& vertex_hashes,
  const boost::dynamic_bitset<>& vertex_mask,
  const Graph& graph,
  std::uint8_t n_updates,
  const std::function<
    std::uint64_t(const Edge&)>& edge_hasher = nullptr) {
  // Note that this hash is NOT the same (neither in value nor conceptually) as
  // the RDKit's Morgan fingerprints. Our hashes only capture information of a
  // specified subgraph, whereas the RDKit's implementation captures information
  // of farther atoms.
  if (vertex_mask.none()) {
    return 0;
  };
  std::vector<std::uint64_t> morgan_hashes (vertex_hashes);
  for (std::uint8_t n = 0; n < n_updates; ++n) {
    MorganHashUpdate(morgan_hashes, vertex_mask, graph, edge_hasher);
  };
  std::vector<std::uint64_t> mask_hashes;
  mask_hashes.reserve(vertex_mask.count());
  for (std::size_t vertex = vertex_mask.find_first();
    vertex != boost::dynamic_bitset<>::npos;
    vertex = vertex_mask.find_next(vertex)) {
    mask_hashes.push_back(morgan_hashes[vertex]);
  };
  std::sort(mask_hashes.begin(), mask_hashes.end());
  return boost::hash_range(mask_hashes.begin(), mask_hashes.end());
};


template <class Graph, class Edge>
std::uint64_t CollapseHash(
  std::size_t vertex,
  const std::vector<std::uint64_t>& vertex_hashes,
  const std::vector<std::uint8_t>& distances_to_root,
  const boost::dynamic_bitset<>& vertex_mask,
  const Graph& graph,
  const std::function<
    std::uint64_t(const Edge&)>& edge_hasher = nullptr) {
  // Collapses the hashes of neighboring vertices at a distance greater or equal
  // to that of the root vertex onto the vertex's hash.
  // We include equidistant neighbors to capture intra-mask cycles.
  std::uint8_t vertex_distance = distances_to_root[vertex];
  std::vector<std::uint64_t> neighbor_hashes;
  neighbor_hashes.reserve(boost::out_degree(vertex, graph));
  for (std::size_t neighbor : boost::make_iterator_range(
    boost::adjacent_vertices(vertex, graph))) {
    std::uint8_t neighbor_distance = distances_to_root[neighbor];
    if (!vertex_mask[neighbor] || neighbor_distance < vertex_distance) {
      continue;
    };
    std::uint64_t neighbor_hash = vertex_hashes[neighbor];
    if (edge_hasher) {
      auto [edge, _] = boost::edge(vertex, neighbor, graph);
      boost::hash_combine(neighbor_hash, edge_hasher(graph[edge]));
    };
    neighbor_hashes.push_back(neighbor_hash);
  };
  std::uint64_t collapsed_hash = vertex_hashes[vertex];
  std::sort(neighbor_hashes.begin(), neighbor_hashes.end());
  boost::hash_combine(collapsed_hash,
    boost::hash_range(neighbor_hashes.begin(), neighbor_hashes.end()));
  return collapsed_hash;
};

template <class Graph, class Edge>
std::uint64_t CollapsingHash(
  std::size_t root,
  const std::vector<std::uint64_t>& vertex_hashes,
  const std::vector<std::uint8_t>& distances_to_root,
  const boost::dynamic_bitset<>& vertex_mask,
  const Graph& graph,
  const std::function<
    std::uint64_t(const Edge&)>& edge_hasher = nullptr) {
  // Similar to a Morgan hash, but the final value depends on the root vertex.
  // Hence two identical subgraphs may have different hashes depending on which 
  // vertex they are rooted.
  if (vertex_mask.none()) {
    return 0;
  };
  // Sort the vertices by their distance to the root vertex.
  std::vector<std::pair<std::size_t, std::uint8_t>> vertices_by_distance;
  vertices_by_distance.reserve(vertex_mask.count());
  for (std::size_t vertex = vertex_mask.find_first();
    vertex != boost::dynamic_bitset<>::npos;
    vertex = vertex_mask.find_next(vertex)) {
    vertices_by_distance.emplace_back(vertex, distances_to_root[vertex]);
  };
  std::sort(vertices_by_distance.begin(), vertices_by_distance.end(),
    [](const auto& p1, const auto& p2) {
      return p1.second > p2.second;
    }
  );
  // Iterate over the vertices by descending distance and update their hashes by
  // collapsing the hashes of neighboring vertices onto them.
  std::uint8_t prev_distance = vertices_by_distance[0].second;
  std::vector<std::uint64_t> prev_vertex_hashes (vertex_hashes);
  std::vector<std::uint64_t> current_vertex_hashes (vertex_hashes);
  for (const auto& [vertex, distance] : vertices_by_distance) {
    // We update the hashes of all equidistant vertices simultaneously. Since
    // equidistant vertices influence each others hashes (see CollapseHash())
    // this ensures that the update is order independent.
    if (distance < prev_distance) {
      prev_vertex_hashes = current_vertex_hashes;
      prev_distance = distance;
    };
    current_vertex_hashes[vertex] = CollapseHash(vertex,
      prev_vertex_hashes, distances_to_root, vertex_mask, graph, edge_hasher);
  };
  return current_vertex_hashes[root];
};


std::uint64_t BondTypeAsHash(const RDKit::Bond* bond) {
  return bond->getBondType();
};


template <>
struct std::hash<RDKit::ROMol> {
  std::size_t operator()(const RDKit::ROMol& molecule) const {
    std::vector<std::uint32_t> invariants = Invariants(molecule);
    RDKit::SparseIntVect<std::uint32_t>* fingerprint =
      RDKit::MorganFingerprints::getFingerprint(molecule, 2, &invariants);
    std::size_t hash = 0;
    for (auto [feature_id, count] : fingerprint->getNonzeroElements()) {
      boost::hash_combine(hash, feature_id);
      boost::hash_combine(hash, count);
    };
    delete fingerprint;
    return hash;
  };

  std::vector<std::uint32_t> Invariants(const RDKit::ROMol& molecule) const {
    // RDKit::ROMols store pre-computed properties. These properties (including
    // number of implicit hydrogens, valence, aromaticity, stereochemistry and
    // ring membership) may be invalidated upon molecule edition. However, the
    // properties aren't necessarily erased and may be accessed by other down-
    // stream functions. The user is supposed to sanitize the molecule after
    // modifying it, which recalculates these properties, but doing so in a loop
    // is expensive and finicky since the calculation may fail. Instead we avoid
    // using these properties altogether and roll our own invariants.
    std::size_t n_atoms = molecule.getNumAtoms();
    std::vector<std::uint32_t> invariants (n_atoms);
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      invariants[atom_idx] = AtomKeyHash(molecule.getAtomWithIdx(atom_idx));
    };
    return invariants;
  };
};


std::size_t HashMolecule(const RDKit::ROMol& molecule) {
  static std::hash<RDKit::ROMol> molecule_hasher;
  return molecule_hasher(molecule);
};

#endif // !_MOLECULE_HASH_HPP_