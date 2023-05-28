#ifndef _SUBGRAPH_HPP_
#define _SUBGRAPH_HPP_

#include "Valence.hpp"
#include "RandomSampling.hpp"
#include "MolecularGraphSearch.hpp"

struct AttachmentPoint {
  enum class Definition {AllAtoms, Valence, BrokenBonds};
  std::size_t atom_idx = 0;
  std::size_t type = 0;
};


struct SubgraphRequests {
  // These settings are merely requests. An attempt is made to fulfill them,
  // but no promises are made regarding the actual outcome.
  AttachmentPoint::Definition attachment_point_definition =
    AttachmentPoint::Definition::BrokenBonds;
  std::size_t n_attachment_points = std::numeric_limits<std::size_t>::max();
  std::size_t size = std::numeric_limits<std::size_t>::max();
  bool allow_breaking_cycles = false;
  bool one_attachment_point_per_slot = false;
  bool pad_attachment_points_if_insufficient = false;
  bool try_to_satisfy_valences_when_padding = true;
  bool shuffle_attachment_points = true;
};


std::vector<AttachmentPoint> AllAtomsAsAttachmentPoints(
  const boost::dynamic_bitset<>& atoms_mask) {
  std::vector<AttachmentPoint> attachment_points;
  attachment_points.reserve(atoms_mask.count());
  for (std::size_t atom_idx = atoms_mask.find_first();
    atom_idx != boost::dynamic_bitset<>::npos;
    atom_idx = atoms_mask.find_next(atom_idx)) {
    attachment_points.emplace_back(atom_idx, 1);
  };
  return attachment_points;
};

std::vector<AttachmentPoint> ValenceAttachmentPoints(
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask,
  bool available_valence_as_type = false,
  bool one_per_available_valence = false) {
  std::vector<AttachmentPoint> attachment_points;
  attachment_points.reserve(atoms_mask.count());
  for (std::size_t atom_idx = atoms_mask.find_first();
    atom_idx != boost::dynamic_bitset<>::npos;
    atom_idx = atoms_mask.find_next(atom_idx)) {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    int available_valence = AvailableValence(atom);
    if (available_valence <= 0) {
      continue;
    };
    // Available valence can be up to INT_MAX. We cap it to 4 to avoid excessive
    // memory consumption.
    available_valence = available_valence > 4 ? 4 : available_valence;
    std::size_t type = available_valence_as_type ? available_valence : 1;
    std::size_t n = one_per_available_valence ? available_valence : 1;
    for (std::size_t i = 0; i < n; ++i) {
      attachment_points.emplace_back(atom_idx, type);
    };
  };
  return attachment_points; 
};

std::vector<AttachmentPoint> BrokenBondsAttachmentPoints(
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask,
  bool bond_type_as_type = false,
  bool one_per_bond_order = false) {
  std::vector<AttachmentPoint> attachment_points;
  for (std::size_t atom_idx = atoms_mask.find_first();
    atom_idx != boost::dynamic_bitset<>::npos;
    atom_idx = atoms_mask.find_next(atom_idx)) {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
      std::size_t neighbor_idx = bond->getOtherAtomIdx(atom_idx);
      if (atoms_mask[neighbor_idx]) {
        continue;
      };
      std::size_t type = bond_type_as_type ? bond->getBondType() : 1;
      std::size_t n = one_per_bond_order ? bond->getBondTypeAsDouble() : 1;
      for (std::size_t i = 0; i < n; ++i) {
        attachment_points.emplace_back(atom_idx, type);
      };
    };
  };
  return attachment_points;
};


void PadAttachmentPoints(
  std::vector<AttachmentPoint>& attachment_points,
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask,
  const SubgraphRequests& requests,
  std::mt19937& prng) {
  if (attachment_points.size() >= requests.n_attachment_points || 
    requests.n_attachment_points == std::numeric_limits<std::size_t>::max() ||
    atoms_mask.none()) {
    return;
  };
  std::size_t n_extra_attachment_points = 
    requests.n_attachment_points - attachment_points.size();
  // We might have to restrict our choice to atoms having available valences.
  if (requests.try_to_satisfy_valences_when_padding) {
    std::vector<std::size_t> eligible_atom_indices;
    eligible_atom_indices.reserve(atoms_mask.count());
    for (std::size_t atom_idx = atoms_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = atoms_mask.find_next(atom_idx)) {
      const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
      int available_valence = AvailableValence(atom);
      available_valence = available_valence > 4 ? 4 : available_valence;
      for (int i = 0; i < available_valence; ++i) {
        eligible_atom_indices.push_back(atom_idx);
      };
    };
    // If there aren't sufficient atoms with available valences we sample 
    // atoms randomly.
    std::size_t n_eligible_atoms = eligible_atom_indices.size();
    if (n_eligible_atoms < n_extra_attachment_points) {
      for (std::size_t i = 0; 
        i < n_extra_attachment_points - n_eligible_atoms; ++i) {
        eligible_atom_indices.push_back(Sample(atoms_mask, prng));
      };
    };
    // Randomly pick atoms from the eligible ones.
    std::shuffle(
      eligible_atom_indices.begin(), eligible_atom_indices.end(), prng);
    for (std::size_t i = 0; i < n_extra_attachment_points; ++i) {
      attachment_points.emplace_back(eligible_atom_indices[i], 1);
    };
  // If not, pick the atoms randomly.
  } else {
    for (std::size_t i = 0; i < n_extra_attachment_points; ++i) {
      attachment_points.emplace_back(Sample(atoms_mask, prng), 1);
    };
  };
};


std::vector<AttachmentPoint> AttachmentPoints(
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask,
  const SubgraphRequests& requests,
  std::mt19937& prng) {
  // Generate attachment points based on the provided definition.
  std::vector<AttachmentPoint> attachment_points;
  switch (requests.attachment_point_definition) {
    case AttachmentPoint::Definition::AllAtoms:
      attachment_points = AllAtomsAsAttachmentPoints(atoms_mask);
      break;
    case AttachmentPoint::Definition::Valence:
      attachment_points = ValenceAttachmentPoints(
        molecule, atoms_mask, false, requests.one_attachment_point_per_slot);
      break;
    case AttachmentPoint::Definition::BrokenBonds:
      attachment_points = BrokenBondsAttachmentPoints(
        molecule, atoms_mask, false, requests.one_attachment_point_per_slot);
  };
  // If a specific number of attachment points was request but not generated we
  // may add extra attachment points by randomly sampling atoms.
  if (requests.pad_attachment_points_if_insufficient) {
    PadAttachmentPoints(attachment_points, molecule, atoms_mask, requests, prng);
  };
  if (requests.shuffle_attachment_points) {
    std::shuffle(attachment_points.begin(), attachment_points.end(), prng);
  };
  return attachment_points;
};


boost::dynamic_bitset<> RandomSubgraph(
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask,
  const SubgraphRequests& requests,
  std::mt19937& prng) {
  if (atoms_mask.none()) {
    return boost::dynamic_bitset<>();
  };
  // When attachment points are defined at broken bonds finding a subgraph with
  // a desired number of attachment points is challenging. If all broken bonds 
  // are acyclic it is likely that the resulting fragment has at least N 
  // attachment points if the subgraph covers an atom of degree N or larger. The 
  // exception to this rule is if the walk covers a whole branch of the graph.
  boost::dynamic_bitset<> start_atoms_mask (atoms_mask);
  if (requests.attachment_point_definition == 
   AttachmentPoint::Definition::BrokenBonds &&
    requests.n_attachment_points < std::numeric_limits<std::size_t>::max()) {
    // If rings aren't broken walking through them will neither 
    // increase nor decrease the number of attachment points.
    if (!requests.allow_breaking_cycles) {
      const RDKit::RingInfo* ring_info = molecule.getRingInfo();
      if (!ring_info->isInitialized()) {
        RDKit::MolOps::findSSSR(molecule);
      };
      for (std::size_t atom_idx = start_atoms_mask.find_first();
        atom_idx != boost::dynamic_bitset<>::npos;
        atom_idx = start_atoms_mask.find_next(atom_idx)) {
        if (ring_info->numAtomRings(atom_idx)) {
          start_atoms_mask.reset(atom_idx);
        };
      };
    };
    if (start_atoms_mask.none()) {
      start_atoms_mask = atoms_mask;
    };
    // Dangling attachment points are acceptable, but undesirable. Hence, as a 
    // heuristic to maximize the probability of generating a fragment with N 
    // attachment points we start the walk at an atom with degree close to N.
    std::vector<std::pair<std::size_t, int>> ddiffs; // Degree differences
    ddiffs.reserve(start_atoms_mask.count());
    for (std::size_t atom_idx = start_atoms_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = start_atoms_mask.find_next(atom_idx)) {
      const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
      ddiffs.emplace_back(
        atom_idx, atom->getDegree() - requests.n_attachment_points);
    };
    std::sort(ddiffs.begin(), ddiffs.end(),
      [](const auto& p1, const auto& p2) {
        if (p1.second >= 0 && p2.second >= 0) {
          return p1.second < p2.second;
        };
        return p1.second > p2.second;
      });
    // If some atoms have degree N we only consider those.
    int best_degree_difference = ddiffs.front().second;
    for (const auto& [atom_idx, degree_difference] : ddiffs) {
      if (degree_difference != best_degree_difference) {
        start_atoms_mask.reset(atom_idx);
      };
    };
  };

  MolecularGraphDiscoveryJournal journal = MoleculeRandomWalk(
    molecule, Sample(start_atoms_mask, prng), requests.size, prng, 
    ~atoms_mask, boost::dynamic_bitset<>(molecule.getNumBonds()),
    !requests.allow_breaking_cycles);

  return journal.GetDiscoveredAtomsMask();
};

boost::dynamic_bitset<> RandomSubgraph(
  const RDKit::ROMol& molecule,
  const SubgraphRequests& requests,
  std::mt19937& prng) {
  boost::dynamic_bitset<> atoms_mask (molecule.getNumAtoms());
  atoms_mask.set();
  return RandomSubgraph(molecule, std::move(atoms_mask), requests, prng);
};


RDKit::RWMol SubgraphToMolecule(
  const RDKit::ROMol& molecule,
  const boost::dynamic_bitset<>& atoms_mask) {
  RDKit::RWMol subgraph_molecule (molecule);
  std::size_t n_atoms = molecule.getNumAtoms();
  for (std::size_t atom_idx = n_atoms; atom_idx-- > 0;) {
    if (!atoms_mask[atom_idx]) {
      subgraph_molecule.removeAtom(atom_idx);
    };
  };
  return subgraph_molecule;
};

#endif // _SUBGRAPH_HPP_