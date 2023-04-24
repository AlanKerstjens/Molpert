#pragma once
#ifndef _CIRCULAR_ATOMIC_ENVIRONMENT_HPP_
#define _CIRCULAR_ATOMIC_ENVIRONMENT_HPP_

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "MoleculeHash.hpp"
#include "MolecularGraphSearch.hpp"

typedef std::uint64_t EnvironmentKey;
typedef std::pair<std::uint64_t, std::uint64_t> EnvironmentKeyChange;
static const EnvironmentKey NULL_ENVIRONMENT_KEY (0);

struct CircularAtomicEnvironment {
  const RDKit::ROMol* molecule;
  const RDKit::Atom* root_atom;
  const std::uint8_t radius;
  const std::size_t n_atoms; // In molecule, not in environment.
  boost::dynamic_bitset<> atom_mask;
  boost::dynamic_bitset<> bond_mask;
  std::vector<std::uint8_t> atom_distances;

public:
  CircularAtomicEnvironment(const RDKit::Atom* root_atom, std::uint8_t radius) :
    molecule(&root_atom->getOwningMol()), root_atom(root_atom),
    radius(radius),
    n_atoms(molecule->getNumAtoms()),
    atom_mask(n_atoms), bond_mask(molecule->getNumBonds()),
    atom_distances(n_atoms, 0xFF) {
    MolecularGraphDiscoveryJournal journal = MoleculeBFS(
      *molecule,
      root_atom->getIdx(),
      radius,
      boost::dynamic_bitset<>(atom_mask.size()),
      boost::dynamic_bitset<>(bond_mask.size()),
      NeverMetCriterion,
      true);
    atom_mask = journal.GetDiscoveredAtomsMask();
    bond_mask = journal.GetDiscoveredBondsMask();
    for (auto [atom_idx, distance] : journal.GetAtomDiscoveries()) {
      atom_distances[atom_idx] = distance;
    };
  };

  std::uint64_t Hash(const std::vector<std::uint64_t>& atom_hashes) const {
    return MorganHash<RDKit::MolGraph, RDKit::Bond*>(
      atom_hashes,
      atom_mask,
      molecule->getTopology(),
      (radius + 1) / 2,
      BondTypeAsHash);
  };

  EnvironmentKey Key(const std::vector<std::uint64_t>& atom_hashes) const {
    return Hash(atom_hashes);
  };

  EnvironmentKey Key() const {
    return Hash(AtomKeyHashes(*molecule));
  };

  std::vector<std::size_t> AtomsAtDistance(std::uint8_t distance) const {
    assert(distance <= radius);
    std::vector<std::size_t> atom_indices;
    atom_indices.reserve(atom_mask.count());
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      if (atom_distances[atom_idx] == distance) {
        atom_indices.push_back(atom_idx);
      };
    };
    return atom_indices;
  };

  std::string SMILES() const {
    RDKit::RWMol environment;
    std::unordered_map<std::size_t, std::size_t> atom_map;
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      if (atom_mask[atom_idx]) {
        // atom shouldn't be const as this would cause an incorrect 
        // RWMol::addAtom overload resolution due to a pointer to bool implicit
        // conversion. However, I can't legally get a non-const RDKit::Atom* 
        // because molecule is const. I can guarantee this code doesn't modify 
        // molecule so I'll cast the constness away.
        RDKit::Atom* atom = const_cast<RDKit::Atom*>(
          molecule->getAtomWithIdx(atom_idx));
        std::size_t env_atom_idx = environment.addAtom(atom, false, false);
        atom_map.emplace(atom_idx, env_atom_idx);
      };
    };
    for (std::size_t bond_idx = 0; bond_idx < bond_mask.size(); ++bond_idx) {
      if (bond_mask[bond_idx]) {
        const RDKit::Bond* bond = molecule->getBondWithIdx(bond_idx);
        environment.addBond(
          atom_map[bond->getBeginAtomIdx()], atom_map[bond->getEndAtomIdx()],
          bond->getBondType());
      };
    };
    return RDKit::MolToSmiles(environment);
  };
};


class CircularAtomicEnvironmentGenerator {
  std::uint8_t environment_radius;
  AtomsHasher atoms_hasher;
  const RDKit::ROMol* hashed_molecule = nullptr; // Cache
  std::vector<std::uint64_t> atom_hashes; // Cache

private:
  void CalculateAtomHashes(const RDKit::ROMol* molecule) {
    atom_hashes = atoms_hasher(*molecule);
    hashed_molecule = molecule;
  };

public:
  CircularAtomicEnvironmentGenerator(
    std::uint8_t environment_radius = 2,
    const AtomsHasher& atoms_hasher = AtomKeyHashes) :
    environment_radius(environment_radius),
    atoms_hasher(atoms_hasher) {};

  CircularAtomicEnvironment operator()(const RDKit::Atom* atom) const {
    return CircularAtomicEnvironment(atom, environment_radius);
  };

  EnvironmentKey Key(const RDKit::Atom* atom) {
    CircularAtomicEnvironment environment (atom, environment_radius);
    if (environment.molecule != hashed_molecule) {
      CalculateAtomHashes(environment.molecule);
    };
    return environment.Hash(atom_hashes);
  };

  std::uint8_t GetEnvironmentRadius() const {
    return environment_radius;
  };
};


std::pair<std::vector<unsigned>, std::vector<unsigned>> 
CircularEnvironmentOverlap(
  const std::vector<CircularAtomicEnvironment>& environments) {
  const CircularAtomicEnvironment& representative_environment = 
    environments.at(0);
  const RDKit::ROMol* molecule = representative_environment.molecule;
  unsigned radius = representative_environment.radius;
  std::size_t n_atoms = molecule->getNumAtoms();
  std::size_t n_bonds = molecule->getNumBonds();
  std::vector<unsigned> atom_overlap (n_atoms, 0);
  std::vector<unsigned> bond_overlap (n_bonds, 0);
  for (const CircularAtomicEnvironment& environment : environments) {
    assert(molecule == environment.molecule && radius == environment.radius);
    for (std::size_t atom_idx = environment.atom_mask.find_first();
      atom_idx != boost::dynamic_bitset<>::npos;
      atom_idx = environment.atom_mask.find_next(atom_idx)) {
      ++atom_overlap[atom_idx];
    };
    for (std::size_t bond_idx = environment.bond_mask.find_first();
      bond_idx != boost::dynamic_bitset<>::npos;
      bond_idx = environment.bond_mask.find_next(bond_idx)) {
      ++bond_overlap[bond_idx];
    };
  };
  return {std::move(atom_overlap), std::move(bond_overlap)};
};

#endif // !_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_
