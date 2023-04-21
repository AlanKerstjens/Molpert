#pragma once
#ifndef _ADD_NITRO_HPP_
#define _ADD_NITRO_HPP_

#include "MolecularPerturbations.hpp"

class AddNitro : public TopologicalPerturbation {
public:
  AddNitro(const RDKit::ROMol& molecule, AtomIdx atom_idx) {
    std::size_t n_atoms = molecule.getNumAtoms();
    AtomIdx max_atom_idx = n_atoms ? n_atoms - 1 : 0;
    Tag max_atom_tag = GetMaxAtomTag(molecule);
    Tag max_bond_tag = GetMaxBondTag(molecule);
    atom_constructions.emplace(++max_atom_tag, 7, 1);  // N+, max_atom_idx + 1
    atom_constructions.emplace(++max_atom_tag, 8, -1); // O-, max_atom_idx + 2
    atom_constructions.emplace(++max_atom_tag, 8);     // O, max_atom_idx + 3
    bond_constructions.emplace(++max_bond_tag,
      atom_idx, max_atom_idx + 1, RDKit::Bond::SINGLE);
    bond_constructions.emplace(++max_bond_tag,
      max_atom_idx + 1, max_atom_idx + 2, RDKit::Bond::SINGLE);
    bond_constructions.emplace(++max_bond_tag,
      max_atom_idx + 1, max_atom_idx + 3, RDKit::Bond::DOUBLE);
  };
};

#endif // !_ADD_NITRO_HPP_
