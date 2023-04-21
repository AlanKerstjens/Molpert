#pragma once
#ifndef _MOLECULAR_TAGS_HPP_
#define _MOLECULAR_TAGS_HPP_

#include <GraphMol/ROMol.h>

// The RDKit uses std::vector as the underlying representation for the molecular
// graph adjacency list. This means that atom and bond deletions invalidate atom
// and bond indices respectively. To keep track of atoms and bonds despite their
// index invalidation we associate them with a tag.

typedef std::size_t Tag;
static const std::string tag_property_name = "Tag";
static const std::string max_atom_tag_property_name = "XAT";
static const std::string max_bond_tag_property_name = "XBT";

bool AtomsAreTagged(const RDKit::ROMol& molecule) {
  return molecule.hasProp(max_atom_tag_property_name);
};

bool BondsAreTagged(const RDKit::ROMol& molecule) {
  return molecule.hasProp(max_bond_tag_property_name);
};

Tag GetTag(const RDKit::Atom* atom) {
  return atom->getProp<Tag>(tag_property_name);
};

Tag GetTag(const RDKit::Bond* bond) {
  return bond->getProp<Tag>(tag_property_name);
};

Tag GetAtomTag(const RDKit::ROMol& molecule, std::size_t atom_idx) {
  return GetTag(molecule.getAtomWithIdx(atom_idx));
};

Tag GetBondTag(const RDKit::ROMol& molecule, std::size_t bond_idx) {
  return GetTag(molecule.getBondWithIdx(bond_idx));
};

std::pair<Tag, bool> GetTagIfPresent(const RDKit::Atom* atom) {
  Tag atom_tag = 0;
  bool present = atom->getPropIfPresent<Tag>(tag_property_name, atom_tag);
  return {atom_tag, present};
};

std::pair<Tag, bool> GetTagIfPresent(const RDKit::Bond* bond) {
  Tag bond_tag = 0;
  bool present = bond->getPropIfPresent<Tag>(tag_property_name, bond_tag);
  return {bond_tag, present};
};

Tag GetMaxAtomTag(const RDKit::ROMol& molecule) {
  return molecule.getProp<Tag>(max_atom_tag_property_name);
};

Tag GetMaxBondTag(const RDKit::ROMol& molecule) {
  return molecule.getProp<Tag>(max_bond_tag_property_name);
};

std::pair<Tag, bool> GetMaxAtomTag(
  RDKit::ROMol& molecule,
  Tag tag_if_missing = 0) {
  Tag max_atom_tag = tag_if_missing;
  bool present = molecule.getPropIfPresent(
    max_atom_tag_property_name, max_atom_tag);
  if (!present) {
    molecule.setProp<Tag>(max_atom_tag_property_name, tag_if_missing);
  };
  return {max_atom_tag, present};
};

std::pair<Tag, bool> GetMaxBondTag(
  RDKit::ROMol& molecule,
  Tag tag_if_missing = 0) {
  Tag max_bond_tag = tag_if_missing;
  bool present = molecule.getPropIfPresent(
    max_bond_tag_property_name, max_bond_tag);
  if (!present) {
    molecule.setProp<Tag>(max_bond_tag_property_name, tag_if_missing);
  };
  return {max_bond_tag, present};
};

// There is no computationally efficient way to retrieve an atom or bond with
// a specific tag. The best we can do is loop over all atoms or bonds until
// we find the appropiate one. Be mindful of this inefficiency and use these
// functions sparingly.
const RDKit::Atom* GetAtomWithTag(const RDKit::ROMol& molecule, Tag atom_tag) {
  for (const RDKit::Atom* atom : molecule.atoms()) {
    if (atom->getProp<Tag>(tag_property_name) == atom_tag) {
      return atom;
    };
  };
  return nullptr;
};

const RDKit::Bond* GetBondWithTag(const RDKit::ROMol& molecule, Tag bond_tag) {
  for (const RDKit::Bond* bond : molecule.bonds()) {
    if (bond->getProp<Tag>(tag_property_name) == bond_tag) {
      return bond;
    };
  };
  return nullptr;
};

Tag TagAtom(RDKit::ROMol& molecule, RDKit::Atom* atom) {
  auto [max_atom_tag, was_present] = GetMaxAtomTag(molecule, 0);
  if (was_present) {
    atom->setProp<Tag>(tag_property_name, ++max_atom_tag);
    molecule.setProp<Tag>(max_atom_tag_property_name, max_atom_tag);
  } else {
    atom->setProp<Tag>(tag_property_name, max_atom_tag);
  };
  return max_atom_tag;
};

void TagAtom(RDKit::ROMol& molecule, RDKit::Atom* atom, Tag atom_tag) {
  auto [max_atom_tag, was_present] = GetMaxAtomTag(molecule, atom_tag);
  if (was_present) {
    if (atom_tag <= max_atom_tag) {
      throw std::runtime_error("Atom tag may have been used previously.");
    };
    molecule.setProp<Tag>(max_atom_tag_property_name, atom_tag);
  };
  atom->setProp<Tag>(tag_property_name, atom_tag);
};

Tag TagBond(RDKit::ROMol& molecule, RDKit::Bond* bond) {
  auto [max_bond_tag, was_present] = GetMaxBondTag(molecule, 0);
  if (was_present) {
    bond->setProp<Tag>(tag_property_name, ++max_bond_tag);
    molecule.setProp<Tag>(max_bond_tag_property_name, max_bond_tag);
  } else {
    bond->setProp<Tag>(tag_property_name, max_bond_tag);
  };
  return max_bond_tag;
};

void TagBond(RDKit::ROMol& molecule, RDKit::Bond* bond, Tag bond_tag) {
  auto [max_bond_tag, was_present] = GetMaxBondTag(molecule, bond_tag);
  if (was_present) {
    if (bond_tag <= max_bond_tag) {
      throw std::runtime_error("Bond tag may have been used previously.");
    };
    molecule.setProp<Tag>(max_bond_tag_property_name, bond_tag);
  };
  bond->setProp<Tag>(tag_property_name, bond_tag);
};

void TagAtoms(
  RDKit::ROMol& molecule,
  bool skip_if_tagged = false,
  Tag start_tag = 0) {
  if (skip_if_tagged && AtomsAreTagged(molecule)) {
    return;
  };
  for (RDKit::Atom* atom : molecule.atoms()) {
    atom->setProp<Tag>(tag_property_name, start_tag++);
  };
  molecule.setProp<Tag>(max_atom_tag_property_name,
    start_tag ? start_tag - 1 : 0);
};

void TagBonds(
  RDKit::ROMol& molecule,
  bool skip_if_tagged = false,
  Tag start_tag = 0) {
  if (skip_if_tagged && BondsAreTagged(molecule)) {
    return;
  };
  for (RDKit::Bond* bond : molecule.bonds()) {
    bond->setProp<Tag>(tag_property_name, start_tag++);
  };
  molecule.setProp<Tag>(max_bond_tag_property_name,
    start_tag ? start_tag - 1 : 0);
};

void TagMolecule(
  RDKit::ROMol& molecule,
  bool skip_if_tagged = false,
  Tag start_atom_tag = 0,
  Tag start_bond_tag = 0) {
  TagAtoms(molecule, skip_if_tagged, start_atom_tag);
  TagBonds(molecule, skip_if_tagged, start_bond_tag);
};

void UntagAtom(RDKit::Atom* atom) {
  atom->clearProp(tag_property_name);
};

void UntagBond(RDKit::Bond* bond) {
  bond->clearProp(tag_property_name);
};

void UntagAtoms(RDKit::ROMol& molecule) {
  for (RDKit::Atom* atom : molecule.atoms()) {
    atom->clearProp(tag_property_name);
  };
  molecule.clearProp(max_atom_tag_property_name);
};

void UntagBonds(RDKit::ROMol& molecule) {
  for (RDKit::Bond* bond : molecule.bonds()) {
    bond->clearProp(tag_property_name);
  };
  molecule.clearProp(max_bond_tag_property_name);
};

void UntagMolecule(RDKit::ROMol& molecule) {
  UntagAtoms(molecule);
  UntagBonds(molecule);
};

#endif // !_MOLECULAR_TAGS_HPP_
