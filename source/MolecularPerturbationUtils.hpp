#pragma once
#ifndef _MOLECULAR_PERTURBATION_UTILS_HPP_
#define _MOLECULAR_PERTURBATION_UTILS_HPP_

#include "Valence.hpp"
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <bitset>

struct ElementAromaticitySpecification {
  std::bitset<119> aromaticity_mask;

  ElementAromaticitySpecification() {
    // Potentially aromatic elements according to the OpenSMILES specification.
    aromaticity_mask.set(5);
    aromaticity_mask.set(6);
    aromaticity_mask.set(7);
    aromaticity_mask.set(8);
    aromaticity_mask.set(15);
    aromaticity_mask.set(16);
    aromaticity_mask.set(33);
    aromaticity_mask.set(34);
  };

  bool operator()(std::size_t atomic_number) const {
    return aromaticity_mask[atomic_number];
  };
};

static const ElementAromaticitySpecification
  open_smiles_aromaticity_specification;

void CorrectElementAromaticity(
  RDKit::RWMol& molecule,
  const ElementAromaticitySpecification&
    aromaticity_specification = open_smiles_aromaticity_specification) {
  // Molecule generators may label any ring atom as aromatic, regardless of the
  // atomic number. This isn't necessarily an issue. However, the OpenSMILES
  // specification only allows certain elements to be aromatic. Hence, the
  // SMILES of our designed molecules may not be parseable. If you need to use
  // SMILES later on you have two options: try to kekulize the molecule (and
  // probably fail) or simply label the problematic atoms as non-aromatic.
  for (RDKit::Atom* atom : molecule.atoms()) {
    if (atom->getIsAromatic() &&
      !aromaticity_specification(atom->getAtomicNum())) {
      atom->setIsAromatic(false);
    };
  };
};

void CorrectHydrogenCount(
  RDKit::RWMol& molecule,
  bool update_property_cache = true) {
  for (RDKit::Atom* atom : molecule.atoms()) {
    auto [n_hydrogens, can_be_corrected] = MinNumHydrogensForValidValence(atom);
    if (can_be_corrected) {
      atom->setNumExplicitHs(n_hydrogens);
    };
  };
  if (update_property_cache) {
    molecule.updatePropertyCache();
  };
};

void PartialSanitization(
  RDKit::RWMol& molecule,
  bool kekulize = false,
  bool aromatize = true) {
  // Sanitize the molecule without strict valence checks.
  unsigned sanitization_flags =
    RDKit::MolOps::SanitizeFlags::SANITIZE_ALL ^
    RDKit::MolOps::SanitizeFlags::SANITIZE_PROPERTIES;
  // Some incorrect aromatic molecules aren't kekulizable. It may be desirable
  // to avoid kekulization altogether. Note that if we don't kekulize we can't
  // search for radicals either.
  if (!kekulize) {
    sanitization_flags ^= RDKit::MolOps::SanitizeFlags::SANITIZE_KEKULIZE;
    sanitization_flags ^= RDKit::MolOps::SanitizeFlags::SANITIZE_FINDRADICALS;
  };
  // If we are dealing with unkekulizable molecules we might not want to
  // attempt to assign their aromaticity either.
  if (!aromatize) {
    sanitization_flags ^= RDKit::MolOps::SanitizeFlags::SANITIZE_SETAROMATICITY;
    sanitization_flags ^= RDKit::MolOps::SanitizeFlags::SANITIZE_SETCONJUGATION;
  };
  unsigned failure_flag;
  RDKit::MolOps::sanitizeMol(molecule, failure_flag, sanitization_flags);
};

RDKit::RWMOL_SPTR PartiallySanitizedMoleculeFromSMILES(const std::string& smiles) {
  // Parse the SMILES without sanitizing the resulting molecule.
  RDKit::RWMOL_SPTR molecule (RDKit::SmilesToMol(smiles, 0, false));
  // SMILES parsing can fail due to a variety of reasons.
  if (!molecule) {
    throw RDKit::SmilesParseException("Couldn't parse SMILES " + smiles);
  };
  // Sanitize the molecule.
  PartialSanitization(*molecule);
  return molecule;
};

std::string UnsanitizedMoleculeSMILES(const RDKit::ROMol& molecule) {
  // RDKit::ROMols store pre-computed properties. These properties (including
  // number of implicit hydrogens, valence, aromaticity, stereochemistry and
  // ring membership) may be invalidated upon molecule edition. However, the
  // properties aren't necessarily erased and may be accessed by other down-
  // stream functions. The user is supposed to sanitize the molecule after
  // modifying it, which recalculates these properties, but doing so in a loop
  // is expensive and finicky since the calculation may fail. Instead we avoid
  // using these properties altogether.
  return RDKit::MolToSmiles(
    molecule,
    false,  // Don't include stereochemistry since it's likely invalid.
    false,  // Don't kekulize since it requires aromaticity and valence info.
    -1,
    false); // Don't canonizalize since that requires valid stereochemistry info.
};

#endif // !_MOLECULAR_PERTURBATION_UTILS_HPP_
