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
    if (!aromaticity_specification(atom->getAtomicNum())) {
      atom->setIsAromatic(false);
    };
  };
};

void Aromatize(
  RDKit::RWMol& molecule,
  const boost::dynamic_bitset<>& bond_mask) {
  for (std::size_t bond_idx = bond_mask.find_first();
    bond_idx != boost::dynamic_bitset<>::npos;
    bond_idx = bond_mask.find_next(bond_idx)) {
    RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    RDKit::Atom* begin_atom = bond->getBeginAtom();
    RDKit::Atom* end_atom = bond->getEndAtom();
    bond->setBondType(RDKit::Bond::AROMATIC);
    bond->setIsAromatic(true);
    begin_atom->setIsAromatic(true);
    end_atom->setIsAromatic(true);
  };
};

std::pair<RDKit::Bond*, int> LeastAvailableBond(
  RDKit::ROMol& molecule,
  const std::vector<int>& available_valences,
  const boost::dynamic_bitset<>& bond_mask) {
  RDKit::Bond* lav_bond = nullptr;
  int lav = std::numeric_limits<int>::max();
  // Note for the future: if available valence ties are causing problems the
  // ties could be resolve with the highest atom degree.
  for (std::size_t bond_idx = bond_mask.find_first();
    bond_idx != boost::dynamic_bitset<>::npos;
    bond_idx = bond_mask.find_next(bond_idx)) {
    RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    int bav = available_valences[bond->getBeginAtomIdx()];
    int eav = available_valences[bond->getEndAtomIdx()];
    int av = bav < eav ? bav : eav;
    // Bonds with negative available valence can't be handled.
    if (av < 0) {
      continue;
    };
    if (av < lav) {
      lav_bond = bond;
      lav = av;
    };
  };
  return {lav_bond, lav};
};

boost::dynamic_bitset<> NextBonds(
  RDKit::RWMol& molecule,
  const RDKit::Bond* bond,
  const boost::dynamic_bitset<>& bond_mask) {
  boost::dynamic_bitset<> next_bonds (bond_mask.size());
  for (RDKit::Bond* next_bond : molecule.atomBonds(bond->getBeginAtom())) {
    std::size_t next_bond_idx = next_bond->getIdx();
    if (bond_mask[next_bond_idx]) {
      next_bonds.set(next_bond_idx);
    };
  };
  for (RDKit::Bond* next_bond : molecule.atomBonds(bond->getEndAtom())) {
    std::size_t next_bond_idx = next_bond->getIdx();
    if (bond_mask[next_bond_idx]) {
      next_bonds.set(next_bond_idx);
    };
  };
  next_bonds.reset(bond->getIdx());
  return next_bonds;
};

void Kekulize(
  RDKit::RWMol& molecule,
  std::vector<int>& available_valences,
  boost::dynamic_bitset<>& bond_mask) {
  if (bond_mask.none()) {
    return;
  };
  boost::dynamic_bitset<> bonds = bond_mask;
  auto [bond, availability] = LeastAvailableBond(
      molecule, available_valences, bonds);
  // If all bonds have negative availability abort.
  if (!bond) {
    return;
  };
  boost::dynamic_bitset<> next_bonds (bonds.size());
  next_bonds.set(bond->getIdx());
  RDKit::Bond::BondType bond_type = RDKit::Bond::SINGLE; 
  while (!next_bonds.none()) {
    if (!bond) {
      return;
    };
    bond->setBondType(bond_type);
    bond->setIsAromatic(false);
    switch (bond_type) {
      case RDKit::Bond::SINGLE:
        bond_type = RDKit::Bond::DOUBLE;
        break;
      case RDKit::Bond::DOUBLE:
        --available_valences[bond->getBeginAtomIdx()];
        --available_valences[bond->getEndAtomIdx()];
        bond_type = RDKit::Bond::SINGLE;
    };
    bonds.reset(bond->getIdx());
    next_bonds = NextBonds(molecule, bond, bonds);
    if (next_bonds.none() && !bonds.none()) {
      std::tie(bond, availability) = LeastAvailableBond(
        molecule, available_valences, bonds);
      if (bond) {
        next_bonds.set(bond->getIdx());
      };
      bond_type = RDKit::Bond::SINGLE;
      continue;
    };
    std::tie(bond, availability) = LeastAvailableBond(
      molecule, available_valences, next_bonds);
  };
};

void FindSSSRIfNotInitialized(const RDKit::ROMol& molecule) {
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  if (!ring_info->isInitialized()) {
    RDKit::MolOps::findSSSR(molecule);
  };
};

struct RingMask {
  boost::dynamic_bitset<> bonds;
  boost::dynamic_bitset<> aromatic_bonds;
  bool potentially_aromatic = true;

  RingMask(const RDKit::ROMol& molecule, const std::vector<int>& bond_ring) :
    bonds(molecule.getNumBonds()),
    aromatic_bonds(molecule.getNumBonds()),
    potentially_aromatic(bond_ring.size() == 5 || bond_ring.size() == 6) {
    for (std::size_t bond_idx : bond_ring) {
      const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
      const RDKit::Atom* begin_atom = bond->getBeginAtom();
      const RDKit::Atom* end_atom = bond->getEndAtom();
      bonds.set(bond_idx);
      if (bond->getBondType() == RDKit::Bond::AROMATIC) {
        aromatic_bonds.set(bond_idx);
      };
      // Check if the ring contains elements that can't be aromatic.
      if (potentially_aromatic && (
        !open_smiles_aromaticity_specification(begin_atom->getAtomicNum()) ||
        !open_smiles_aromaticity_specification(end_atom->getAtomicNum()))) {
        potentially_aromatic = false;
      };
    };
  };

  double AromaticFraction() const {
    return aromatic_bonds.count() / static_cast<double>(bonds.count());
  };
};

void CorrectAromaticity(RDKit::RWMol& molecule) {
  std::size_t n_atoms = molecule.getNumAtoms(), n_bonds = molecule.getNumBonds();
  boost::dynamic_bitset<> acyclic_bonds_to_kekulize (n_bonds);
  boost::dynamic_bitset<> cyclic_bonds_to_kekulize (n_bonds);
  boost::dynamic_bitset<> aromatized_bonds (n_bonds);
  acyclic_bonds_to_kekulize.set();
  // Calculate the atom's available valences.
  std::vector<int> available_valences (n_atoms, 0);
  for (const RDKit::Atom* atom : molecule.atoms()) {
    int available_valence = AvailableValence(
      atom,
      false, // include_hydrogens
      true   // aromatic_is_single
    );
    // If the available valence is zero all aromatic bonds must be made 
    // single to be valence compliant. If the valence is negative there is 
    // nothing we can do to make the atom valence compliant, but we can at 
    // least not make the issue worse.
    if (available_valence <= 0) {
      for (RDKit::Bond* bond : molecule.atomBonds(atom)) {
        if (bond->getBondType() == RDKit::Bond::AROMATIC) {
          bond->setBondType(RDKit::Bond::SINGLE);
        };
      };
    };
    available_valences[atom->getIdx()] = available_valence;
  };
  // Create bond masks for the rings of the SSSR.
  const RDKit::RingInfo* ring_info = molecule.getRingInfo();
  if (!ring_info->isInitialized()) {
    RDKit::MolOps::findSSSR(molecule);
  };
  std::vector<RingMask> ring_masks;
  ring_masks.reserve(ring_info->numRings());
  for (const std::vector<int>& ring : ring_info->bondRings()) {
    ring_masks.emplace_back(molecule, ring);
  };
  // Sort the rings according to their aromaticity.
  std::sort(ring_masks.begin(), ring_masks.end(), 
    [] (const RingMask& r1, const RingMask& r2) {
      if (r1.potentially_aromatic != r2.potentially_aromatic) {
        return r1.potentially_aromatic > r2.potentially_aromatic;
      };
      return r1.AromaticFraction() > r2.AromaticFraction();
    });
  // Iterate over the rings in order of descending aromaticity.
  for (RingMask& ring_mask : ring_masks) {
    // Add to the aromatic bonds mask bonds that have already been aromatized.
    // This, coupled with the above iteration order, allows aromaticity to 
    // propagate through fused ring systems.
    ring_mask.aromatic_bonds |= ring_mask.bonds & aromatized_bonds;
    acyclic_bonds_to_kekulize -= ring_mask.bonds;
    double aromatic_fraction = ring_mask.AromaticFraction();
    if (ring_mask.potentially_aromatic) {
      // If the ring is fully aromatic leave it be.
      if (aromatic_fraction >= 1.0) {
        aromatized_bonds |= ring_mask.bonds;
        continue;
      // If the ring is partially aromatic, and at least half of the bonds
      // are aromatic, aromatize the ring.
      } else if (aromatic_fraction >= 0.5) {
        Aromatize(molecule, ring_mask.bonds);
        aromatized_bonds |= ring_mask.bonds;
        continue;
      };
    };
    // If the ring can't be aromatic, but less than half of its bonds are
    // labeled as aromatic, we flag the aromatic bonds for kekulization. 
    // We don't kekulize the  ring immediately because there may be other rings 
    // with lower available valences that have priority.
    cyclic_bonds_to_kekulize |= ring_mask.aromatic_bonds;
  };
  // Kekulize the aromatic bonds of non-aromatic cycles.
  // Aromatization has preference over kekulization.
  // We don't want kekulization to affect already aromatized cycles.
  cyclic_bonds_to_kekulize -= aromatized_bonds;
  Kekulize(molecule, available_valences, cyclic_bonds_to_kekulize);
  // Kekulize acyclic aromatic bonds.
  for (std::size_t bond_idx = acyclic_bonds_to_kekulize.find_first();
    bond_idx != boost::dynamic_bitset<>::npos;
    bond_idx = acyclic_bonds_to_kekulize.find_next(bond_idx)) {
    const RDKit::Bond* bond = molecule.getBondWithIdx(bond_idx);
    if (bond->getBondType() != RDKit::Bond::AROMATIC) {
      acyclic_bonds_to_kekulize.reset(bond_idx);
    };
  };
  Kekulize(molecule, available_valences, acyclic_bonds_to_kekulize);
  // Use the new bond types to clean up atom aromaticity labels.
  for (RDKit::Bond* bond : molecule.bonds()) {
    bool is_aromatic = bond->getBondType() == RDKit::Bond::AROMATIC;
    bond->setIsAromatic(is_aromatic);
    bond->getBeginAtom()->setIsAromatic(is_aromatic);
    bond->getEndAtom()->setIsAromatic(is_aromatic);
  };
};

void CorrectHydrogenCounts(RDKit::RWMol& molecule) {
  // Allow the RDKit to recalculate the number of implicit hydrogens. For this 
  // to work properly we must remove all explicit hydrogens and ensure that
  // all atoms are labeled as potentially implicit hydrogen baring.
  for (RDKit::Atom* atom : molecule.atoms()) {
    atom->setNumExplicitHs(0);
    atom->setNoImplicit(false);
  };
  molecule.updatePropertyCache(false);
  // If the number of implicit hydrogens isn't sufficient to satisfy an
  // element's valence we add extra explicit hydrogens. In theory
  // updatePropertyCache should do this for us, but it doesn't always work.
  for (RDKit::Atom* atom : molecule.atoms()) {
    auto [min_n_hydrogens, correctable] = MinNumHydrogensForValidValence(atom);
    if (correctable && min_n_hydrogens > atom->getTotalNumHs()) {
      atom->setNumExplicitHs(min_n_hydrogens);
    };
  };
};

void RemoveBondStereochemistry(RDKit::RWMol& molecule) {
  RDKit::MolOps::clearSingleBondDirFlags(molecule);
  for (RDKit::Bond* bond : molecule.bonds()) {
    bond->setStereo(RDKit::Bond::STEREONONE);
  };
};

void PartialSanitization(
  RDKit::RWMol& molecule) {
  // CorrectAromaticity handles partial aromaticity. This includes promoting 
  // partially aromatic rings to fully aromatic and kekulizing the remainder of 
  // the aromatic bonds. However, it doesn't perceive aromaticity. For example, 
  // C1=CC=CC=C1 is kept kekulized. This task is delegated to the RDKit.
  CorrectAromaticity(molecule);
  CorrectHydrogenCounts(molecule);
  // Remove bond stereochemistry. The RDKit gets REALLY confused if we keep it.
  RemoveBondStereochemistry(molecule);
  // Sanitize the molecule skipping all steps that we have done manually.
  static const unsigned sanitization_flags =
    RDKit::MolOps::SanitizeFlags::SANITIZE_ALL ^
    RDKit::MolOps::SanitizeFlags::SANITIZE_PROPERTIES ^
    RDKit::MolOps::SanitizeFlags::SANITIZE_KEKULIZE ^
    RDKit::MolOps::SanitizeFlags::SANITIZE_FINDRADICALS;
  unsigned failure_flag;
  RDKit::MolOps::sanitizeMol(molecule, failure_flag, sanitization_flags);
};

RDKit::RWMOL_SPTR UnsanitizedMoleculeFromSMILES(const std::string& smiles) {
  // Parse the SMILES without sanitizing the resulting molecule.
  RDKit::RWMOL_SPTR molecule (RDKit::SmilesToMol(smiles, 0, false));
  // SMILES parsing can fail due to a variety of reasons.
  if (!molecule) {
    throw RDKit::SmilesParseException("Couldn't parse SMILES " + smiles);
  };
  return molecule;
};

std::string UnsanitizedMoleculeToSMILES(const RDKit::ROMol& molecule) {
  // RDKit::ROMols store pre-computed properties. These properties (including
  // number of implicit hydrogens, valence, aromaticity, stereochemistry and
  // ring membership) may be invalidated upon molecule edition. However, the
  // properties aren't necessarily erased and may be accessed by other down-
  // stream functions. The user is supposed to sanitize the molecule after
  // modifying it, which recalculates these properties, but doing so in a loop
  // is expensive and finicky since the calculation may fail. Instead we avoid
  // using these properties altogether.
  // NOTE: It seems that sometimes the RDKit somehow still tries to canonicalize 
  // the molecules, despite being instructed not to do so?
  return RDKit::MolToSmiles(
    molecule,
    false,  // Don't include stereochemistry since it's likely invalid.
    false,  // Don't kekulize since it requires aromaticity and valence info.
    -1,
    false); // Don't canonizalize since that requires valid stereochemistry info.
};

#endif // !_MOLECULAR_PERTURBATION_UTILS_HPP_
