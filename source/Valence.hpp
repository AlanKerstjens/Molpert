#ifndef _VALENCE_HPP_
#define _VALENCE_HPP_

#include <GraphMol/ROMol.h>
#include <GraphMol/PeriodicTable.h>

double ExplicitValence(
  const RDKit::Atom* atom,
  bool include_hydrogens = true,
  bool aromatic_is_single = false) {
  // The way in which the RDKit calculates valence is weird and causes many
  // headaches when working with edited molecules. We define valence as the
  // sum of bond orders. Note that we don't include formal charge in the
  // calculation to avoid cyclic dependencies between atom properties. We
  // define formal charge as a less significant property and don't want it to
  // influence the more significant property of valence.
  const RDKit::ROMol& molecule = atom->getOwningMol();
  double valence = include_hydrogens ? atom->getNumExplicitHs() : 0.0;
  for (const RDKit::Bond* bond : molecule.atomBonds(atom)) {
    if (aromatic_is_single && bond->getBondType() == RDKit::Bond::AROMATIC) {
      valence += 1.0;
      continue;
    };
    valence += bond->getValenceContrib(atom);
  };
  // We return a double, which may be truncated to an int. This would be
  // undesirable behaviour if our double approached an integer from below (e.g.
  // 5.999999), but assuming a 64-bit IEEE 754 floating point representation we
  // should be able to perfectly represent integers up to 53-bits. Given that
  // valence values tend to be small this should be plenty to avoid the problem.
  return valence;
};

int IntValenceDiff(double valence_diff) {
  // When the valence of an atom increases or decreases by a float, truncation
  // may hide the real change in valence. For example, a change of +0.5 would be
  // truncated to 0. We solve this by rounding away from zero.
  // Adding an aromatic bond changes the valence by +1.5, rounded to +2.
  // Removing an aromatic bond changes the valence by -1.5, rounded to -2.
  // As an additional benefit this cancels out the error due to calculating base
  // valences through truncation. For example, suppose an atom had a valence of
  // +1.5, truncated to +1. Adding an aromatic bond would bring the valence to,
  // +3.0, which is correctly reflected by adding the ceiling of +1.5 (+2) to +1.
  return valence_diff < 0 ? std::floor(valence_diff) : std::ceil(valence_diff);
};

double BondTypeValenceContribution(RDKit::Bond::BondType bond_type) {
  // Mimics the behaviour of RDKit::Bond::getValenceContrib, but can be called
  // without RDKit::Atom and RDKit::Bond instances.
  switch (bond_type) {
    case RDKit::Bond::UNSPECIFIED:
    case RDKit::Bond::IONIC:
    case RDKit::Bond::HYDROGEN:
    case RDKit::Bond::ZERO:
      return 0.0;
    case RDKit::Bond::SINGLE:
      return 1.0;
    case RDKit::Bond::DOUBLE:
      return 2.0;
    case RDKit::Bond::TRIPLE:
      return 3.0;
    case RDKit::Bond::QUADRUPLE:
      return 4.0;
    case RDKit::Bond::QUINTUPLE:
      return 5.0;
    case RDKit::Bond::HEXTUPLE:
      return 6.0;
    case RDKit::Bond::AROMATIC:
    case RDKit::Bond::ONEANDAHALF:
      return 1.5;
    case RDKit::Bond::TWOANDAHALF:
      return 2.5;
    case RDKit::Bond::THREEANDAHALF:
      return 3.5;
    case RDKit::Bond::FOURANDAHALF:
      return 4.5;
    case RDKit::Bond::FIVEANDAHALF:
      return 5.5;
    default:
      throw std::runtime_error(
        "Unknown valence contribution for bond type " + bond_type);
  };
};

int MaxAllowedValence(unsigned atomic_number) {
  static const RDKit::PeriodicTable* periodic_table =
    RDKit::PeriodicTable::getTable();
  const std::vector<int>& allowed_valences =
    periodic_table->getValenceList(atomic_number);
  return allowed_valences.back();
};

int AvailableValence(
  const RDKit::Atom* atom,
  bool include_hydrogens = true,
  bool aromatic_is_single = false) {
  int max_allowed_valence = MaxAllowedValence(atom->getAtomicNum());
  if (max_allowed_valence < 0) {
    return std::numeric_limits<int>::max();
  };
  return max_allowed_valence - ExplicitValence(
    atom, include_hydrogens, aromatic_is_single);
};

bool IsValenceBelowMax(unsigned atomic_number, int valence) {
  if (valence < 0) {
    return false;
  };
  if (valence == 0) {
    return true;
  };
  // Each atomic number is assigned a maximum valence. So long as the valence
  // falls below this value it's assumed that implicit hydrogens can be added to
  // reach some valid valence, which may be the maximum or some lower value.
  // Note that you may still be responsible for determining the appropiate
  // number of hydrogens to add (usually the lowest number you get away with).
  // As an example, the maximum valid valence of sulphur is 6, yet a single
  // explicit sulfur atom should be padded with hydrogens to become [SH2]
  // (valence = 2) as opposed to [SH6] (valence = 6).
  int max_allowed_valence = MaxAllowedValence(atomic_number);
  // In the RDKit::PeriodicTable elements for which no valence information is
  // available are assigned a valence list of {-1}. It should be interpreted as
  // any (positive) valence being valid.
  if (max_allowed_valence < 0) {
    return true;
  };
  return valence <= max_allowed_valence;
};

// Calculate the minimum number of explicit hydrogens necessary for the atom's
// valence to be valid.
std::pair<unsigned, bool> MinNumHydrogensForValidValence(
  const RDKit::Atom* atom) {
  static const RDKit::PeriodicTable* periodic_table =
    RDKit::PeriodicTable::getTable();
  const std::vector<int>& allowed_valences =
    periodic_table->getValenceList(atom->getAtomicNum());
  int max_allowed_valence = allowed_valences.back();
  if (max_allowed_valence < 0) {
    return {0, true};
  };
  int valence = ExplicitValence(atom, false);
  if (valence > max_allowed_valence) {
    return {0, false};
  };
  for (int allowed_valence : allowed_valences) {
    if (valence <= allowed_valence) {
      return {allowed_valence - valence, true};
    };
  };
  return {0, false}; // We should never get here. Just to shut up the compiler.
};

#endif // !_VALENCE_HPP_
