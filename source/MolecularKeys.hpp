#pragma once
#ifndef _MOLECULAR_KEYS_HPP_
#define _MOLECULAR_KEYS_HPP_

#include "Valence.hpp"

struct AtomKey {
  // The order of the key variables isn't capricious. They are listed according
  // to the variable's permanence, in descending order. Degree (D) and valence
  // (V) depend on the surrounding chemical environment and can't be changed by
  // modifying the subject atom. Atomic number (Z), formal charge (Q) and number
  // of hydrogens (H) can be modified directly by accessing the corresponding
  // RDKit::Atom property, but not all combinations of these three properties
  // are allowed. As an example, the allowed formal charges depend on the atomic
  // number and vice versa. However, it's of interest to avoid cyclic dependencies
  // between properties as they could lead to algorithm soft-lock. Instead, we
  // consider that number of hydrogens depends on formal charge (DVZQ), which in
  // turn depends on atomic number (DVZ), which in turn depends on the degree and
  // valence (DV). Inevitably this means that when changing one of the more
  // permanent properties (e.g. atomic number) we may invalidate the accompanying
  // less permanent properties (e.g. formal charge), but we can easily fix this
  // down the line.

  typedef std::uint8_t D;
  typedef std::pair <std::uint8_t, std::int8_t> DV;
  typedef std::tuple<std::uint8_t, std::int8_t, std::uint8_t> DVZ;
  typedef std::tuple<std::uint8_t, std::int8_t, std::uint8_t, std::int8_t> DVZQ;

  // Ordered according to severity.
  enum class Error {
    NONE = 0,
    DVZQH = 1,
    DVZQ = 2,
    DVZ = 3,
    DV = 4,
    D = 5
  };

  std::uint8_t degree = 0;        // D (explicit)
  std::int8_t valence = 0;        // V (explicit)
  std::uint8_t atomic_number = 0; // Z
  std::int8_t formal_charge = 0;  // Q
  std::uint8_t n_hydrogens = 0;   // H (explicit)

  AtomKey() = default;

  AtomKey(
    std::uint8_t degree,
    std::int8_t valence,
    std::uint8_t atomic_number,
    std::int8_t formal_charge,
    std::uint8_t n_hydrogens) :
    degree(degree),
    valence(valence),
    atomic_number(atomic_number),
    formal_charge(formal_charge),
    n_hydrogens(n_hydrogens) {};

  AtomKey(
    const RDKit::Atom* atom) :
    degree(atom->getDegree()),
    valence(ExplicitValence(atom)),
    atomic_number(atom->getAtomicNum()),
    formal_charge(atom->getFormalCharge()),
    n_hydrogens(atom->getNumExplicitHs()) {};

  AtomKey(
    const RDKit::Atom* atom,
    double valence) :
    degree(atom->getDegree()),
    valence(valence),
    atomic_number(atom->getAtomicNum()),
    formal_charge(atom->getFormalCharge()),
    n_hydrogens(atom->getNumExplicitHs()) {};

  friend auto operator<=>(const AtomKey&, const AtomKey&) = default;
  friend std::hash<AtomKey>;
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& archive, const unsigned version = 20221125) {
    archive & degree & valence & atomic_number & formal_charge & n_hydrogens;
  };

  D d() const {
    return degree;
  };

  DV dv() const {
    return {degree, valence};
  };

  DVZ dvz() const {
    return {degree, valence, atomic_number};
  };

  DVZQ dvzq() const {
    return {degree, valence, atomic_number, formal_charge};
  };

  std::string str() const {
    std::stringstream stream;
    stream << "("
           << static_cast<unsigned>(degree) << ","
           << static_cast<int>(valence) << ","
           << static_cast<unsigned>(atomic_number) << ","
           << static_cast<int>(formal_charge) << ","
           << static_cast<unsigned>(n_hydrogens) << ")";
    return stream.str();
  };
};

// std::hash specialization.
template <>
struct std::hash<AtomKey> {
  std::size_t operator()(const AtomKey& atom_key) const {
    std::size_t hash = atom_key.degree;
    boost::hash_combine(hash, atom_key.valence);
    boost::hash_combine(hash, atom_key.atomic_number);
    boost::hash_combine(hash, atom_key.formal_charge);
    boost::hash_combine(hash, atom_key.n_hydrogens);
    return hash;
  };
};

// boost::hash_value specialization. We need this to use boost::hash.
std::size_t hash_value(const AtomKey& atom_key) {
  static std::hash<AtomKey> atom_key_hasher;
  return atom_key_hasher(atom_key);
};

// operator<< overloads for AtomKey partial keys
std::ostream& operator<<(std::ostream& os, const AtomKey::DV& dv) {
  os << "(" << (unsigned) dv.first << ","
            << (int) dv.second << ")";
  return os;
};

std::ostream& operator<<(std::ostream& os, const AtomKey::DVZ& dvz) {
  os << "(" << (unsigned) std::get<0>(dvz) << ","
            << (int) std::get<1>(dvz) << ","
            << (unsigned) std::get<2>(dvz) << ")";
  return os;
};

std::ostream& operator<<(std::ostream& os, const AtomKey::DVZQ& dvzq) {
  os << "(" << (unsigned) std::get<0>(dvzq) << ","
            << (int) std::get<1>(dvzq) << ","
            << (unsigned) std::get<2>(dvzq) << ","
            << (int) std::get<3>(dvzq) << ")";
  return os;
};

static const AtomKey NULL_ATOM_KEY (0, 0, 0, 0, 0);


struct BondKey {
  typedef std::pair<AtomKey, AtomKey> K1K2;

  enum class Error {
    NONE = 0,
    K1K2B = 1,
    K1K2 = 2
  };

  AtomKey atom_key1;                                          // K1
  AtomKey atom_key2;                                          // K2
  RDKit::Bond::BondType bond_type = RDKit::Bond::UNSPECIFIED; // B

  BondKey() = default;

  BondKey(
    const AtomKey& ak1,
    const AtomKey& ak2,
    RDKit::Bond::BondType bond_type,
    bool canonicalize = true) :
    atom_key1(ak1),
    atom_key2(ak2),
    bond_type(bond_type) {
    // Sort the atom keys for canonicalization purposes.
    if (canonicalize) {
      if (atom_key2 < atom_key1) {
        std::swap(atom_key1, atom_key2);
      };
    };
  };

  BondKey(
    const RDKit::Atom* atom1,
    const RDKit::Atom* atom2,
    RDKit::Bond::BondType bond_type,
    bool canonicalize = true) :
    atom_key1(atom1),
    atom_key2(atom2),
    bond_type(bond_type) {
    if (canonicalize) {
      if (atom_key2 < atom_key1) {
        std::swap(atom_key1, atom_key2);
      };
    };
  };

  BondKey(
    const RDKit::Bond* bond,
    bool canonicalize = true) :
    atom_key1(bond->getBeginAtom()),
    atom_key2(bond->getEndAtom()),
    bond_type(bond->getBondType()) {
    if (canonicalize) {
      if (atom_key2 < atom_key1) {
        std::swap(atom_key1, atom_key2);
      };
    };
  };

  friend auto operator<=>(const BondKey&, const BondKey&) = default;
  friend std::hash<BondKey>;
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& archive, const unsigned version = 20220902) {
    archive & atom_key1 & atom_key2 & bond_type;
  };

  K1K2 k1k2() const {
    return {atom_key1, atom_key2};
  };

  std::string str() const {
    std::stringstream stream;
    stream << atom_key1.str()
           << "-[" << bond_type << "]->"
           << atom_key2.str();
    return stream.str();
  };
};

template <>
struct std::hash<BondKey> {
  std::hash<AtomKey> atom_key_hasher;
  std::size_t operator()(const BondKey& bond_key) const {
    std::size_t hash = atom_key_hasher(bond_key.atom_key1);
    boost::hash_combine(hash, atom_key_hasher(bond_key.atom_key2));
    boost::hash_combine(hash, bond_key.bond_type);
    return hash;
  };
};

std::size_t hash_value(const BondKey& bond_key) {
  static std::hash<BondKey> bond_key_hasher;
  return bond_key_hasher(bond_key);
};

std::ostream& operator<<(std::ostream& os, const BondKey::K1K2& k1k2) {
  os << k1k2.first.str() << "->"
     << k1k2.second.str();
  return os;
};

static const BondKey NULL_BOND_KEY (
  NULL_ATOM_KEY, NULL_ATOM_KEY, RDKit::Bond::UNSPECIFIED);


struct MolecularKeys {
  std::vector<AtomKey> atom_keys;
  std::vector<BondKey> bond_keys;

  MolecularKeys() = default;
  MolecularKeys(const RDKit::ROMol& molecule) {
    atom_keys.reserve(molecule.getNumAtoms());
    bond_keys.reserve(molecule.getNumBonds());
    for (const RDKit::Atom* atom : molecule.atoms()) {
      atom_keys.emplace_back(atom);
    };
    for (const RDKit::Bond* bond : molecule.bonds()) {
      bond_keys.emplace_back(
        atom_keys[bond->getBeginAtomIdx()],
        atom_keys[bond->getEndAtomIdx()],
        bond->getBondType());
    };
  };

  std::size_t size() const {
    return atom_keys.size() + bond_keys.size();
  };
};


typedef std::pair<AtomKey, AtomKey> AtomKeyChange;
typedef std::pair<BondKey, BondKey> BondKeyChange;

#endif // !_MOLECULAR_KEYS_HPP_
