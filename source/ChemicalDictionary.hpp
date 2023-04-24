#pragma once
#ifndef _CHEMICAL_DICTIONARY_HPP_
#define _CHEMICAL_DICTIONARY_HPP_

#include "MolecularKeys.hpp"
#include "CircularAtomicEnvironment.hpp"
#include "boost_serialization_tuple.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>

class ChemicalDictionary {
public:
  typedef std::map<AtomKey, std::uint64_t> AtomDictionary;
  typedef std::map<BondKey, std::uint64_t> BondDictionary;

  typedef std::unordered_map<AtomKey::D, std::uint64_t,
    boost::hash<AtomKey::D>> DDictionary;
  typedef std::unordered_map<AtomKey::DV, std::uint64_t,
    boost::hash<AtomKey::DV>> DVDictionary;
  typedef std::unordered_map<AtomKey::DVZ, std::uint64_t,
    boost::hash<AtomKey::DVZ>> DVZDictionary;
  typedef std::unordered_map<AtomKey::DVZQ, std::uint64_t,
    boost::hash<AtomKey::DVZQ>> DVZQDictionary;
  typedef std::unordered_map<BondKey::K1K2, std::uint64_t,
    boost::hash<BondKey::K1K2>> K1K2Dictionary;

  // Vector values are sorted according to frequency in descending order.
  // This allows us to prioritize the most frequent values fulfilling a query.
  typedef std::unordered_map<
    AtomKey::DV, std::vector<std::pair<std::uint8_t, std::uint64_t>>,
    boost::hash<AtomKey::DV>> DV_Z;
  typedef std::unordered_map<
    AtomKey::DVZ, std::vector<std::pair<std::int8_t, std::uint64_t>>,
    boost::hash<AtomKey::DVZ>> DVZ_Q;
  typedef std::unordered_map<
    AtomKey::DVZQ, std::vector<std::pair<std::uint8_t, std::uint64_t>>,
    boost::hash<AtomKey::DVZQ>> DVZQ_H;
  typedef std::unordered_map<
    BondKey::K1K2, std::vector<std::pair<RDKit::Bond::BondType, std::uint64_t>>,
    boost::hash<BondKey::K1K2>> K1K2_B;

  // EnvironmentDictionary stores only CircularAtomicEnvironment hashes as
  // opposed to the actual object. This is because CircularAtomicEnvironment
  // references atoms and bonds of a molecule, and the moment said molecule
  // is destructed these references become invalid. Storing molecule copies
  // would be expensive. The only advantage would be to test for equality with
  // operator==, but since this would also be expensive and the probability of
  // a hash collision is small we skip this.
  typedef std::unordered_map<EnvironmentKey, std::uint64_t> EnvironmentDictionary;

private:
  AtomDictionary atom_dictionary;
  BondDictionary bond_dictionary;

  DDictionary d_dictionary;
  DVDictionary dv_dictionary;
  DVZDictionary dvz_dictionary;
  DVZQDictionary dvzq_dictionary;
  K1K2Dictionary k1k2_dictionary;

  DV_Z dv_z;
  DVZ_Q dvz_q;
  DVZQ_H dvzq_h;
  K1K2_B k1k2_b;

  EnvironmentDictionary environment_dictionary;
  CircularAtomicEnvironmentGenerator environment_generator;
  unsigned environment_radius = 2;

  std::uint64_t total_atom_frequency = 0, total_bond_frequency = 0;

public:
  std::uint64_t foreign_d_frequency_threshold = 1;
  std::uint64_t foreign_dv_frequency_threshold = 1;
  std::uint64_t foreign_dvz_frequency_threshold = 1;
  std::uint64_t foreign_dvzq_frequency_threshold = 1;
  std::uint64_t foreign_k1k2_frequency_threshold = 1;
  std::uint64_t foreign_atom_frequency_threshold = 1;
  std::uint64_t foreign_bond_frequency_threshold = 1;
  std::uint64_t foreign_environment_frequency_threshold = 1;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& archive, const unsigned version = 20230424) {
    // NOTE: We can't serialize atoms_hasher because it's a functional.
    // If you are not using the default hasher you should ensure you construct
    // deserialized instances with the same hasher of the serialized instance.
    archive &
      atom_dictionary & bond_dictionary & d_dictionary & dv_dictionary &
      dvz_dictionary & dvzq_dictionary & k1k2_dictionary &
      dv_z & dvz_q & dvzq_h & k1k2_b &
      environment_dictionary & environment_radius &
      total_atom_frequency & total_bond_frequency &
      foreign_d_frequency_threshold & foreign_dv_frequency_threshold &
      foreign_dvz_frequency_threshold & foreign_dvzq_frequency_threshold &
      foreign_k1k2_frequency_threshold & foreign_atom_frequency_threshold &
      foreign_bond_frequency_threshold & foreign_environment_frequency_threshold;
  };

  template <class First, class Second>
  std::vector<First> FirstOfPairsVector(
    const std::vector<std::pair<First, Second>>& pairs_vector) const {
    std::size_t n = pairs_vector.size();
    std::vector<First> firsts (n);
    for (std::size_t i = 0; i < n; ++i) {
      firsts[i] = pairs_vector[i].first;
    };
    return firsts;
  };

  template <class Dictionary, class Key>
  void IncrementDictionaryValue(
    Dictionary& dictionary,const Key& key) {
    auto [it, emplaced] = dictionary.emplace(key, 1);
    if (!emplaced) {
      ++(it->second);
    };
  };

public:
  ChemicalDictionary(unsigned environment_radius = 2) :
    environment_radius(environment_radius),
    environment_generator(environment_radius) {};
  ChemicalDictionary(
    const std::string& path) {
    Load(path);
    environment_generator = CircularAtomicEnvironmentGenerator(
      environment_radius);
  };

  void AddMolecule(const RDKit::ROMol& molecule) {
    MolecularKeys molecular_keys (molecule);
    total_atom_frequency += molecular_keys.atom_keys.size();
    total_bond_frequency += molecular_keys.bond_keys.size();
    for (const AtomKey& atom_key : molecular_keys.atom_keys) {
      IncrementDictionaryValue(d_dictionary, atom_key.d());
      IncrementDictionaryValue(dv_dictionary, atom_key.dv());
      IncrementDictionaryValue(dvz_dictionary, atom_key.dvz());
      IncrementDictionaryValue(dvzq_dictionary, atom_key.dvzq());
      IncrementDictionaryValue(atom_dictionary, std::move(atom_key));
    };
    for (const BondKey& bond_key : molecular_keys.bond_keys) {
      IncrementDictionaryValue(k1k2_dictionary, bond_key.k1k2());
      IncrementDictionaryValue(bond_dictionary, std::move(bond_key));
    };
    for (const RDKit::Atom* atom : molecule.atoms()) {
      IncrementDictionaryValue(
        environment_dictionary, environment_generator.Key(atom));
    };
  };

  void BuildPartialAtomKeyDictionaries() {
    auto second_gt = [] (const auto& p1, const auto& p2) {
      return p1.second > p2.second;
    };

    dv_z.clear();
    dvz_q.clear();
    dvzq_h.clear();

    AtomDictionary::const_iterator
      begin_it = atom_dictionary.cbegin(),
      end_it = atom_dictionary.cend();
    const auto& [first_atom_key, first_atom_frequency] = *begin_it;
    ++begin_it;

    AtomKey::DV prev_dv = first_atom_key.dv();
    AtomKey::DVZ prev_dvz = first_atom_key.dvz();
    AtomKey::DVZQ prev_dvzq = first_atom_key.dvzq();

    std::vector<std::pair<std::uint8_t, std::uint64_t>> dv_zf;
    std::vector<std::pair<std::int8_t, std::uint64_t>> dvz_qf;
    std::vector<std::pair<std::uint8_t, std::uint64_t>> dvzq_hf;

    dv_zf.emplace_back(first_atom_key.atomic_number, first_atom_frequency);
    dvz_qf.emplace_back(first_atom_key.formal_charge, first_atom_frequency);
    dvzq_hf.emplace_back(first_atom_key.n_hydrogens, first_atom_frequency);

    for (auto it = begin_it; it != end_it; ++it) {
      const auto& [atom_key, frequency] = *it;

      AtomKey::DV dv = atom_key.dv();
      AtomKey::DVZ dvz = atom_key.dvz();
      AtomKey::DVZQ dvzq = atom_key.dvzq();

      if (dv == prev_dv) {
        auto& [z, f] = dv_zf.back();
        if (atom_key.atomic_number == z) {
          f += frequency;
        } else {
          dv_zf.emplace_back(atom_key.atomic_number, frequency);
        };
      } else {
        std::sort(dv_zf.begin(), dv_zf.end(), second_gt);
        dv_z.emplace(prev_dv, std::move(dv_zf));
        dv_zf.clear();
        dv_zf.emplace_back(atom_key.atomic_number, frequency);
      };

      if (dvz == prev_dvz) {
        auto& [q, f] = dvz_qf.back();
        if (atom_key.formal_charge == q) {
          f += frequency;
        } else {
          dvz_qf.emplace_back(atom_key.formal_charge, frequency);
        };
      } else {
        std::sort(dvz_qf.begin(), dvz_qf.end(), second_gt);
        dvz_q.emplace(prev_dvz, std::move(dvz_qf));
        dvz_qf.clear();
        dvz_qf.emplace_back(atom_key.formal_charge, frequency);
      };

      if (dvzq == prev_dvzq) {
        dvzq_hf.emplace_back(atom_key.n_hydrogens, frequency);
      } else {
        std::sort(dvzq_hf.begin(), dvzq_hf.end(), second_gt);
        dvzq_h.emplace(prev_dvzq, std::move(dvzq_hf));
        dvzq_hf.clear();
        dvzq_hf.emplace_back(atom_key.n_hydrogens, frequency);
      };

      prev_dv = dv;
      prev_dvz = dvz;
      prev_dvzq = dvzq;
    };

    std::sort(dv_zf.begin(), dv_zf.end(), second_gt);
    std::sort(dvz_qf.begin(), dvz_qf.end(), second_gt);
    std::sort(dvzq_hf.begin(), dvzq_hf.end(), second_gt);
    dv_z.emplace(prev_dv, std::move(dv_zf));
    dvz_q.emplace(prev_dvz, std::move(dvz_qf));
    dvzq_h.emplace(prev_dvzq, std::move(dvzq_hf));
  };

  void BuildPartialBondKeyDictionary() {
    auto second_gt = [] (const auto& p1, const auto& p2) {
      return p1.second > p2.second;
    };

    k1k2_b.clear();

    BondDictionary::const_iterator
      begin_it = bond_dictionary.cbegin(),
      end_it = bond_dictionary.cend();
    const auto& [first_bond_key, first_bond_frequency] = *begin_it;
    ++begin_it;

    BondKey::K1K2 prev_k1k2 = first_bond_key.k1k2();

    std::vector<std::pair<RDKit::Bond::BondType, std::uint64_t>> k1k2_bf;
    k1k2_bf.emplace_back(first_bond_key.bond_type, first_bond_frequency);

    for (auto it = begin_it; it != end_it; ++it) {
      const auto& [bond_key, frequency] = *it;
      BondKey::K1K2 k1k2 = bond_key.k1k2();
      if (k1k2 == prev_k1k2) {
        k1k2_bf.emplace_back(bond_key.bond_type, frequency);
      } else {
        std::sort(k1k2_bf.begin(), k1k2_bf.end(), second_gt);
        k1k2_b.emplace(prev_k1k2, std::move(k1k2_bf));
        k1k2_bf.clear();
        k1k2_bf.emplace_back(bond_key.bond_type, frequency);
      };
      prev_k1k2 = k1k2;
    };

    std::sort(k1k2_bf.begin(), k1k2_bf.end(), second_gt);
    k1k2_b.emplace(prev_k1k2, std::move(k1k2_bf));
  };

  void BuildPartialKeyDictionaries() {
    BuildPartialAtomKeyDictionaries();
    BuildPartialBondKeyDictionary();
  };

  std::uint64_t AtomFrequency(const AtomKey& atom_key) const {
    AtomDictionary::const_iterator it = atom_dictionary.find(atom_key);
    return it == atom_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t BondFrequency(const BondKey& bond_key) const {
    BondDictionary::const_iterator it = bond_dictionary.find(bond_key);
    return it == bond_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t DFrequency(const AtomKey::D& d) const {
    DDictionary::const_iterator it = d_dictionary.find(d);
    return it == d_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t DVFrequency(const AtomKey::DV& dv) const {
    DVDictionary::const_iterator it = dv_dictionary.find(dv);
    return it == dv_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t DVZFrequency(const AtomKey::DVZ& dvz) const {
    DVZDictionary::const_iterator it = dvz_dictionary.find(dvz);
    return it == dvz_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t DVZQFrequency(const AtomKey::DVZQ& dvzq) const {
    DVZQDictionary::const_iterator it = dvzq_dictionary.find(dvzq);
    return it == dvzq_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t K1K2Frequency(const BondKey::K1K2& k1k2) const {
    K1K2Dictionary::const_iterator it = k1k2_dictionary.find(k1k2);
    return it == k1k2_dictionary.cend() ? 0 : it->second;
  };

  std::uint64_t EnvironmentFrequency(
    const EnvironmentKey& environment_key) const {
    EnvironmentDictionary::const_iterator
      it = environment_dictionary.find(environment_key);
    return it == environment_dictionary.cend() ? 0 : it->second;
  };

  std::vector<std::uint8_t> Z_DV(const AtomKey::DV& dv) const {
    return FirstOfPairsVector(dv_z.at(dv));
  };

  std::vector<std::int8_t> Q_DVZ(const AtomKey::DVZ& dvz) const {
    return FirstOfPairsVector(dvz_q.at(dvz));
  };

  std::vector<std::uint8_t> H_DVZQ(const AtomKey::DVZQ& dvzq) const {
    return FirstOfPairsVector(dvzq_h.at(dvzq));
  };

  std::vector<RDKit::Bond::BondType> B_K1K2(const BondKey::K1K2& k1k2) const {
    return FirstOfPairsVector(k1k2_b.at(k1k2));
  };

  bool IsForeignAtom(const AtomKey& atom_key) const {
    return AtomFrequency(atom_key) < foreign_atom_frequency_threshold;
  };

  bool IsForeignAtom(const RDKit::Atom* atom) const {
    return IsForeignAtom(AtomKey(atom));
  };

  bool IsForeignBond(const BondKey& bond_key) const {
    return BondFrequency(bond_key) < foreign_bond_frequency_threshold;
  };

  bool IsForeignBond(const RDKit::Bond* bond) const {
    return IsForeignBond(BondKey(bond));
  };

  CircularAtomicEnvironment Environment(const RDKit::Atom* atom) const {
    return environment_generator(atom);
  };

  bool IsForeignEnvironment(const EnvironmentKey& environment_key) const {
    return EnvironmentFrequency(environment_key) <
      foreign_environment_frequency_threshold;
  };

  bool IsForeignEnvironment(const RDKit::Atom* atom) {
    return IsForeignEnvironment(environment_generator.Key(atom));
  };

  AtomKey::Error AtomKeyError(const AtomKey& atom_key) const {
    if (!IsForeignAtom(atom_key)) {
      return AtomKey::Error::NONE;
    };
    if (DFrequency(atom_key.d()) < foreign_d_frequency_threshold) {
      return AtomKey::Error::D;
    };
    if (DVFrequency(atom_key.dv()) < foreign_dv_frequency_threshold) {
      return AtomKey::Error::DV;
    };
    if (DVZFrequency(atom_key.dvz()) < foreign_dvz_frequency_threshold) {
      return AtomKey::Error::DVZ;
    };
    if (DVZQFrequency(atom_key.dvzq()) < foreign_dvzq_frequency_threshold) {
      return AtomKey::Error::DVZQ;
    };
    return AtomKey::Error::DVZQH;
  };

  BondKey::Error BondKeyError(const BondKey& bond_key) const {
    if (!IsForeignBond(bond_key)) {
      return BondKey::Error::NONE;
    };
    if (K1K2Frequency(bond_key.k1k2()) < foreign_k1k2_frequency_threshold) {
      return BondKey::Error::K1K2;
    };
    return BondKey::Error::K1K2B;
  };

  std::uint64_t GetTotalAtomFrequency() const {
    return total_atom_frequency;
  };

  std::uint64_t GetTotalBondFrequency() const {
    return total_bond_frequency;
  };

  const AtomDictionary& GetAtomDictionary() const {
    return atom_dictionary;
  };

  const BondDictionary& GetBondDictionary() const {
    return bond_dictionary;
  };

  const EnvironmentDictionary& GetEnvironmentDictionary() const {
    return environment_dictionary;
  };

  const CircularAtomicEnvironmentGenerator& GetEnvironmentGenerator() const {
    return environment_generator;
  };

  void Save(const std::string& path) const {
    std::ofstream output_stream (path, std::ofstream::binary);
    boost::archive::binary_oarchive archive (output_stream);
    archive << *this;
    output_stream.close();
  };

  void Load(const std::string& path) {
    std::ifstream input_stream (path, std::ifstream::binary);
    boost::archive::binary_iarchive archive (input_stream);
    archive >> *this;
    input_stream.close();
  };
};

#endif // !_CHEMICAL_DICTIONARY_HPP_
