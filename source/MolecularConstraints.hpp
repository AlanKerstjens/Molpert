#pragma once
#ifndef _MOLECULAR_CONSTRAINTS_HPP_
#define _MOLECULAR_CONSTRAINTS_HPP_

#include "MolecularPerturbations.hpp"
#include "ChemicalDictionary.hpp"
#include <boost/core/null_deleter.hpp>

class MolecularConstraints {
public:
  // Given a molecular key, a constraint evaluates if the key is allowed (true)
  // or not (false).
  typedef std::function<bool(const AtomKey&)> AtomConstraint;
  typedef std::function<bool(const BondKey&)> BondConstraint;
  typedef std::function<bool(const EnvironmentKey&)> EnvironmentConstraint;
  // Built-in constraint types.
  enum class AtomConstraintType {Null, Valence, AtomKey};
  enum class BondConstraintType {Null, BondKey};
  enum class EnvironmentConstraintType {Null, EnvironmentKey};
  // Constraint generators compare two molecular keys (prior and posterior) and,
  // if pertinent, generate new constraints for the posterior key. The return
  // type is a std::optional that should be empty if the constraint hasn't
  // changed. If the constraint has been erased it should be set to nullptr.
  typedef std::function<
    std::optional<AtomConstraint>(const AtomKeyChange&)
  > AtomConstraintGenerator;
  typedef std::function<
    std::optional<BondConstraint>(const BondKeyChange&)
  > BondConstraintGenerator;
  typedef std::function<
    std::optional<EnvironmentConstraint>(const EnvironmentKeyChange&)
  > EnvironmentConstraintGenerator;

private:
  std::map<Tag, AtomConstraint> atom_constraints;
  std::map<Tag, BondConstraint> bond_constraints;
  std::map<Tag, EnvironmentConstraint> environment_constraints;
  AtomConstraintGenerator atom_constraint_generator;
  BondConstraintGenerator bond_constraint_generator;
  EnvironmentConstraintGenerator environment_constraint_generator;
  CircularAtomicEnvironmentGenerator environment_generator;
  unsigned min_cycle_size = 3;
  unsigned max_cycle_size = std::numeric_limits<unsigned>::max();

private:
  template <class Constraint>
  bool SetConstraint(
    Tag tag,
    const Constraint& constraint,
    std::map<Tag, Constraint>& constraints,
    bool replace = true) {
    // If we are provided a null constraint we either clear the existing one or
    // do nothing at all. This ensures we never store null constraints.
    if (!constraint) {
      if (replace) {
        return constraints.erase(tag);
      };
      return false;
    };
    // Try to insert a new constraint.
    auto [it, emplaced] = constraints.emplace(tag, constraint);
    if (emplaced) {
      return true;
    };
    if (!replace) {
      return false;
    };
    // If we already had a constraint for the tag, if possible, replace it.
    it->second = constraint;
    return true;
  };

  template <class Constraint>
  const Constraint* GetConstraint(
    Tag tag,
    const std::map<Tag, Constraint>& constraints) const {
    typename std::map<Tag, Constraint>::const_iterator it = constraints.find(tag);
    // Return a null constraint if the tag has no constraints.
    return it != constraints.cend() ? &it->second : nullptr;
  };

  std::pair<std::shared_ptr<const AtomConstraint>, bool> UpdatedAtomConstraint(
    Tag atom_tag,
    const AtomKeyChange& atom_key_change) const {
    std::shared_ptr<const AtomConstraint> updated_atom_constraint;
    // If no constraint generator was defined we can't update anything.
    // Return the current constraint.
    if (!atom_constraint_generator) {
      updated_atom_constraint.reset(
        GetConstraint(atom_tag, atom_constraints), boost::null_deleter());
      return {updated_atom_constraint, false};
    };
    // Try to generate a new constraint.
    std::optional<AtomConstraint> new_atom_constraint =
      atom_constraint_generator(atom_key_change);
    // If no new constraint was generated we don't update anything.
    // Return the current constraint.
    if (!new_atom_constraint) {
      updated_atom_constraint.reset(
        GetConstraint(atom_tag, atom_constraints), boost::null_deleter());
      return {updated_atom_constraint, false};
    };
    // If a non-null constraint was generated take ownership of it.
    if (*new_atom_constraint) {
      updated_atom_constraint.reset(
        new AtomConstraint(std::move(*new_atom_constraint)));
    };
    // Return the new constraint.
    return {updated_atom_constraint, true};
  };

  std::pair<std::shared_ptr<const BondConstraint>, bool> UpdatedBondConstraint(
    Tag bond_tag,
    const BondKeyChange& bond_key_change) const {
    std::shared_ptr<const BondConstraint> updated_bond_constraint;
    if (!bond_constraint_generator) {
      updated_bond_constraint.reset(
        GetConstraint(bond_tag, bond_constraints), boost::null_deleter());
      return {updated_bond_constraint, false};
    };
    std::optional<BondConstraint> new_bond_constraint =
      bond_constraint_generator(bond_key_change);
    if (!new_bond_constraint) {
      updated_bond_constraint.reset(
        GetConstraint(bond_tag, bond_constraints), boost::null_deleter());
      return {updated_bond_constraint, false};
    };
    if (*new_bond_constraint) {
      updated_bond_constraint.reset(
        new BondConstraint(std::move(*new_bond_constraint)));
    };
    return {updated_bond_constraint, true};
  };

  std::pair<std::shared_ptr<const EnvironmentConstraint>, bool> 
  UpdatedEnvironmentConstraint(
    Tag atom_tag,
    const EnvironmentKeyChange& environment_key_change) const {
    std::shared_ptr<const EnvironmentConstraint> updated_environment_constraint;
    if (!environment_constraint_generator) {
      updated_environment_constraint.reset(
        GetConstraint(atom_tag, environment_constraints), boost::null_deleter());
      return {updated_environment_constraint, false};
    };
    std::optional<EnvironmentConstraint> new_environment_constraint =
      environment_constraint_generator(environment_key_change);
    if (!new_environment_constraint) {
      updated_environment_constraint.reset(
        GetConstraint(atom_tag, environment_constraints), boost::null_deleter());
      return {updated_environment_constraint, false};
    };
    if (*new_environment_constraint) {
      updated_environment_constraint.reset(
        new EnvironmentConstraint(std::move(*new_environment_constraint)));
    };
    return {updated_environment_constraint, true};
  };

  bool CyclicityConstraintsSatisfied(
    const MolecularGraphProjection& projection) const {
    if (HasCyclicityConstraints()) {
      std::vector<boost::dynamic_bitset<>> mcb = projection.MinimumCycleBasis();
      for (const boost::dynamic_bitset<>& cycle : mcb) {
        std::size_t cycle_size = cycle.count();
        if (cycle_size < min_cycle_size || cycle_size > max_cycle_size) {
          return false;
        };
      };
    };
    return true;
  };

public:
  MolecularConstraints() = default;
  MolecularConstraints(
    const AtomConstraintGenerator& atom_constraint_generator,
    const BondConstraintGenerator& bond_constraint_generator,
    const EnvironmentConstraintGenerator& environment_constraint_generator,
    unsigned environment_radius = 2) :
    atom_constraint_generator(atom_constraint_generator),
    bond_constraint_generator(bond_constraint_generator),
    environment_constraint_generator(environment_constraint_generator),
    environment_generator(environment_radius, LoneAtomHashes) {};

  void GenerateAtomConstraints(
    const RDKit::ROMol& molecule,
    const AtomConstraintGenerator& atom_constraint_gen) {
    if (!atom_constraint_gen) {
      return;
    };
    atom_constraints.clear();
    for (const RDKit::Atom* atom : molecule.atoms()) {
      Tag atom_tag = GetTag(atom);
      std::optional<AtomConstraint> atom_constraint = atom_constraint_gen(
        {NULL_ATOM_KEY, AtomKey(atom)});
      if (!atom_constraint) {
        continue;
      };
      SetConstraint(atom_tag, std::move(*atom_constraint), atom_constraints);
    };
  };

  void GenerateBondConstraints(
    const RDKit::ROMol& molecule,
    const BondConstraintGenerator& bond_constraint_gen) {
    if (!bond_constraint_gen) {
      return;
    };
    bond_constraints.clear();
    for (const RDKit::Bond* bond : molecule.bonds()) {
      Tag bond_tag = GetTag(bond);
      std::optional<BondConstraint> bond_constraint = bond_constraint_gen(
        {NULL_BOND_KEY, BondKey(bond)});
      if (!bond_constraint) {
        continue;
      };
      SetConstraint(bond_tag, std::move(*bond_constraint), bond_constraints);
    };
  };

  void GenerateEnvironmentConstraints(
    const RDKit::ROMol& molecule,
    const EnvironmentConstraintGenerator& environment_constraint_gen) {
    if (!environment_constraint_gen) {
      return;
    };
    environment_constraints.clear();
    for (const RDKit::Atom* atom : molecule.atoms()) {
      Tag atom_tag = GetTag(atom);
      std::optional<EnvironmentConstraint> environment_constraint = 
        environment_constraint_gen(
          {NULL_ENVIRONMENT_KEY, environment_generator.Key(atom)});
      if (!environment_constraint) {
        continue;
      };
      SetConstraint(
        atom_tag, std::move(*environment_constraint), environment_constraints);
    };
  };

  void GenerateConstraints(const RDKit::ROMol& molecule) {
    GenerateAtomConstraints(molecule, atom_constraint_generator);
    GenerateBondConstraints(molecule, bond_constraint_generator);
    GenerateEnvironmentConstraints(molecule, environment_constraint_generator);
  };

  bool SetAtomConstraint(
    Tag atom_tag,
    const AtomConstraint& atom_constraint,
    bool replace = true) {
    return SetConstraint(atom_tag, atom_constraint, atom_constraints, replace);
  };

  bool SetBondConstraint(
    Tag bond_tag,
    const BondConstraint& bond_constraint,
    bool replace = true) {
    return SetConstraint(bond_tag, bond_constraint, bond_constraints, replace);
  };

  bool SetEnvironmentConstraint(
    Tag atom_tag,
    const EnvironmentConstraint& environment_constraint,
    bool replace = true) {
    return SetConstraint(
      atom_tag, environment_constraint, environment_constraints, replace);
  };

  void SetMinCycleSize(unsigned new_min_cycle_size) {
    min_cycle_size = new_min_cycle_size;
  };

  void SetMaxCycleSize(unsigned new_max_cycle_size) {
    max_cycle_size = new_max_cycle_size;
  };

  bool UpdateAtomConstraint(Tag atom_tag, const AtomKeyChange& atom_key_change) {
    if (!atom_constraint_generator) {
      return false;
    };
    std::optional<AtomConstraint> new_atom_constraint = 
      atom_constraint_generator(atom_key_change);
    if (!new_atom_constraint) {
      return false;
    };
    return SetConstraint(
      atom_tag, std::move(*new_atom_constraint), atom_constraints);
  };

  bool UpdateBondConstraint(Tag bond_tag, const BondKeyChange& bond_key_change) {
    if (!bond_constraint_generator) {
      return false;
    };
    std::optional<BondConstraint> new_bond_constraint =
      bond_constraint_generator(bond_key_change);
    if (!new_bond_constraint) {
      return false;
    };
    return SetConstraint(
      bond_tag, std::move(*new_bond_constraint), bond_constraints);
  };

  bool UpdateEnvironmentConstraint(
    Tag atom_tag, const EnvironmentKeyChange& environment_key_change) {
    if (!environment_constraint_generator) {
      return false;
    };
    std::optional<EnvironmentConstraint> new_environment_constraint = 
      environment_constraint_generator(environment_key_change);
    if (!new_environment_constraint) {
      return false;
    };
    return SetConstraint(
      atom_tag, std::move(*new_environment_constraint), environment_constraints);
  };

  bool UpdateConstraints(
    const RDKit::ROMol& molecule,
    const MolecularPerturbation& perturbation) {
    if (!HasConstraintGenerators()) {
      return false;
    };
    MolecularGraphProjection projection (molecule);
    perturbation.ProjectMolecularGraph(projection);
    std::size_t n_updated_keys = 0;
    if (atom_constraint_generator || bond_constraint_generator) {
      auto [atom_key_changes, bond_key_changes] =
        projection.MolecularKeyChanges(molecule, !!bond_constraint_generator);
      for (const auto& [atom_tag, atom_key_change] : atom_key_changes) {
        n_updated_keys += UpdateAtomConstraint(atom_tag, atom_key_change);
      };
      for (const auto& [bond_tag, bond_key_change] : bond_key_changes) {
        n_updated_keys += UpdateBondConstraint(bond_tag, bond_key_change);
      };
    };
    if (environment_constraint_generator) {
      auto environment_key_changes = projection.EnvironmentKeyChanges(
        molecule, environment_generator.GetEnvironmentRadius());
      for (const auto& [atom_tag, ekc] : environment_key_changes) {
        n_updated_keys += UpdateEnvironmentConstraint(atom_tag, ekc);
      };
    };
    return n_updated_keys;
  };

  bool ClearAtomConstraint(Tag atom_tag) {
    return atom_constraints.erase(atom_tag);
  };

  bool ClearBondConstraint(Tag bond_tag) {
    return bond_constraints.erase(bond_tag);
  };

  bool ClearEnvironmentConstraint(Tag atom_tag) {
    return environment_constraints.erase(atom_tag);
  };

  void ClearAtomConstraints() {
    atom_constraints.clear();
  };

  void ClearBondConstraints() {
    bond_constraints.clear();
  };

  void ClearEnvironmentConstraints() {
    environment_constraints.clear();
  };

  void ClearCyclicityConstraints() {
    min_cycle_size = 3;
    max_cycle_size = std::numeric_limits<unsigned>::max();
  };

  void ClearConstraints() {
    atom_constraints.clear();
    bond_constraints.clear();
    environment_constraints.clear();
    ClearCyclicityConstraints();
  };

  void Clear() {
    atom_constraint_generator = nullptr;
    bond_constraint_generator = nullptr;
    environment_constraint_generator = nullptr;
    ClearConstraints();
  };

  const AtomConstraint* GetAtomConstraint(Tag atom_tag) const {
    return GetConstraint(atom_tag, atom_constraints);
  };

  const AtomConstraint* GetAtomConstraint(const RDKit::Atom* atom) const {
    auto [atom_tag, atom_is_tagged] = GetTagIfPresent(atom);
    if (!atom_is_tagged) {
      return nullptr;
    };
    return GetConstraint(atom_tag, atom_constraints);
  };

  const BondConstraint* GetBondConstraint(Tag bond_tag) const {
    return GetConstraint(bond_tag, bond_constraints);
  };

  const BondConstraint* GetBondConstraint(const RDKit::Bond* bond) const {
    auto [bond_tag, bond_is_tagged] = GetTagIfPresent(bond);
    if (!bond_is_tagged) {
      return nullptr;
    };
    return GetConstraint(bond_tag, bond_constraints);
  };

  const EnvironmentConstraint* GetEnvironmentConstraint(Tag atom_tag) const {
    return GetConstraint(atom_tag, environment_constraints);
  };

  const EnvironmentConstraint* GetEnvironmentConstraint(
    const RDKit::Atom* atom) const {
    auto [atom_tag, atom_is_tagged] = GetTagIfPresent(atom);
    if (!atom_is_tagged) {
      return nullptr;
    };
    return GetConstraint(atom_tag, environment_constraints);
  };

  unsigned GetMinCycleSize() const {
    return min_cycle_size;
  };

  unsigned GetMaxCycleSize() const {
    return max_cycle_size;
  };

  unsigned GetEnvironmentRadius() const {
    return environment_generator.GetEnvironmentRadius();
  };

  bool IsAllowed(Tag atom_tag, const AtomKey& atom_key) const {
    const AtomConstraint* atom_constraint = GetAtomConstraint(atom_tag);
    // If the atom has no constraints, everything is allowed.
    // If it does have a constraint, check if the AtomKey satisfies it.
    return !atom_constraint || (*atom_constraint)(atom_key);
  };

  bool IsAllowed(Tag bond_tag, const BondKey& bond_key) const {
    const BondConstraint* bond_constraint = GetBondConstraint(bond_tag);
    return !bond_constraint || (*bond_constraint)(bond_key);
  };

  bool IsAllowed(Tag atom_tag, const EnvironmentKey& environment_key) const {
    const EnvironmentConstraint* environment_constraint = 
      GetEnvironmentConstraint(atom_tag);
    return !environment_constraint || (*environment_constraint)(environment_key);
  };

  bool IsAllowed(Tag atom_tag, const AtomKeyChange& atom_key_change) const {
    auto [atom_constraint, updated] = 
      UpdatedAtomConstraint(atom_tag, atom_key_change);
    return !atom_constraint || (*atom_constraint)(atom_key_change.second);
  };

  bool IsAllowed(Tag bond_tag, const BondKeyChange& bond_key_change) const {
    auto [bond_constraint, updated] = 
      UpdatedBondConstraint(bond_tag, bond_key_change);
    return !bond_constraint || (*bond_constraint)(bond_key_change.second);
  };

  bool IsAllowed(
    Tag atom_tag, const EnvironmentKeyChange& environment_key_change) const {
    auto [environment_constraint, updated] = 
      UpdatedEnvironmentConstraint(atom_tag, environment_key_change);
    return !environment_constraint ||
      (*environment_constraint)(environment_key_change.second);
  };

  bool IsAllowed(
    const RDKit::ROMol& molecule,
    const MolecularPerturbation& perturbation) const {
    if (!(*this)) {
      return true;
    };
    MolecularGraphProjection projection (molecule);
    perturbation.ProjectMolecularGraph(projection);
    if (MustCheckAtomConstraints() || MustCheckBondConstraints()) {
      auto [atom_key_changes, bond_key_changes] =
        projection.MolecularKeyChanges(molecule, MustCheckBondConstraints());
      for (const auto& [atom_tag, atom_key_change] : atom_key_changes) {
        if (!IsAllowed(atom_tag, atom_key_change)) {
          return false;
        };
      };
      for (const auto& [bond_tag, bond_key_change] : bond_key_changes) {
        if (!IsAllowed(bond_tag, bond_key_change)) {
          return false;
        };
      };
    };
    if (MustCheckEnvironmentConstraints()) {
      auto environment_key_changes = projection.EnvironmentKeyChanges(
        molecule, environment_generator.GetEnvironmentRadius());
      for (const auto& [atom_tag, ekc] : environment_key_changes) {
        if (!IsAllowed(atom_tag, ekc)) {
          return false;
        };
      };
    };
    return CyclicityConstraintsSatisfied(projection);
  };

  bool UpdateIfAllowed(Tag atom_tag, const AtomKeyChange& atom_key_change) {
    auto [atom_constraint, updated] = 
      UpdatedAtomConstraint(atom_tag, atom_key_change);
    bool is_allowed = 
      !atom_constraint || (*atom_constraint)(atom_key_change.second);
    if (is_allowed && updated) {
      SetConstraint(atom_tag, std::move(*atom_constraint), atom_constraints);
    };
    return is_allowed;
  };

  bool UpdateIfAllowed(Tag bond_tag, const BondKeyChange& bond_key_change) {
    auto [bond_constraint, updated] = 
      UpdatedBondConstraint(bond_tag, bond_key_change);
    bool is_allowed = 
      !bond_constraint || (*bond_constraint)(bond_key_change.second);
    if (is_allowed && updated) {
      SetConstraint(bond_tag, std::move(*bond_constraint), bond_constraints);
    };
    return is_allowed;
  };

  bool UpdateIfAllowed(
    Tag atom_tag, const EnvironmentKeyChange& environment_key_change) {
    auto [environment_constraint, updated] = 
      UpdatedEnvironmentConstraint(atom_tag, environment_key_change);
    bool is_allowed = !environment_constraint ||
      (*environment_constraint)(environment_key_change.second);
    if (is_allowed && updated) {
      SetConstraint(
        atom_tag, std::move(*environment_constraint), environment_constraints);
    };
    return is_allowed;
  };

  bool UpdateIfAllowed(
    const RDKit::ROMol& molecule,
    const MolecularPerturbation& perturbation) {
    if (!HasConstraintGenerators()) {
      return IsAllowed(molecule, perturbation);
    };
    std::vector<std::pair<Tag, std::shared_ptr<const AtomConstraint>>>
      updated_atom_constraints;
    std::vector<std::pair<Tag, std::shared_ptr<const BondConstraint>>>
      updated_bond_constraints;
    std::vector<std::pair<Tag, std::shared_ptr<const EnvironmentConstraint>>>
      updated_environment_constraints;
    MolecularGraphProjection projection (molecule);
    perturbation.ProjectMolecularGraph(projection);
    if (MustCheckAtomConstraints() || MustCheckBondConstraints()) {
      auto [atom_key_changes, bond_key_changes] =
        projection.MolecularKeyChanges(molecule, MustCheckBondConstraints());
      for (const auto& [atom_tag, atom_key_change] : atom_key_changes) {
        auto [atom_constraint, updated] = 
          UpdatedAtomConstraint(atom_tag, atom_key_change);
        if (atom_constraint && !(*atom_constraint)(atom_key_change.second)) {
          return false;
        };
        if (updated) {
          updated_atom_constraints.emplace_back(
            atom_tag, std::move(atom_constraint));
        };
      };
      for (const auto& [bond_tag, bond_key_change] : bond_key_changes) {
        auto [bond_constraint, updated] = 
          UpdatedBondConstraint(bond_tag, bond_key_change);
        if (bond_constraint && !(*bond_constraint)(bond_key_change.second)) {
          return false;
        };
        if (updated) {
          updated_bond_constraints.emplace_back(
            bond_tag, std::move(bond_constraint));
        };
      };
    };
    if (MustCheckEnvironmentConstraints()) {
      auto environment_key_changes = projection.EnvironmentKeyChanges(
        molecule, environment_generator.GetEnvironmentRadius());
      for (const auto& [atom_tag, ekc] : environment_key_changes) {
        auto [environment_constraint, updated] =
          UpdatedEnvironmentConstraint(atom_tag, ekc);
        if (environment_constraint && !(*environment_constraint)(ekc.second)) {
          return false;
        };
        if (updated) {
          updated_environment_constraints.emplace_back(
            atom_tag, std::move(environment_constraint));
        };
      };
    };
    if (!CyclicityConstraintsSatisfied(projection)) {
      return false;
    };
    for (const auto& [atom_tag, atom_constraint] : updated_atom_constraints) {
      SetConstraint(atom_tag, std::move(*atom_constraint), atom_constraints);
    };
    for (const auto& [bond_tag, bond_constraint] : updated_bond_constraints) {
      SetConstraint(bond_tag, std::move(*bond_constraint), bond_constraints);
    };
    for (const auto& [atom_tag, envc] : updated_environment_constraints) {
      SetConstraint(atom_tag, std::move(*envc), environment_constraints);
    };
    return true;
  };

  bool HasAtomConstraintGenerator() const {
    return !!atom_constraint_generator;
  };

  bool HasBondConstraintGenerator() const {
    return !!bond_constraint_generator;
  };

  bool HasEnvironmentConstraintGenerator() const {
    return !!environment_constraint_generator;
  };

  bool HasConstraintGenerators() const {
    return atom_constraint_generator || 
      bond_constraint_generator ||
      environment_constraint_generator;
  };

  bool MustCheckAtomConstraints() const {
    return !atom_constraints.empty() || atom_constraint_generator;
  };

  bool MustCheckBondConstraints() const {
    return !bond_constraints.empty() || bond_constraint_generator;
  };

  bool MustCheckEnvironmentConstraints() const {
    return !environment_constraints.empty() || environment_constraint_generator;
  };

  const CircularAtomicEnvironmentGenerator& GetEnvironmentGenerator() const {
    return environment_generator;
  };

  std::size_t Size() const {
    return atom_constraints.size() + 
      bond_constraints.size() + 
      environment_constraints.size();
  };

  bool HasCyclicityConstraints() const {
    return min_cycle_size > 3 || 
      max_cycle_size < std::numeric_limits<unsigned>::max();
  };

  explicit operator bool() const {
    return MustCheckAtomConstraints() ||
      MustCheckBondConstraints() ||
      MustCheckEnvironmentConstraints() ||
      HasCyclicityConstraints();
  };
};


bool ValenceConstraint(const AtomKey& atom_key) {
  return IsValenceBelowMax(atom_key.atomic_number, atom_key.valence);
};

std::optional<MolecularConstraints::AtomConstraint> ValenceConstraintGenerator(
  const AtomKeyChange& atom_key_change) {
  // If the atomic number didn't change neither did the valence constraint.
  if (atom_key_change.first.atomic_number ==
    atom_key_change.second.atomic_number) {
    return std::nullopt;
  };
  // Null atoms have no constraints.
  if (atom_key_change.second == NULL_ATOM_KEY) {
    return nullptr;
  } else {
    return ValenceConstraint;
  };
};


class AtomKeyConstraint {
  const ChemicalDictionary* dictionary;
public:
  AtomKeyConstraint(const ChemicalDictionary* dictionary) : 
    dictionary(dictionary) {};
  
  bool operator()(const AtomKey& atom_key) const {
    // Null keys correspond to deletions and are always allowed, unless 
    // otherwise specified.
    if (atom_key == NULL_ATOM_KEY) {
      return true;
    };
    return !dictionary->IsForeignAtom(atom_key);
  };
};

class AtomKeyConstraintGenerator {
  const ChemicalDictionary* dictionary;
public:
  AtomKeyConstraintGenerator(const ChemicalDictionary* dictionary) :
    dictionary(dictionary) {};

  std::optional<MolecularConstraints::AtomConstraint> operator()(
    const AtomKeyChange& atom_key_change) const {
    if (atom_key_change.first == atom_key_change.second) {
      return std::nullopt;
    };
    if (atom_key_change.second == NULL_ATOM_KEY) {
      return nullptr;
    } else {
      return AtomKeyConstraint(dictionary);
    };
  };
};


class BondKeyConstraint {
  const ChemicalDictionary* dictionary;
public:
  BondKeyConstraint(const ChemicalDictionary* dictionary) : 
    dictionary(dictionary) {};
  
  bool operator()(const BondKey& bond_key) const {
    if (bond_key == NULL_BOND_KEY) {
      return true;
    };
    return !dictionary->IsForeignBond(bond_key);
  };
};

class BondKeyConstraintGenerator {
  const ChemicalDictionary* dictionary;
public:
  BondKeyConstraintGenerator(const ChemicalDictionary* dictionary) :
    dictionary(dictionary) {};

  std::optional<MolecularConstraints::BondConstraint> operator()(
    const BondKeyChange& bond_key_change) const {
    if (bond_key_change.first == bond_key_change.second) {
      return std::nullopt;
    };
    if (bond_key_change.second == NULL_BOND_KEY) {
      return nullptr;
    } else {
      return BondKeyConstraint(dictionary);
    };
  };
};


class EnvironmentKeyConstraint {
  const ChemicalDictionary* dictionary;
public:
  EnvironmentKeyConstraint(const ChemicalDictionary* dictionary) : 
    dictionary(dictionary) {};
  
  bool operator()(const EnvironmentKey& environment_key) const {
    if (environment_key == NULL_ENVIRONMENT_KEY) {
      return true;
    };
    return !dictionary->IsForeignEnvironment(environment_key);
  };
};

class EnvironmentKeyConstraintGenerator {
  const ChemicalDictionary* dictionary;
public:
  EnvironmentKeyConstraintGenerator(const ChemicalDictionary* dictionary) :
    dictionary(dictionary) {};

  std::optional<MolecularConstraints::EnvironmentConstraint> operator()(
    const EnvironmentKeyChange& environment_key_change) const {
    if (environment_key_change.first == environment_key_change.second) {
      return std::nullopt;
    };
    if (environment_key_change.second == NULL_ENVIRONMENT_KEY) {
      return nullptr;
    } else {
      return EnvironmentKeyConstraint(dictionary);
    };
  };
};

#endif // !_MOLECULAR_CONSTRAINTS_HPP_
