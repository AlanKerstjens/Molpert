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
  // The bool in the value is a "static" flag.
  // When set the constraint won't be settable or updateable.
  std::map<Tag, std::pair<AtomConstraint, bool>> atom_constraints;
  std::map<Tag, std::pair<BondConstraint, bool>> bond_constraints;
  std::map<Tag, std::pair<EnvironmentConstraint, bool>> environment_constraints;
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
    std::map<Tag, std::pair<Constraint, bool>>& constraints,
    bool replace = true,
    bool make_static = false) const {
    // Check if we already have a constraint for this Tag.
    auto it = constraints.find(tag);
    // If we don't and the constraint is not null we set it.
    if (it == constraints.cend()) {
      if (!constraint) {
        return false;
      };
      constraints.try_emplace(tag, constraint, make_static);
      return true;
    };
    // If we do, check if we are allowed to change it.
    if (!replace || it->second.second) {
      return false;
    };
    // If we are but the new constraint is null we erase the existing one.
    if (!constraint) {
      constraints.erase(it);
      return true;
    };
    // Otherwise replace it with the new constraint.
    it->second = {constraint, make_static};
    return true;
  };

  template <class Constraint>
  const Constraint* GetConstraint(
    Tag tag,
    const std::map<Tag, std::pair<Constraint, bool>>& constraints) const {
    auto it = constraints.find(tag);
    // Return a null constraint if the tag has no constraints.
    return it != constraints.cend() ? &it->second.first : nullptr;
  };

  template <class Constraint, class KeyChange, class ConstraintGenerator>
  std::pair<std::shared_ptr<const Constraint>, bool> UpdatedConstraint(
    Tag tag,
    const std::map<Tag, std::pair<Constraint, bool>>& constraints,
    const KeyChange& key_change,
    const ConstraintGenerator& constraint_generator) const {
    std::shared_ptr<const Constraint> updated_constraint;
    // Check if we already have a constraint for this tag.
    const Constraint* constraint = nullptr;
    bool is_static = false;
    auto it = constraints.find(tag);
    if (it != constraints.cend()) {
      constraint = &it->second.first;
      is_static = it->second.second;
    };
    // If the existing constraint is static or we don't have a constraint
    // generator we return the existing constraint, which may be null.
    if (is_static || !constraint_generator) {
      updated_constraint.reset(constraint, boost::null_deleter());
      return {updated_constraint, false};
    };
    // If not, we can try generating a new constraint.
    std::optional<Constraint> new_constraint = constraint_generator(key_change);
    // If no new constraint was generated we don't update anything.
    if (!new_constraint) {
      updated_constraint.reset(constraint, boost::null_deleter());
      return {updated_constraint, false};
    };
    // If a non-null constraint was generated take ownership of it.
    if (*new_constraint) {
      updated_constraint.reset(new Constraint(std::move(*new_constraint)));
    };
    // Return the new constraint.
    return {updated_constraint, true};
  };

  template <class Constraint>
  bool ClearConstraint(
    Tag tag,
    std::map<Tag, std::pair<Constraint, bool>>& constraints,
    bool clear_static = false) const {
    auto it = constraints.find(tag);
    if (it == constraints.cend()) {
      return false;
    };
    bool is_static = it->second.second;
    if (is_static && !clear_static) {
      return false;
    };
    constraints.erase(it);
    return true;
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
    environment_generator(environment_radius) {};

  bool SetAtomConstraint(
    Tag atom_tag,
    const AtomConstraint& atom_constraint,
    bool replace = true,
    bool make_static = false) {
    return SetConstraint(
      atom_tag, atom_constraint, atom_constraints, replace, make_static);
  };

  bool SetBondConstraint(
    Tag bond_tag,
    const BondConstraint& bond_constraint,
    bool replace = true,
    bool make_static = false) {
    return SetConstraint(
      bond_tag, bond_constraint, bond_constraints, replace, make_static);
  };

  bool SetEnvironmentConstraint(
    Tag atom_tag,
    const EnvironmentConstraint& environment_constraint,
    bool replace = true,
    bool make_static = false) {
    return SetConstraint(atom_tag, 
      environment_constraint, environment_constraints, replace, make_static);
  };

  void SetMinCycleSize(unsigned new_min_cycle_size) {
    min_cycle_size = new_min_cycle_size;
  };

  void SetMaxCycleSize(unsigned new_max_cycle_size) {
    max_cycle_size = new_max_cycle_size;
  };

  void GenerateAtomConstraint(const RDKit::Atom* atom) {
    if (!atom_constraint_generator) {
      return;
    };
    std::optional<AtomConstraint> atom_constraint = atom_constraint_generator(
      {NULL_ATOM_KEY, AtomKey(atom)});
    if (!atom_constraint) {
      return;
    };
    SetAtomConstraint(GetTag(atom), std::move(*atom_constraint));
  };

  void GenerateAtomConstraints(const RDKit::ROMol& molecule) {
    if (!atom_constraint_generator) {
      return;
    };
    for (const RDKit::Atom* atom : molecule.atoms()) {
      GenerateAtomConstraint(atom);
    };
  };

  void GenerateBondConstraint(const RDKit::Bond* bond) {
    if (!bond_constraint_generator) {
      return;
    };
    std::optional<BondConstraint> bond_constraint = bond_constraint_generator(
      {NULL_BOND_KEY, BondKey(bond)});
    if (!bond_constraint) {
      return;
    };
    SetBondConstraint(GetTag(bond), std::move(*bond_constraint));
  };

  void GenerateBondConstraints(const RDKit::ROMol& molecule) {
    if (!bond_constraint_generator) {
      return;
    };
    for (const RDKit::Bond* bond : molecule.bonds()) {
      GenerateBondConstraint(bond);
    };
  };

  void GenerateEnvironmentConstraint(
    const RDKit::ROMol& molecule,
    const RDKit::Atom* atom,
    bool recalculate_atom_hashes = false) {
    if (!environment_constraint_generator) {
      return;
    };
    if (recalculate_atom_hashes) {
      environment_generator.CalculateAtomHashes(molecule);
    };
    std::optional<EnvironmentConstraint> environment_constraint = 
      environment_constraint_generator(
        {NULL_ENVIRONMENT_KEY, environment_generator.Key(atom)});
    if (!environment_constraint) {
      return;
    };
    SetEnvironmentConstraint(GetTag(atom), std::move(*environment_constraint));
  };

  void GenerateEnvironmentConstraints(const RDKit::ROMol& molecule) {
    if (!environment_constraint_generator) {
      return;
    };
    environment_generator.CalculateAtomHashes(molecule);
    for (const RDKit::Atom* atom : molecule.atoms()) {
      GenerateEnvironmentConstraint(molecule, atom, false);
    };
  };

  void GenerateConstraints(const RDKit::ROMol& molecule) {
    GenerateAtomConstraints(molecule);
    GenerateBondConstraints(molecule);
    GenerateEnvironmentConstraints(molecule);
  };

  std::pair<std::shared_ptr<const AtomConstraint>, bool>
  UpdatedAtomConstraint(
    Tag atom_tag, const AtomKeyChange& atom_key_change) const {
    return UpdatedConstraint(
      atom_tag, atom_constraints, atom_key_change, atom_constraint_generator);
  };

  std::pair<std::shared_ptr<const BondConstraint>, bool>
  UpdatedBondConstraint(
    Tag bond_tag, const BondKeyChange& bond_key_change) const {
    return UpdatedConstraint(
      bond_tag, bond_constraints, bond_key_change, bond_constraint_generator);
  };

  std::pair<std::shared_ptr<const EnvironmentConstraint>, bool>
  UpdatedEnvironmentConstraint(
    Tag atom_tag, const EnvironmentKeyChange& environment_key_change) const {
    return UpdatedConstraint(atom_tag, environment_constraints, 
      environment_key_change, environment_constraint_generator);
  };

  bool UpdateAtomConstraint(
    Tag atom_tag, const AtomKeyChange& atom_key_change) {
    auto [atom_constraint, updated] = UpdatedAtomConstraint(
      atom_tag, atom_key_change);
    if (!updated) {
      return false;
    };
    return SetAtomConstraint(atom_tag, std::move(*atom_constraint));
  };

  bool UpdateBondConstraint(
    Tag bond_tag, const BondKeyChange& bond_key_change) {
    auto [bond_constraint, updated] = UpdatedBondConstraint(
      bond_tag, bond_key_change);
    if (!updated) {
      false;
    };
    return SetBondConstraint(bond_tag, std::move(*bond_constraint));
  };

  bool UpdateEnvironmentConstraint(
    Tag atom_tag, const EnvironmentKeyChange& environment_key_change) {
    auto [environment_constraint, updated] = UpdatedEnvironmentConstraint(
      atom_tag, environment_key_change);
    if (!updated) {
      false;
    };
    return SetEnvironmentConstraint(
      atom_tag, std::move(*environment_constraint));
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

  bool ClearAtomConstraint(Tag atom_tag, bool clear_static = false) {
    return ClearConstraint(atom_tag, atom_constraints, clear_static);
  };

  bool ClearBondConstraint(Tag bond_tag, bool clear_static = false) {
    return ClearConstraint(bond_tag, bond_constraints, clear_static);
  };

  bool ClearEnvironmentConstraint(Tag atom_tag, bool clear_static = false) {
    return ClearConstraint(atom_tag, environment_constraints, clear_static);
  };

  void ClearAtomConstraints(bool clear_static = false) {
    if (clear_static) {
      return atom_constraints.clear();
    };
    for (auto it = atom_constraints.cbegin(); it != atom_constraints.cend();) {
      bool is_static = it->second.second;
      if (!is_static) {
        it = atom_constraints.erase(it);
      } else {
        ++it;
      };
    };
  };

  void ClearBondConstraints(bool clear_static = false) {
    if (clear_static) {
      return bond_constraints.clear();
    };
    for (auto it = bond_constraints.cbegin(); it != bond_constraints.cend();) {
      bool is_static = it->second.second;
      if (!is_static) {
        it = bond_constraints.erase(it);
      } else {
        ++it;
      };
    };
  };

  void ClearEnvironmentConstraints(bool clear_static = false) {
    if (clear_static) {
      return environment_constraints.clear();
    };
    for (auto it = environment_constraints.cbegin();
      it != environment_constraints.cend();) {
      bool is_static = it->second.second;
      if (!is_static) {
        it = environment_constraints.erase(it);
      } else {
        ++it;
      };
    };
  };

  void ClearCyclicityConstraints() {
    min_cycle_size = 3;
    max_cycle_size = std::numeric_limits<unsigned>::max();
  };

  void ClearConstraints(
    bool clear_static = false,
    bool clear_cyclicity_constraints = true) {
    ClearAtomConstraints(clear_static);
    ClearBondConstraints(clear_static);
    ClearEnvironmentConstraints(clear_static);
    if (clear_cyclicity_constraints) {
      ClearCyclicityConstraints();  
    };
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
    auto [atom_constraint, updated] = UpdatedAtomConstraint(
      atom_tag, atom_key_change);
    return !atom_constraint || (*atom_constraint)(atom_key_change.second);
  };

  bool IsAllowed(Tag bond_tag, const BondKeyChange& bond_key_change) const {
    auto [bond_constraint, updated] = UpdatedBondConstraint(
      bond_tag, bond_key_change);
    return !bond_constraint || (*bond_constraint)(bond_key_change.second);
  };

  bool IsAllowed(
    Tag atom_tag, const EnvironmentKeyChange& environment_key_change) const {
    auto [environment_constraint, updated] = UpdatedEnvironmentConstraint(
      atom_tag, environment_key_change);
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
      if (MustCheckAtomConstraints()) {
        for (const auto& [atom_tag, atom_key_change] : atom_key_changes) {
          if (!IsAllowed(atom_tag, atom_key_change)) {
            return false;
          };
        };
      };
      if (MustCheckBondConstraints()) {
        for (const auto& [bond_tag, bond_key_change] : bond_key_changes) {
          if (!IsAllowed(bond_tag, bond_key_change)) {
            return false;
          };
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
      SetAtomConstraint(atom_tag, std::move(*atom_constraint));
    };
    return is_allowed;
  };

  bool UpdateIfAllowed(Tag bond_tag, const BondKeyChange& bond_key_change) {
    auto [bond_constraint, updated] = 
      UpdatedBondConstraint(bond_tag, bond_key_change);
    bool is_allowed = 
      !bond_constraint || (*bond_constraint)(bond_key_change.second);
    if (is_allowed && updated) {
      SetBondConstraint(bond_tag, std::move(*bond_constraint));
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
      SetEnvironmentConstraint(atom_tag, std::move(*environment_constraint));
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
      if (MustCheckAtomConstraints()) {
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
      };
      if (MustCheckBondConstraints()) {
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
      SetAtomConstraint(atom_tag, std::move(*atom_constraint));
    };
    for (const auto& [bond_tag, bond_constraint] : updated_bond_constraints) {
      SetBondConstraint(bond_tag, std::move(*bond_constraint));
    };
    for (const auto& [atom_tag, envc] : updated_environment_constraints) {
      SetEnvironmentConstraint(atom_tag, std::move(*envc));
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
