#ifndef _MOLECULAR_CONSTRAINTS_HPP_
#define _MOLECULAR_CONSTRAINTS_HPP_

#include "MolecularPerturbations.hpp"
#include <boost/core/null_deleter.hpp>

#ifdef BUILD_PYTHON_BINDINGS

#include <boost/python.hpp>

// Our constraints may be C++ functions (i.e. std::function) or Python functions 
// (i.e. a callable boost::python::object). They take as input a pointer, namely
// a const RDKit::Atom* or const RDKit::Bond*. To quote the Boost.Python
// documentation: "Normally, when passing pointers to Python callbacks, the 
// pointee is copied to ensure that the Python object never holds a dangling 
// reference". This means that the RDKit::Atom and RDKit::Bond copy constructors 
// are invoked. These copy constructors DO NOT copy the pointer to the owning 
// molecule, which makes the copies virtually useless. "To specify that the new 
// Python object should merely contain a copy of a pointer p, the user can pass 
// ptr(p) instead of passing p directly". This wrapper class allows us to store 
// a Python callback and to pass the pointer as described.
template <class T>
class ConstraintWrapper {
  typedef std::function<bool(const T&)> CPPConstraint;
  typedef boost::python::object PythonConstraint;

  std::variant<CPPConstraint, PythonConstraint> constraint;

public:
  ConstraintWrapper(
    const CPPConstraint& function) : 
    constraint(function) {
  };
  ConstraintWrapper(
    const PythonConstraint& object) {
    if (object.is_none()) {
      constraint = nullptr;
      return;
    };
    if (!PyCallable_Check(object.ptr())) {
      PyErr_SetString(PyExc_TypeError, "Expected a callable object");
      throw boost::python::error_already_set();
    };
    constraint = object;
  };

  bool operator()(const T& t) const {
    if (std::holds_alternative<CPPConstraint>(constraint)) {
      const CPPConstraint& c = std::get<CPPConstraint>(constraint);
      return c(t);
    };
    const PythonConstraint& c = std::get<PythonConstraint>(constraint);
    boost::python::object result = c(boost::python::ptr(t));
    return boost::python::extract<bool>(result);
  };

  explicit operator bool() const {
    if (std::holds_alternative<CPPConstraint>(constraint)) {
      const CPPConstraint& c = std::get<CPPConstraint>(constraint);
      return !!c;
    };
    const PythonConstraint& c = std::get<PythonConstraint>(constraint);
    return !c.is_none();
  };

};

typedef ConstraintWrapper<const RDKit::Atom*> AtomConstraint;
typedef ConstraintWrapper<const RDKit::Bond*> BondConstraint;
typedef ConstraintWrapper<const RDKit::ROMol*> MoleculeConstraint;

#else

// Evaluates an atom, bond or molecule and returns a boolean indicating if it's
// allowed (true) or not (false).
typedef std::function<bool(const RDKit::Atom*)> AtomConstraint;
typedef std::function<bool(const RDKit::Bond*)> BondConstraint;
typedef std::function<bool(const RDKit::ROMol*)> MoleculeConstraint;

#endif

// Constraint generators compare two atoms or bonds (prior and  posterior) and 
// generate new constraints for the posterior one. The return type is a 
// Constraint that should evaluate to false (nullptr or None) if no constraint 
// applies or the constraint has been erased.
typedef std::function<AtomConstraint(const RDKit::Atom*, const RDKit::Atom*)
> AtomConstraintGenerator;
typedef std::function<BondConstraint(const RDKit::Bond*, const RDKit::Bond*)
> BondConstraintGenerator;


class MolecularConstraints {
  // The bool in the value is a "static" flag.
  // When set the constraint won't be settable or updateable.
  std::map<Tag, std::pair<AtomConstraint, bool>> atom_constraints;
  std::map<Tag, std::pair<BondConstraint, bool>> bond_constraints;
  std::vector<std::pair<MoleculeConstraint, bool>> molecule_constraints;
  AtomConstraintGenerator atom_constraint_generator;
  BondConstraintGenerator bond_constraint_generator;

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

  template <class T, class Constraint, class ConstraintGenerator>
  std::tuple<Tag, std::shared_ptr<const Constraint>, bool> UpdatedConstraint(
    const T* prior,
    const T* posterior,
    const std::map<Tag, std::pair<Constraint, bool>>& constraints,
    const ConstraintGenerator& constraint_generator) const {
    assert(prior || posterior); // If both are null something went wrong.
    std::shared_ptr<const Constraint> updated_constraint;
    Tag tag = prior ? GetTag(prior) : GetTag(posterior);
    // Null posteriors don't get constraints.
    if (!posterior) {
      return {tag, updated_constraint, true};
    };
    // Check if we already had a constraint for the prior.
    const Constraint* constraint = nullptr;
    bool is_static = false;
    if (prior) {
      auto it = constraints.find(tag);
      if (it != constraints.cend()) {
        constraint = &it->second.first;
        is_static = it->second.second;
      };
    };
    // If the existing constraint is static or we don't have a constraint
    // generator we return the existing constraint, which may be null.
    if (is_static || !constraint_generator) {
      updated_constraint.reset(constraint, boost::null_deleter());
      return {tag, updated_constraint, false};
    };
    // If not, we can try generating a new constraint.
    Constraint new_constraint = constraint_generator(prior, posterior);
    // If a non-null constraint was generated take ownership of it.
    if (new_constraint) {
      updated_constraint.reset(new Constraint(std::move(new_constraint)));
    };
    // Return the new constraint.
    return {tag, updated_constraint, true};
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

public:
  MolecularConstraints() = default;
  MolecularConstraints(
    const AtomConstraintGenerator& atom_constraint_generator) :
    atom_constraint_generator(atom_constraint_generator) {};
  MolecularConstraints(
    const BondConstraintGenerator& bond_constraint_generator) :
    bond_constraint_generator(bond_constraint_generator) {};
  MolecularConstraints(
    const AtomConstraintGenerator& atom_constraint_generator,
    const BondConstraintGenerator& bond_constraint_generator) :
    atom_constraint_generator(atom_constraint_generator),
    bond_constraint_generator(bond_constraint_generator) {};

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

  void SetMoleculeConstraint(
    const MoleculeConstraint& molecule_constraint,
    bool make_static = false) {
    if (molecule_constraint) {
      molecule_constraints.emplace_back(molecule_constraint, make_static);
    };
  };

  void SetAtomConstraintGenerator(
    const AtomConstraintGenerator& new_atom_constraint_generator) {
    atom_constraint_generator = new_atom_constraint_generator;
  };

  void SetBondConstraintGenerator(
    const BondConstraintGenerator& new_bond_constraint_generator) {
    bond_constraint_generator = new_bond_constraint_generator;
  };

  void GenerateAtomConstraint(const RDKit::Atom* atom) {
    if (!atom_constraint_generator) {
      return;
    };
    AtomConstraint atom_constraint = atom_constraint_generator(nullptr, atom);
    if (!atom_constraint) {
      return;
    };
    SetAtomConstraint(GetTag(atom), std::move(atom_constraint));
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
    BondConstraint bond_constraint = bond_constraint_generator(nullptr, bond);
    if (!bond_constraint) {
      return;
    };
    SetBondConstraint(GetTag(bond), std::move(bond_constraint));
  };

  void GenerateBondConstraints(const RDKit::ROMol& molecule) {
    if (!bond_constraint_generator) {
      return;
    };
    for (const RDKit::Bond* bond : molecule.bonds()) {
      GenerateBondConstraint(bond);
    };
  };

  void GenerateConstraints(const RDKit::ROMol& molecule) {
    GenerateAtomConstraints(molecule);
    GenerateBondConstraints(molecule);
  };

  std::tuple<Tag, std::shared_ptr<const AtomConstraint>, bool>
  UpdatedAtomConstraint(
    const RDKit::Atom* prior_atom,
    const RDKit::Atom* posterior_atom) const {
    return UpdatedConstraint(prior_atom, posterior_atom, 
      atom_constraints, atom_constraint_generator);
  };

  std::tuple<Tag, std::shared_ptr<const BondConstraint>, bool>
  UpdatedBondConstraint(
    const RDKit::Bond* prior_bond, 
    const RDKit::Bond* posterior_bond) const {
    return UpdatedConstraint(prior_bond, posterior_bond,
      bond_constraints, bond_constraint_generator);
  };

  bool UpdateAtomConstraint(
    const RDKit::Atom* prior_atom,
    const RDKit::Atom* posterior_atom) {
    auto [atom_tag, atom_constraint, updated] = UpdatedAtomConstraint(
      prior_atom, posterior_atom);
    if (!updated) {
      return false;
    };
    if (!atom_constraint) {
      return ClearAtomConstraint(atom_tag);
    };
    return SetAtomConstraint(atom_tag, std::move(*atom_constraint));
  };

  bool UpdateBondConstraint(
    const RDKit::Bond* prior_bond,
    const RDKit::Bond* posterior_bond) {
    auto [bond_tag, bond_constraint, updated] = UpdatedBondConstraint(
      prior_bond, posterior_bond);
    if (!updated) {
      false;
    };
    if (!bond_constraint) {
      return ClearBondConstraint(bond_tag);
    };
    return SetBondConstraint(bond_tag, std::move(*bond_constraint));
  };

  bool UpdateConstraints(
    const MolecularPerturbation& perturbation,
    const RDKit::ROMol& molecule) {
    if (!HasConstraintGenerators()) {
      return false;
    };
    RDKit::RWMol perturbed_molecule = perturbation(molecule);
    std::size_t n_updated_keys = 0;
    if (atom_constraint_generator) {
      auto atom_pairs = AtomPairs(molecule, perturbed_molecule);
      for (auto [prior_atom, posterior_atom] : atom_pairs) {
        n_updated_keys += UpdateAtomConstraint(prior_atom, posterior_atom);
      };
    };
    if (bond_constraint_generator) {
      auto bond_pairs = BondPairs(molecule, perturbed_molecule);
      for (auto [prior_bond, posterior_bond] : bond_pairs) {
        n_updated_keys += UpdateBondConstraint(prior_bond, posterior_bond);
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

  bool ClearMoleculeConstraint(
    std::size_t constraint_idx, bool clear_static = false) {
    if (constraint_idx >= molecule_constraints.size()) {
      return false;
    };
    if (clear_static || !molecule_constraints[constraint_idx].second) {
      molecule_constraints.erase(molecule_constraints.begin() + constraint_idx);
      return true;
    };
    return false;
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

  void ClearMoleculeConstraints(bool clear_static = false) {
    if (clear_static) {
      return molecule_constraints.clear();
    };
    auto it = std::remove_if(
      molecule_constraints.begin(), molecule_constraints.end(),
      [](const std::pair<MoleculeConstraint, bool>& constraint) {
        return !constraint.second;
      });
    molecule_constraints.erase(it, molecule_constraints.end());
  };

  void ClearConstraints(bool clear_static = false) {
    ClearAtomConstraints(clear_static);
    ClearBondConstraints(clear_static);
    ClearMoleculeConstraints(clear_static);
  };

  void Clear(bool clear_static = false) {
    atom_constraint_generator = nullptr;
    bond_constraint_generator = nullptr;
    ClearConstraints(clear_static);
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

  const MoleculeConstraint* GetMoleculeConstraint(
    std::size_t constraint_idx) const {
    return constraint_idx >= molecule_constraints.size() ? 
      nullptr : &molecule_constraints[constraint_idx].first;
  };

  bool IsAllowed(Tag atom_tag, const RDKit::Atom* atom) const {
    const AtomConstraint* atom_constraint = GetAtomConstraint(atom_tag);
    // If the atom has no constraints, everything is allowed.
    // If it does have a constraint, check if the AtomKey satisfies it.
    return !atom_constraint || (*atom_constraint)(atom);
  };

  bool IsAllowed(Tag bond_tag, const RDKit::Bond* bond) const {
    const BondConstraint* bond_constraint = GetBondConstraint(bond_tag);
    return !bond_constraint || (*bond_constraint)(bond);
  };

  bool IsAllowed(
    const RDKit::Atom* prior_atom,
    const RDKit::Atom* posterior_atom) const {
    auto [atom_tag, atom_constraint, updated] = UpdatedAtomConstraint(
      prior_atom, posterior_atom);
    return !atom_constraint || (*atom_constraint)(posterior_atom);
  };

  bool IsAllowed(
    const RDKit::Bond* prior_bond,
    const RDKit::Bond* posterior_bond) const {
    auto [bond_tag, bond_constraint, updated] = UpdatedBondConstraint(
      prior_bond, posterior_bond);
    return !bond_constraint || (*bond_constraint)(posterior_bond);
  };

  bool IsAllowed(
    const MolecularPerturbation& perturbation,
    const RDKit::ROMol& molecule) const {
    if (!(*this)) {
      return true;
    };
    RDKit::RWMol perturbed_molecule = perturbation(molecule);
    if (MustCheckAtomConstraints()) {
      auto atom_pairs = AtomPairs(molecule, perturbed_molecule);
      for (auto [prior_atom, posterior_atom] : atom_pairs) {
        if (!IsAllowed(prior_atom, posterior_atom)) {
          return false;
        };
      };
    };
    if (MustCheckBondConstraints()) {
      auto bond_pairs = BondPairs(molecule, perturbed_molecule);
      for (auto [prior_bond, posterior_bond] : bond_pairs) {
        if (!IsAllowed(prior_bond, posterior_bond)) {
          return false;
        };
      };
    };
    if (MustCheckMoleculeConstraints()) {
      for (const auto& [molecule_constraint, _] : molecule_constraints) {
        if (!molecule_constraint(&perturbed_molecule)) {
          return false;
        };
      };
    };
    return true;
  };

  bool UpdateIfAllowed(
    const RDKit::Atom* prior_atom,
    const RDKit::Atom* posterior_atom) {
    auto [atom_tag, atom_constraint, updated] = 
      UpdatedAtomConstraint(prior_atom, posterior_atom);
    if (!atom_constraint) {
      ClearAtomConstraint(atom_tag);
      return true;
    };
    if ((*atom_constraint)(posterior_atom)) {
      if (updated) {
        SetAtomConstraint(atom_tag, std::move(*atom_constraint));
      };
      return true;
    };
    return false;
  };

  bool UpdateIfAllowed(
    const RDKit::Bond* prior_bond,
    const RDKit::Bond* posterior_bond) {
    auto [bond_tag, bond_constraint, updated] = 
      UpdatedBondConstraint(prior_bond, posterior_bond);
    if (!bond_constraint) {
      ClearBondConstraint(bond_tag);
      return true;
    };
    if ((*bond_constraint)(posterior_bond)) {
      if (updated) {
        SetBondConstraint(bond_tag, std::move(*bond_constraint));
      };
      return true;
    };
    return false;
  };

  bool UpdateIfAllowed(
    const MolecularPerturbation& perturbation,
    const RDKit::ROMol& molecule) {
    if (!HasConstraintGenerators()) {
      return IsAllowed(perturbation, molecule);
    };
    std::vector<std::pair<Tag, std::shared_ptr<const AtomConstraint>>>
      updated_atom_constraints;
    std::vector<std::pair<Tag, std::shared_ptr<const BondConstraint>>>
      updated_bond_constraints;
    RDKit::RWMol perturbed_molecule = perturbation(molecule);
    if (MustCheckAtomConstraints()) {
      auto atom_pairs = AtomPairs(molecule, perturbed_molecule);
      for (auto [prior_atom, posterior_atom] : atom_pairs) {
        auto [atom_tag, atom_constraint, updated] = 
          UpdatedAtomConstraint(prior_atom, posterior_atom);
        if (atom_constraint && !(*atom_constraint)(posterior_atom)) {
          return false;
        };
        if (updated) {
          updated_atom_constraints.emplace_back(
            atom_tag, std::move(atom_constraint));
        };
      };       
    };
    if (MustCheckBondConstraints()) {
      auto bond_pairs = BondPairs(molecule, perturbed_molecule);
      for (auto [prior_bond, posterior_bond] : bond_pairs) {
        auto [bond_tag, bond_constraint, updated] = 
          UpdatedBondConstraint(prior_bond, posterior_bond);
        if (bond_constraint && !(*bond_constraint)(posterior_bond)) {
          return false;
        };
        if (updated) {
          updated_bond_constraints.emplace_back(
            bond_tag, std::move(bond_constraint));
        };
      };       
    };
    if (MustCheckMoleculeConstraints()) {
      for (const auto& [molecule_constraint, _] : molecule_constraints) {
        if (!molecule_constraint(&perturbed_molecule)) {
          return false;
        };
      };
    };
    for (const auto& [atom_tag, atom_constraint] : updated_atom_constraints) {
      if (atom_constraint) {
        SetAtomConstraint(atom_tag, std::move(*atom_constraint));
      } else {
        ClearAtomConstraint(atom_tag);
      };
    };
    for (const auto& [bond_tag, bond_constraint] : updated_bond_constraints) {
      if (bond_constraint) {
        SetBondConstraint(bond_tag, std::move(*bond_constraint));
      } else {
        ClearBondConstraint(bond_tag);
      };
    };
    return true;
  };

  bool HasAtomConstraintGenerator() const {
    return !!atom_constraint_generator;
  };

  bool HasBondConstraintGenerator() const {
    return !!bond_constraint_generator;
  };

  bool HasConstraintGenerators() const {
    return atom_constraint_generator || 
      bond_constraint_generator;
  };

  bool MustCheckAtomConstraints() const {
    return !atom_constraints.empty() || atom_constraint_generator;
  };

  bool MustCheckBondConstraints() const {
    return !bond_constraints.empty() || bond_constraint_generator;
  };

  bool MustCheckMoleculeConstraints() const {
    return !molecule_constraints.empty();
  };

  std::size_t Size() const {
    return atom_constraints.size() + 
      bond_constraints.size() + 
      molecule_constraints.size();
  };

  explicit operator bool() const {
    return MustCheckAtomConstraints() ||
      MustCheckBondConstraints() ||
      MustCheckMoleculeConstraints();
  };
};

#endif // !_MOLECULAR_CONSTRAINTS_HPP_
