#pragma once
#ifndef _PY_MOLECULAR_PERTURBATIONS_HPP_
#define _PY_MOLECULAR_PERTURBATIONS_HPP_

#include "MolecularPerturbations.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void ExecutePerturbation(
  const MolecularPerturbation& perturbation, RDKit::ROMol& molecule) {
  RDKit::RWMol& rwmol = static_cast<RDKit::RWMol&> (molecule);
  perturbation(rwmol);
};

python::object PythonPerturbation(const MolecularPerturbation* perturbation) {
  switch (perturbation->GetType()) {
    case MolecularPerturbation::Type::AtomicNumberChange_t:
      return python::object(
        static_cast<const AtomicNumberChange*>(perturbation));
    case MolecularPerturbation::Type::FormalChargeChange_t:
      return python::object(
        static_cast<const FormalChargeChange*>(perturbation));
    case MolecularPerturbation::Type::ExplicitHydrogensChange_t:
      return python::object(
        static_cast<const ExplicitHydrogensChange*>(perturbation));
    case MolecularPerturbation::Type::BondTypeChange_t:
      return python::object(
        static_cast<const BondTypeChange*>(perturbation));
    case MolecularPerturbation::Type::AtomInsertion_t:
      return python::object(
        static_cast<const AtomInsertion*>(perturbation));
    case MolecularPerturbation::Type::AtomDeletion_t:
      return python::object(
        static_cast<const AtomDeletion*>(perturbation));
    case MolecularPerturbation::Type::BondConstruction_t:
    case MolecularPerturbation::Type::BondInsertion_t:
      return python::object(
        static_cast<const BondConstruction*>(perturbation));
    case MolecularPerturbation::Type::BondDeletion_t:
      return python::object(
        static_cast<const BondDeletion*>(perturbation));
    case MolecularPerturbation::Type::SubgraphConstruction_t:
      return python::object(
        static_cast<const SubgraphConstruction*>(perturbation));
    case MolecularPerturbation::Type::SubgraphDestruction_t:
      return python::object(
        static_cast<const SubgraphDestruction*>(perturbation));
    case MolecularPerturbation::Type::SubgraphPerturbation_t:
      return python::object(
        static_cast<const SubgraphPerturbation*>(perturbation));
    default:
      return python::object();
  };
};

template<class T>
std::vector<T> to_vector(const python::object& iterable) {
  return std::vector<T>(
    python::stl_input_iterator<T>(iterable), python::stl_input_iterator<T>());
};

std::shared_ptr<AtomInsertion> AtomInsertionFactory(
  const RDKit::ROMol& molecule,
  const python::object& neighbor_atom_indices,
  const python::object& bond_types,
  std::uint8_t atomic_number = 6,
  std::int8_t formal_charge = 0,
  std::uint8_t n_explicit_hydrogens = 0,
  int dropped_atom_idx = -1) {
  return std::shared_ptr<AtomInsertion>(new AtomInsertion(
    molecule,
    to_vector<AtomIdx>(neighbor_atom_indices),
    to_vector<RDKit::Bond::BondType>(bond_types),
    atomic_number, formal_charge, n_explicit_hydrogens, dropped_atom_idx));
};

std::shared_ptr<AtomDeletion> AtomDeletionFactory(
  const RDKit::ROMol& molecule,
  AtomIdx atom_idx,
  AtomIdx reconnection_center_atom_idx,
  const python::object& reconnection_partner_atom_indices,
  const python::object& reconnection_bond_types) {
  return std::shared_ptr<AtomDeletion>(new AtomDeletion(
    molecule, atom_idx, reconnection_center_atom_idx,
    to_vector<AtomIdx>(reconnection_partner_atom_indices),
    to_vector<RDKit::Bond::BondType>(reconnection_bond_types)));
};

std::shared_ptr<SubgraphConstruction> SubgraphConstructionFactory(
  const RDKit::ROMol& target,
  const RDKit::ROMol& source,
  const python::object& subgraph_atom_indices) {
  return std::shared_ptr<SubgraphConstruction>(new SubgraphConstruction(
    target, source, to_vector<AtomIdx>(subgraph_atom_indices)));
};

std::shared_ptr<SubgraphDestruction> SubgraphDestructionFactory(
  const RDKit::ROMol& molecule,
  const python::object& subgraph_atom_indices) {
  return std::shared_ptr<SubgraphDestruction>(new SubgraphDestruction(
    molecule, to_vector<AtomIdx>(subgraph_atom_indices)));
};

std::shared_ptr<SubgraphPerturbation> SubgraphPerturbationFactory(
  const RDKit::ROMol& target,
  const RDKit::ROMol& source,
  const python::object& constructed_subgraph_atom_indices,
  const python::object& destroyed_subgraph_atom_indices) {
  return std::shared_ptr<SubgraphPerturbation>(new SubgraphPerturbation(
    target, source, 
    to_vector<AtomIdx>(constructed_subgraph_atom_indices),
    to_vector<AtomIdx>(destroyed_subgraph_atom_indices)));
};

void WrapMolecularPerturbations() {
  python::scope module_scope = python::scope();

  python::def("ExecutePerturbation", ExecutePerturbation, (
    python::arg("perturbation"),
    python::arg("molecule")));

  // The base class is exposed to enable polymorphism in Python. For the time
  // being inheriting from it on the Python side isn't supported because things
  // get messy. But nothing is more permanent than a temporary solution.
  python::scope base_class_scope =
  python::class_<MolecularPerturbation, boost::noncopyable>(
    "MolecularPerturbation", python::no_init);

  python::enum_<MolecularPerturbation::Type>("Type")
    .value("AtomicNumberChange", MolecularPerturbation::AtomicNumberChange_t)
    .value("FormalChargeChange", MolecularPerturbation::FormalChargeChange_t)
    .value("ExplicitHydrogensChange", MolecularPerturbation::ExplicitHydrogensChange_t)
    .value("BondTypeChange", MolecularPerturbation::BondTypeChange_t)
    .value("AtomConstruction", MolecularPerturbation::AtomConstruction_t)
    .value("AtomDestruction", MolecularPerturbation::AtomDestruction_t)
    .value("BondConstruction", MolecularPerturbation::BondConstruction_t)
    .value("BondDestruction", MolecularPerturbation::BondDestruction_t)
    .value("TopologicalPerturbation", MolecularPerturbation::TopologicalPerturbation_t)
    .value("AtomInsertion", MolecularPerturbation::AtomInsertion_t)
    .value("AtomDeletion", MolecularPerturbation::AtomDeletion_t)
    .value("BondInsertion", MolecularPerturbation::BondInsertion_t)
    .value("BondDeletion", MolecularPerturbation::BondDeletion_t)
    .value("SubgraphConstruction", MolecularPerturbation::SubgraphConstruction_t)
    .value("SubgraphDestruction", MolecularPerturbation::SubgraphDestruction_t)
    .value("SubgraphPerturbation", MolecularPerturbation::SubgraphPerturbation_t);

  python::scope derived_class_scope (module_scope);

  python::class_<AtomicNumberChange, python::bases<MolecularPerturbation>>(
    "AtomicNumberChange",
    python::init<AtomIdx, std::uint8_t>((
      python::arg("atom_idx"),
      python::arg("atomic_number"))))
    .def("ID", &AtomicNumberChange::ID)
    .add_property("type", &AtomicNumberChange::GetType);

  python::class_<FormalChargeChange, python::bases<MolecularPerturbation>>(
    "FormalChargeChange",
    python::init<AtomIdx, std::int8_t>((
      python::arg("atom_idx"),
      python::arg("formal_charge"))))
    .def("ID", &FormalChargeChange::ID)
    .add_property("type", &FormalChargeChange::GetType);

  python::class_<ExplicitHydrogensChange, python::bases<MolecularPerturbation>>(
    "ExplicitHydrogensChange",
    python::init<AtomIdx, std::uint8_t>((
      python::arg("atom_idx"),
      python::arg("n_explicit_hydrogens"))))
    .def("ID", &ExplicitHydrogensChange::ID)
    .add_property("type", &ExplicitHydrogensChange::GetType);

  python::class_<BondTypeChange, python::bases<MolecularPerturbation>>(
    "BondTypeChange",
    python::init<AtomIdx, AtomIdx, RDKit::Bond::BondType>((
      python::arg("begin_atom_idx"),
      python::arg("end_atom_idx"),
      python::arg("bond_type"))))
    .def(python::init<const RDKit::ROMol&, BondIdx, RDKit::Bond::BondType>((
      python::arg("molecule"),
      python::arg("bond_idx"),
      python::arg("bond_type"))))
    .def("ID", &BondTypeChange::ID)
    .add_property("type", &BondTypeChange::GetType);

  python::class_<AtomInsertion, std::shared_ptr<AtomInsertion>, python::bases<MolecularPerturbation>>(
    "AtomInsertion",
    python::init<const RDKit::ROMol&, std::uint8_t, std::int8_t, std::uint8_t>((
      python::arg("molecule"),
      python::arg("atomic_number") = 6,
      python::arg("formal_charge") = 0,
      python::arg("n_explicit_hydrogens") = 0)))
    .def("__init__", python::make_constructor(&AtomInsertionFactory,
      python::default_call_policies(), (
      python::arg("molecule"),
      python::arg("neighbor_atom_indices"),
      python::arg("bond_types"),
      python::arg("atomic_number") = 6,
      python::arg("formal_charge") = 0,
      python::arg("n_explicit_hydrogens") = 0,
      python::arg("dropped_atom_idx") = -1)))
    .def("ID", &AtomInsertion::ID)
    .add_property("type", &AtomInsertion::GetType);

  python::class_<AtomDeletion, std::shared_ptr<AtomDeletion>, python::bases<MolecularPerturbation>>(
    "AtomDeletion",
    python::init<AtomIdx>(
      python::arg("atom_idx")))
    .def(python::init<const RDKit::ROMol&, AtomIdx, AtomIdx, bool, RDKit::Bond::BondType>((
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("reconnection_atom_idx"),
      python::arg("preserve_bond_types") = true,
      python::arg("default_bond_type") = 1)))
    .def("__init__", python::make_constructor(&AtomDeletionFactory,
      python::default_call_policies(), (
      python::arg("molecule"),
      python::arg("atom_idx"),
      python::arg("reconnection_center_atom_idx"),
      python::arg("reconnection_partner_atom_indices"),
      python::arg("reconnection_bond_types"))))
    .def("ID", &AtomDeletion::ID)
    .add_property("type", &AtomDeletion::GetType);

  python::class_<BondInsertion, python::bases<MolecularPerturbation>>(
    "BondInsertion",
    python::init<const RDKit::ROMol&, AtomIdx, AtomIdx, RDKit::Bond::BondType>((
      python::arg("molecule"),
      python::arg("begin_atom_idx"),
      python::arg("end_atom_idx"),
      python::arg("bond_type") = 1)))
    .def("ID", &BondInsertion::ID)
    .add_property("type", &BondInsertion::GetType);

  python::class_<BondDeletion, python::bases<MolecularPerturbation>>(
    "BondDeletion",
    python::init<AtomIdx, AtomIdx>((
      python::arg("begin_atom_idx"),
      python::arg("end_atom_idx"))))
    .def(python::init<const RDKit::ROMol&, BondIdx>((
      python::arg("molecule"),
      python::arg("bond_idx"))))
    .def(python::init<const RDKit::ROMol&, AtomIdx, AtomIdx, AtomIdx, AtomIdx>((
      python::arg("molecule"),
      python::arg("begin_atom_idx"),
      python::arg("end_atom_idx"),
      python::arg("reroute_begin_atom_idx"),
      python::arg("reroute_end_atom_idx"))))
    .def("ID", &BondDeletion::ID)
    .add_property("type", &BondDeletion::GetType);

  python::class_<SubgraphConstruction, std::shared_ptr<SubgraphConstruction>, python::bases<MolecularPerturbation>>(
    "SubgraphConstruction", python::no_init)
    .def("__init__", python::make_constructor(&SubgraphConstructionFactory,
      python::default_call_policies(), (
      python::arg("target_molecule"),
      python::arg("source_molecule"),
      python::arg("subgraph_atom_indices"))))
    .def("ID", &SubgraphConstruction::ID)
    .add_property("type", &SubgraphConstruction::GetType);

  python::class_<SubgraphDestruction, std::shared_ptr<SubgraphDestruction>, python::bases<MolecularPerturbation>>(
    "SubgraphDestruction", python::no_init)
    .def("__init__", python::make_constructor(&SubgraphDestructionFactory,
      python::default_call_policies(), (
      python::arg("molecule"),
      python::arg("subgraph_atom_indices"))))
    .def("ID", &SubgraphDestruction::ID)
    .add_property("type", &SubgraphDestruction::GetType);

  python::class_<SubgraphPerturbation, std::shared_ptr<SubgraphPerturbation>, python::bases<MolecularPerturbation>>(
    "SubgraphPerturbation", python::no_init)
    .def("__init__", python::make_constructor(&SubgraphPerturbationFactory,
      python::default_call_policies(), (
      python::arg("target_molecule"),
      python::arg("source_molecule"),
      python::arg("constructed_subgraph_atom_indices"),
      python::arg("destroyed_subgraph_atom_indices"))))
    .def("ID", &SubgraphPerturbation::ID)
    .add_property("type", &SubgraphPerturbation::GetType);
};

#endif // !_PY_MOLECULAR_PERTURBATIONS_HPP_
