#pragma once
#ifndef _PY_MOLECULAR_CONSTRAINTS_HPP_
#define _PY_MOLECULAR_CONSTRAINTS_HPP_

#include "MolecularConstraints.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

// Credit for the converter voodoo goes to Tanner Sansbury @ StackOverflow.
struct FunctionConverter {
  // Registers converter from a python callable type to the provided type.
  template <class FunctionSignature>
  FunctionConverter& Register() {
    python::converter::registry::push_back(
      &FunctionConverter::Convertible,
      &FunctionConverter::Construct<FunctionSignature>,
      python::type_id<std::function<FunctionSignature>>());
    return *this; // Support chaining.
  };

  // Check if PyObject is callable.
  static void* Convertible(PyObject* object) {
    return PyCallable_Check(object) ? object : NULL;
  };

  // Convert callable PyObject to a C++ std::function.
  template <class FunctionSignature>
  static void Construct(
    PyObject* object,
    python::converter::rvalue_from_python_stage1_data* data) {
    // Object is a borrowed reference, so create a handle indicting it is
    // borrowed for proper reference counting.
    python::handle<> handle (python::borrowed(object));
    // Obtain a handle to the memory block that the converter has allocated
    // for the C++ type.
    typedef std::function<FunctionSignature> Function;
    typedef python::converter::rvalue_from_python_storage<Function> Storage;
    void* storage = reinterpret_cast<Storage*>(data)->storage.bytes;
    // Allocate the C++ type into the converter's memory block, and assign
    // its handle to the converter's convertible variable.
    new (storage) Function(python::object(handle));
    data->convertible = storage;
  };
};

MolecularConstraints::AtomConstraintGenerator AtomConstraintGeneratorFactory(
  MolecularConstraints::AtomConstraintType atom_constraint_type = 
    MolecularConstraints::AtomConstraintType::Null,
  const ChemicalDictionary* dictionary = nullptr) {
  switch (atom_constraint_type) {
    case MolecularConstraints::AtomConstraintType::Null:
      return nullptr;
    case MolecularConstraints::AtomConstraintType::Valence:
      return ValenceConstraintGenerator;
    case MolecularConstraints::AtomConstraintType::AtomKey:
      if (dictionary) {
        return AtomKeyConstraintGenerator(dictionary);
      };
    default:
      return nullptr;
  };
};

MolecularConstraints::BondConstraintGenerator BondConstraintGeneratorFactory(
  MolecularConstraints::BondConstraintType bond_constraint_type = 
    MolecularConstraints::BondConstraintType::Null,
  const ChemicalDictionary* dictionary = nullptr) {
  switch (bond_constraint_type) {
    case MolecularConstraints::BondConstraintType::Null:
      return nullptr;
    case MolecularConstraints::BondConstraintType::BondKey:
      if (dictionary) {
        return BondKeyConstraintGenerator(dictionary);
      };
    default:
      return nullptr;
  };
};

MolecularConstraints::EnvironmentConstraintGenerator
EnvironmentConstraintGeneratorFactory(
  MolecularConstraints::EnvironmentConstraintType environment_constraint_type = 
    MolecularConstraints::EnvironmentConstraintType::Null,
  const ChemicalDictionary* dictionary = nullptr) {
  switch (environment_constraint_type) {
    case MolecularConstraints::EnvironmentConstraintType::Null:
      return nullptr;
    case MolecularConstraints::EnvironmentConstraintType::EnvironmentKey:
      if (dictionary) {
        return EnvironmentKeyConstraintGenerator(dictionary);
      };
    default:
      return nullptr;
  };
};

std::shared_ptr<MolecularConstraints> MolecularConstraintsFactory(
  MolecularConstraints::AtomConstraintType atom_constraint_type = 
    MolecularConstraints::AtomConstraintType::Null,
  MolecularConstraints::BondConstraintType bond_constraint_type = 
    MolecularConstraints::BondConstraintType::Null,
  MolecularConstraints::EnvironmentConstraintType environment_constraint_type = 
    MolecularConstraints::EnvironmentConstraintType::Null,
  const ChemicalDictionary* dictionary = nullptr) {
  auto atom_constraints_generator = AtomConstraintGeneratorFactory(
    atom_constraint_type, dictionary);
  auto bond_constraint_generator = BondConstraintGeneratorFactory(
    bond_constraint_type, dictionary);
  auto environment_constraint_generator = EnvironmentConstraintGeneratorFactory(
    environment_constraint_type, dictionary);
  unsigned environment_radius = 
    dictionary ? dictionary->GetEnvironmentRadius() : 2;
  return std::shared_ptr<MolecularConstraints>(new MolecularConstraints(
    atom_constraints_generator,
    bond_constraint_generator,
    environment_constraint_generator,
    environment_radius));
};

void WrapMolecularConstraints() {
  FunctionConverter()
    .Register<bool(const AtomKey&)>()
    .Register<bool(const BondKey&)>()
    .Register<bool(const EnvironmentKey&)>();

  python::class_ constraints = python::class_<
    MolecularConstraints, std::shared_ptr<MolecularConstraints>>(
    "MolecularConstraints", python::init());

  python::scope constraints_scope (constraints);

  python::enum_<MolecularConstraints::AtomConstraintType>("AtomConstraintType")
    .value("Null", MolecularConstraints::AtomConstraintType::Null)
    .value("Valence", MolecularConstraints::AtomConstraintType::Valence)
    .value("AtomKey", MolecularConstraints::AtomConstraintType::AtomKey);

  python::enum_<MolecularConstraints::BondConstraintType>("BondConstraintType")
    .value("Null", MolecularConstraints::BondConstraintType::Null)
    .value("BondKey", MolecularConstraints::BondConstraintType::BondKey);

  python::enum_<MolecularConstraints::EnvironmentConstraintType>(
    "EnvironmentConstraintType")
    .value("Null", MolecularConstraints::EnvironmentConstraintType::Null)
    .value("EnvironmentKey",
      MolecularConstraints::EnvironmentConstraintType::EnvironmentKey);

  constraints.def("__init__", python::make_constructor(&MolecularConstraintsFactory, 
    python::default_call_policies(), (
    python::arg("atom_constraint_type") = MolecularConstraints::AtomConstraintType::Null, 
    python::arg("bond_constraint_type") = MolecularConstraints::BondConstraintType::Null,
    python::arg("environment_constraint_type") = MolecularConstraints::EnvironmentConstraintType::Null,
    python::arg("dictionary") = python::ptr((const ChemicalDictionary*) nullptr))))
  .def("GenerateAtomConstraint", &MolecularConstraints::GenerateAtomConstraint, (
    python::arg("atom")))
  .def("GenerateAtomConstraints", &MolecularConstraints::GenerateAtomConstraints, (
    python::arg("molecule")))
  .def("GenerateBondConstraint", &MolecularConstraints::GenerateBondConstraint, (
    python::arg("bond")))
  .def("GenerateBondConstraints", &MolecularConstraints::GenerateBondConstraints, (
    python::arg("molecule")))
  .def("GenerateEnvironmentConstraint", &MolecularConstraints::GenerateEnvironmentConstraint, (
    python::arg("molecule"),
    python::arg("atom"),
    python::arg("recalculate_atom_hashes") = false))
  .def("GenerateEnvironmentConstraints", &MolecularConstraints::GenerateEnvironmentConstraints, (
    python::arg("molecule")))
  .def("GenerateConstraints", &MolecularConstraints::GenerateConstraints, (
    python::arg("molecule")))
  .def("SetAtomConstraint", &MolecularConstraints::SetAtomConstraint, (
    python::arg("atom_tag"), 
    python::arg("atom_constraint"), 
    python::arg("replace") = true,
    python::arg("make_static") = false))
  .def("SetBondConstraint", &MolecularConstraints::SetBondConstraint, (
    python::arg("bond_tag"), 
    python::arg("bond_constraint"), 
    python::arg("replace") = true,
    python::arg("make_static") = false))
  .def("SetEnvironmentConstraint", &MolecularConstraints::SetEnvironmentConstraint, (
    python::arg("atom_tag"), 
    python::arg("environment_constraint"), 
    python::arg("replace") = true,
    python::arg("make_static") = false))
  .def("UpdateConstraints", &MolecularConstraints::UpdateConstraints, (
    python::arg("molecule"),
    python::arg("perturbation")))
  .def("ClearAtomConstraint", &MolecularConstraints::ClearAtomConstraint, (
    python::arg("atom_tag"),
    python::arg("clear_static") = false))
  .def("ClearBondConstraint", &MolecularConstraints::ClearBondConstraint, (
    python::arg("bond_tag"),
    python::arg("clear_static") = false))
  .def("ClearEnvironmentConstraint", &MolecularConstraints::ClearEnvironmentConstraint, (
    python::arg("atom_tag"),
    python::arg("clear_static") = false))
  .def("ClearAtomConstraints", &MolecularConstraints::ClearAtomConstraints, (
    python::arg("clear_static") = false))
  .def("ClearBondConstraints", &MolecularConstraints::ClearBondConstraints, (
    python::arg("clear_static") = false))
  .def("ClearEnvironmentConstraints", &MolecularConstraints::ClearEnvironmentConstraints, (
    python::arg("clear_static") = false))
  .def("ClearCyclicityConstraints", &MolecularConstraints::ClearCyclicityConstraints)
  .def("ClearConstraints", &MolecularConstraints::ClearConstraints, (
    python::arg("clear_static") = false,
    python::arg("clear_cyclicity_constraints") = true))
  .def("Clear", &MolecularConstraints::Clear)
  .def<bool (MolecularConstraints::*)(Tag, const AtomKey&) const>(
    "IsAtomKeyAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("atom_tag"),
    python::arg("atom_key")))
  .def<bool (MolecularConstraints::*)(Tag, const BondKey&) const>(
    "IsBondKeyAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("bond_tag"),
    python::arg("bond_key")))
  .def<bool (MolecularConstraints::*)(Tag, const EnvironmentKey&) const>(
    "IsEnvironmentKeyAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("atom_tag"),
    python::arg("environment_key")))
  .def<bool (MolecularConstraints::*)(Tag, const AtomKeyChange&) const>(
    "IsAtomKeyChangeAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("atom_tag"),
    python::arg("atom_key_change")))
  .def<bool (MolecularConstraints::*)(Tag, const BondKeyChange&) const>(
    "IsBondKeyChangeAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("bond_tag"),
    python::arg("bond_key_change")))
  .def<bool (MolecularConstraints::*)(Tag, const EnvironmentKeyChange&) const>(
    "IsEnvironmentKeyChangeAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("atom_tag"),
    python::arg("environment_key_change")))
  .def<bool (MolecularConstraints::*)(const RDKit::ROMol&, const MolecularPerturbation&) const>(
    "IsPerturbationAllowed", &MolecularConstraints::IsAllowed, (
    python::arg("molecule"),
    python::arg("perturbation")))
  .def<bool (MolecularConstraints::*)(const RDKit::ROMol&, const MolecularPerturbation&)>(
    "UpdateConstraintsIfAllowed", &MolecularConstraints::UpdateIfAllowed, (
    python::arg("molecule"),
    python::arg("perturbation")))
  .def("HasAtomConstraintGenerator", &MolecularConstraints::HasAtomConstraintGenerator)
  .def("HasBondConstraintGenerator", &MolecularConstraints::HasBondConstraintGenerator)
  .def("HasEnvironmentConstraintGenerator", &MolecularConstraints::HasEnvironmentConstraintGenerator)
  .def("HasConstraintGenerators", &MolecularConstraints::HasConstraintGenerators)
  .def("HasCyclicityConstraints", &MolecularConstraints::HasCyclicityConstraints)
  .def("__bool__", &MolecularConstraints::operator bool)
  .def("__len__", &MolecularConstraints::Size)
  .add_property("min_n_cycles",
    &MolecularConstraints::GetMinNCycles, 
    &MolecularConstraints::SetMinNCycles)
  .add_property("max_n_cycles",
    &MolecularConstraints::GetMaxNCycles, 
    &MolecularConstraints::SetMaxNCycles)
  .add_property("min_cycle_size",
    &MolecularConstraints::GetMinCycleSize, 
    &MolecularConstraints::SetMinCycleSize)
  .add_property("max_cycle_size",
    &MolecularConstraints::GetMaxCycleSize, 
    &MolecularConstraints::SetMaxCycleSize)
  .add_property("max_atom_cycles_membership",
    &MolecularConstraints::GetMaxAtomCyclesMembership, 
    &MolecularConstraints::SetMaxAtomCyclesMembership)
  .add_property("environment_radius", &MolecularConstraints::GetEnvironmentRadius);
};

#endif // !_PY_MOLECULAR_CONSTRAINTS_HPP_