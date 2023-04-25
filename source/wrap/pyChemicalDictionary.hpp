#pragma once
#ifndef _PY_CHEMICAL_DICTIONARY_HPP_
#define _PY_CHEMICAL_DICTIONARY_HPP_

#include "ChemicalDictionary.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapChemicalDictionary() {
  python::class_<ChemicalDictionary>(
    "ChemicalDictionary",
    python::init<unsigned>((
      python::arg("environment_radius") = 2)))
    .def(python::init<const std::string&>((
      python::arg("path"))))
    .def("AddMolecule", &ChemicalDictionary::AddMolecule, (
      python::arg("molecule")))
    .def("BuildPartialKeyDictionaries", &ChemicalDictionary::BuildPartialKeyDictionaries)
    .def("AtomFrequency", &ChemicalDictionary::AtomFrequency, (
      python::arg("atom_key")))
    .def("BondFrequency", &ChemicalDictionary::BondFrequency, (
      python::arg("bond_key")))
    .def("EnvironmentFrequency", &ChemicalDictionary::EnvironmentFrequency, (
      python::arg("environment_key")))
    .def<bool (ChemicalDictionary::*)(const AtomKey&) const>(
      "IsForeignAtom", &ChemicalDictionary::IsForeignAtom, (
        python::arg("atom_key")))
    .def<bool (ChemicalDictionary::*)(const RDKit::Atom*) const>(
      "IsForeignAtom", &ChemicalDictionary::IsForeignAtom, (
        python::arg("atom")))
    .def<bool (ChemicalDictionary::*)(const BondKey&) const>(
      "IsForeignBond", &ChemicalDictionary::IsForeignBond, (
        python::arg("bond_key")))
    .def<bool (ChemicalDictionary::*)(const RDKit::Bond*) const>(
      "IsForeignBond", &ChemicalDictionary::IsForeignBond, (
        python::arg("bond")))
    .def<bool (ChemicalDictionary::*)(const EnvironmentKey&) const>(
      "IsForeignEnvironment", &ChemicalDictionary::IsForeignEnvironment, (
        python::arg("environment_key")))
    .def<bool (ChemicalDictionary::*)(const RDKit::Atom*)>(
      "IsForeignEnvironment", &ChemicalDictionary::IsForeignEnvironment, (
        python::arg("atom")))
    .def("Save", &ChemicalDictionary::Save, (
      python::arg("path")))
    .def("Load", &ChemicalDictionary::Load, (
      python::arg("path")))
    .add_property("total_atom_frequency", 
      &ChemicalDictionary::GetTotalAtomFrequency)
    .add_property("total_bond_frequency", 
      &ChemicalDictionary::GetTotalBondFrequency)
    .add_property("environment_radius", 
      &ChemicalDictionary::GetEnvironmentRadius)
    .def_readwrite("foreign_d_frequency_threshold", &ChemicalDictionary::foreign_d_frequency_threshold)
    .def_readwrite("foreign_dv_frequency_threshold", &ChemicalDictionary::foreign_dv_frequency_threshold)
    .def_readwrite("foreign_dvz_frequency_threshold", &ChemicalDictionary::foreign_dvz_frequency_threshold)
    .def_readwrite("foreign_dvzq_frequency_threshold", &ChemicalDictionary::foreign_dvzq_frequency_threshold)
    .def_readwrite("foreign_k1k2_frequency_threshold", &ChemicalDictionary::foreign_k1k2_frequency_threshold)
    .def_readwrite("foreign_atom_frequency_threshold", &ChemicalDictionary::foreign_atom_frequency_threshold)
    .def_readwrite("foreign_bond_frequency_threshold", &ChemicalDictionary::foreign_bond_frequency_threshold)
    .def_readwrite("foreign_environment_frequency_threshold", &ChemicalDictionary::foreign_environment_frequency_threshold);
};

#endif //_PY_CHEMICAL_DICTIONARY_HPP_