#ifndef _PY_MOLECULAR_KEYS_HPP_
#define _PY_MOLECULAR_KEYS_HPP_

#include "MolecularKeys.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapMolecularKeys() {
  python::class_<AtomKey>(
    "AtomKey",
    python::init<std::uint8_t, std::int8_t, std::uint8_t, std::int8_t, std::uint8_t>((
      python::arg("degree"), 
      python::arg("valence"),
      python::arg("atomic_number"),
      python::arg("formal_charge"),
      python::arg("n_hydrogens"))))
    .def(python::init<const RDKit::Atom*>((python::arg("atom"))))
    .def_readwrite("degree", &AtomKey::degree)
    .def_readwrite("valence", &AtomKey::valence)
    .def_readwrite("atomic_number", &AtomKey::atomic_number)
    .def_readwrite("formal_charge", &AtomKey::formal_charge)
    .def_readwrite("n_hydrogens", &AtomKey::n_hydrogens)
    .def("__str__", &AtomKey::str)
    .def("__repr__", &AtomKey::str);

  python::class_<BondKey>(
    "BondKey", 
    python::init<const AtomKey&, const AtomKey&, RDKit::Bond::BondType, bool>((
      python::arg("atom_key1"), 
      python::arg("atom_key2"),
      python::arg("bond_type"),
      python::arg("canonicalize") = true)))
    .def(python::init<const RDKit::Atom*, const RDKit::Atom*, RDKit::Bond::BondType, bool>((
      python::arg("atom1"),
      python::arg("atom2"),
      python::arg("bond_type"),
      python::arg("canonicalize") = true)))
    .def(python::init<const RDKit::Bond*, bool>((
      python::arg("bond"),
      python::arg("canonicalize") = true)))
    .def_readwrite("atom_key1", &BondKey::atom_key1)
    .def_readwrite("atom_key2", &BondKey::atom_key2)
    .def_readwrite("bond_type", &BondKey::bond_type)
    .def("__str__", &BondKey::str)
    .def("__repr__", &BondKey::str);
};

#endif // !_PY_MOLECULAR_KEYS_HPP_