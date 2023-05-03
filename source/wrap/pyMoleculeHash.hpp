#pragma once
#ifndef _PY_MOLECULE_HASH_HPP_
#define _PY_MOLECULE_HASH_HPP_

#include "MoleculeHash.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapMoleculeHash() {

  python::def("LoneAtomHash", LoneAtomHash, (python::arg("atom")));
  python::def("RingAwareAtomHash", RingAwareAtomHash, (python::arg("atom")));
  python::def("AtomKeyHash", AtomKeyHash, (python::arg("atom")));

  python::def("HashMolecule", HashMolecule, (python::arg("molecule")));

};

#endif // !_PY_MOLECULE_HASH_HPP_