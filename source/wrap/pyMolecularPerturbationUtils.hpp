#pragma once
#ifndef _PY_MOLECULAR_PERTURBATION_UTILS_HPP_
#define _PY_MOLECULAR_PERTURBATION_UTILS_HPP_

#include "MolecularPerturbationUtils.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrappedPartialSanitization(
  RDKit::ROMol& molecule,
  bool kekulize = false,
  bool aromatize = true) {
  RDKit::RWMol& rwmol = static_cast<RDKit::RWMol&> (molecule);
  PartialSanitization(rwmol, kekulize, aromatize);
};

void WrappedCorrectElementAromaticity(RDKit::ROMol& molecule) {
  RDKit::RWMol& rwmol = static_cast<RDKit::RWMol&> (molecule);
  CorrectElementAromaticity(rwmol);
};

void WrapMolecularPerturbationUtils() {
  python::def("PartialSanitization", WrappedPartialSanitization, (
    python::arg("molecule"), 
    python::arg("kekulize") = false,
    python::arg("aromatize") = true));

  python::def("CorrectElementAromaticity", WrappedCorrectElementAromaticity, (
    python::arg("molecule")));
};

#endif // !_PY_MOLECULAR_PERTURBATION_UTILS_HPP_