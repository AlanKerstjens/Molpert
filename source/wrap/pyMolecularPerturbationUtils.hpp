#pragma once
#ifndef _PY_MOLECULAR_PERTURBATION_UTILS_HPP_
#define _PY_MOLECULAR_PERTURBATION_UTILS_HPP_

#include "MolecularPerturbationUtils.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrappedCorrectAromaticity(RDKit::ROMol& molecule) {
  RDKit::RWMol& rwmol = static_cast<RDKit::RWMol&> (molecule);
  CorrectAromaticity(rwmol);
};

void WrappedCorrectHydrogenCounts(RDKit::ROMol& molecule) {
  RDKit::RWMol& rwmol = static_cast<RDKit::RWMol&> (molecule);
  CorrectHydrogenCounts(rwmol);
};

void WrappedPartialSanitization(RDKit::ROMol& molecule) {
  RDKit::RWMol& rwmol = static_cast<RDKit::RWMol&> (molecule);
  PartialSanitization(rwmol);
};

void WrapMolecularPerturbationUtils() {
  python::def("CorrectAromaticity", WrappedCorrectAromaticity, (
    python::arg("molecule")));

  python::def("CorrectHydrogenCounts", WrappedCorrectHydrogenCounts, (
    python::arg("molecule")));

  python::def("PartialSanitization", WrappedPartialSanitization, (
    python::arg("molecule")));
};

#endif // !_PY_MOLECULAR_PERTURBATION_UTILS_HPP_