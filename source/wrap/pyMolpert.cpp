#include "pyMolecularKeys.hpp"
#include "pyCircularAtomicEnvironment.hpp"
#include "pyMolecularPerturbations.hpp"
#include "pyMolecularPerturbationUtils.hpp"
#include "pyChemicalDictionary.hpp"
#include "pyMolecularConstraints.hpp"
#include "pyMoleculePerturber.hpp"
#include "pyAddons.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

BOOST_PYTHON_MODULE(molpert) {
  python::scope module_scope = python::scope();

  python::def("TagMolecule", TagMolecule, (
    python::arg("molecule"),
    python::arg("skip_if_tagged") = false,
    python::arg("start_atom_tag") = 0,
    python::arg("start_bond_tag") = 0));
  python::def("GetAtomTag", GetAtomTag, (
    python::arg("molecule"),
    python::arg("atom_idx")));
  python::def("GetBondTag", GetBondTag, (
    python::arg("molecule"),
    python::arg("bond_idx")));
  python::def<Tag (*)(const RDKit::Atom*)>("GetAtomTag", GetTag, (
    python::arg("atom")));
  python::def<Tag (*)(const RDKit::Bond*)>("GetBondTag", GetTag, (
    python::arg("bond")));

  python::def("HashMolecule", HashMolecule, (python::arg("molecule")));

  WrapMolecularKeys();
  WrapCircularAtomicEnvironment();
  WrapMolecularPerturbations();
  WrapMolecularPerturbationUtils();
  WrapChemicalDictionary();
  WrapMolecularConstraints();
  WrapMoleculePerturber();
  WrapAddons();
};
