#ifndef _PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_
#define _PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_

#include "CircularAtomicEnvironment.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

EnvironmentKey CircularAtomicEnvironmentKeyWrapper(
  const CircularAtomicEnvironment& environment,
  const python::object& atom_hashes) {
  std::vector<std::uint64_t> ahs = to_vector<std::uint64_t>(atom_hashes);
  if (ahs.size() != environment.molecule->getNumAtoms()) {
    PyErr_SetString(PyExc_TypeError, 
      "Size mismatch between atom hashes and number of atoms");
    throw boost::python::error_already_set();
  };
  return environment.Key(ahs);
};

void WrapCircularAtomicEnvironment() {
  python::class_<CircularAtomicEnvironment>(
    "CircularAtomicEnvironment",
    python::init<const RDKit::Atom*, std::uint8_t>((
      python::arg("atom"), 
      python::arg("radius"))))
    .def<EnvironmentKey (CircularAtomicEnvironment::*)() const>(
      "Key", &CircularAtomicEnvironment::Key)
    .def("Key", CircularAtomicEnvironmentKeyWrapper)
    .def("SMILES", &CircularAtomicEnvironment::SMILES);
};

#endif // !_PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_