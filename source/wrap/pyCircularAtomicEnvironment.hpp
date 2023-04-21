#pragma once
#ifndef _PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_
#define _PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_

#include "CircularAtomicEnvironment.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapCircularAtomicEnvironment() {
  python::class_<std::pair<EnvironmentKey, EnvironmentKey>>(
    "EnvironmentKeyChange", python::init<EnvironmentKey, EnvironmentKey>())
    .def_readwrite("prior", &std::pair<EnvironmentKey, EnvironmentKey>::first)
    .def_readwrite("posterior", &std::pair<EnvironmentKey, EnvironmentKey>::second);

  python::class_<CircularAtomicEnvironment>(
    "CircularAtomicEnvironment",
    python::init<const RDKit::Atom*, std::uint8_t>((
      python::arg("atom"), 
      python::arg("radius"))))
    .def<EnvironmentKey (CircularAtomicEnvironment::*)() const>(
      "Key", &CircularAtomicEnvironment::Key)
    .def("SMILES", &CircularAtomicEnvironment::SMILES);
};

#endif // !_PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_