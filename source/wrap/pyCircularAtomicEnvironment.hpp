#pragma once
#ifndef _PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_
#define _PY_CIRCULAR_ATOMIC_ENVIRONMENT_HPP_

#include "CircularAtomicEnvironment.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapCircularAtomicEnvironment() {
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