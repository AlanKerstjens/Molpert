#pragma once
#ifndef _PY_VALENCE_HPP_
#define _PY_VALENCE_HPP_

#include "Valence.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapValence() {

  python::def("ExplicitValence", ExplicitValence, (
    python::arg("atom"), 
    python::arg("include_hydrogens") = true,
    python::arg("aromatic_is_single") = false));
  python::def("IsValenceBelowMax", IsValenceBelowMax, (
    python::arg("atomic_number"), 
    python::arg("valence")));

};

#endif // !_PY_VALENCE_HPP_