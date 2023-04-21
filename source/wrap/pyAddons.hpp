#pragma once
#ifndef _PY_ADDONS_HPP_
#define _PY_ADDONS_HPP_

#include "Mimicry.hpp"
#include "Crossover.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

void WrapAddons() {
  python::scope module_scope = python::scope();

  python::class_<SubgraphRequests>("SubgraphRequests")
    .def_readwrite("attachment_point_definition", &SubgraphRequests::attachment_point_definition)
    .def_readwrite("n_attachment_points", &SubgraphRequests::n_attachment_points)
    .def_readwrite("size", &SubgraphRequests::size)
    .def_readwrite("allow_breaking_cycles", &SubgraphRequests::allow_breaking_cycles)
    .def_readwrite("one_attachment_point_per_slot", &SubgraphRequests::one_attachment_point_per_slot)
    .def_readwrite("pad_attachment_points_if_insufficient", &SubgraphRequests::pad_attachment_points_if_insufficient)
    .def_readwrite("try_to_satisfy_valences_when_padding", &SubgraphRequests::try_to_satisfy_valences_when_padding)
    .def_readwrite("shuffle_attachment_points", &SubgraphRequests::shuffle_attachment_points);

  python::scope attachment_point_scope =
  python::class_<AttachmentPoint>("AttachmentPoint")
    .def_readwrite("atom_idx", &AttachmentPoint::atom_idx)
    .def_readwrite("type", &AttachmentPoint::type);

  python::enum_<AttachmentPoint::Definition>("Definition")
    .value("AllAtoms", AttachmentPoint::Definition::AllAtoms)
    .value("Valence", AttachmentPoint::Definition::Valence)
    .value("BrokenBonds", AttachmentPoint::Definition::BrokenBonds);
  
  python::scope addons_scope (module_scope);

  python::def("Mimicry", Mimicry, (
    python::arg("target_molecule"),
    python::arg("source_molecule"),
    python::arg("subgraph_requests"),
    python::arg("prng")));

  // I expose a specialized version of std::pair instead of registering a
  // std::pair to Python tuple converter. Ugly, but works.
  python::class_<std::pair<SubgraphPerturbation, SubgraphPerturbation>>(
    "SubgraphPerturbationPair", python::no_init)
    .def_readonly("first", &std::pair<SubgraphPerturbation, SubgraphPerturbation>::first)
    .def_readonly("second", &std::pair<SubgraphPerturbation, SubgraphPerturbation>::second);

  python::def("Crossover", Crossover, (
    python::arg("molecule1"),
    python::arg("molecule2"),
    python::arg("subgraph_requests"),
    python::arg("prng")));
};

#endif // !_PY_ADDONS_HPP_