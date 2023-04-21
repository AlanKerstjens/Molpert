#pragma once
#ifndef _MIMICRY_HPP_
#define _MIMICRY_HPP_

#include "Subgraph.hpp"
#include "MolecularPerturbations.hpp"

std::shared_ptr<SubgraphPerturbation> Mimicry(
  const RDKit::ROMol& target,
  const RDKit::ROMol& source,
  const SubgraphRequests& requests,
  std::mt19937& prng) {
  if (!source.getNumAtoms()) {
    return nullptr;
  };
  boost::dynamic_bitset<> target_subgraph = RandomSubgraph(
    target, requests, prng);
  boost::dynamic_bitset<> inverted_target_subgraph = ~target_subgraph;
  std::vector<AttachmentPoint> inverted_target_attachment_points = 
    AttachmentPoints(target, inverted_target_subgraph, requests, prng);

  SubgraphRequests source_requests (requests);
  source_requests.n_attachment_points = inverted_target_attachment_points.size();
  source_requests.pad_attachment_points_if_insufficient = true;
  source_requests.shuffle_attachment_points = true;

  boost::dynamic_bitset<> source_subgraph = RandomSubgraph(
    source, source_requests, prng);
  std::vector<AttachmentPoint> source_attachment_points = AttachmentPoints(
    source, source_subgraph, source_requests, prng);

  assert(source_attachment_points.size() >= inverted_target_attachment_points.size());

  std::shared_ptr<SubgraphPerturbation> perturbation (new SubgraphPerturbation(
    target, source, source_subgraph, target_subgraph));
  const SubgraphConstruction& subgraph_construction = 
    perturbation->GetSubgraphConstruction();
  const std::unordered_map<AtomIdx, AtomIdx>& atom_mapping = 
    subgraph_construction.GetAtomMapping();
  Tag max_bond_tag = GetMaxBondTag(target) + 
    subgraph_construction.GetBondConstructions().size();
  for (std::size_t i = 0; i < inverted_target_attachment_points.size(); ++i) {
    perturbation->AddBondConstruction(
      ++max_bond_tag,
      inverted_target_attachment_points[i].atom_idx,
      atom_mapping.at(source_attachment_points[i].atom_idx));
  };

  return perturbation;
};

#endif // _MIMICRY_HPP_