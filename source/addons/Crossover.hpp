#ifndef _CROSSOVER_HPP_
#define _CROSSOVER_HPP_

#include "Subgraph.hpp"
#include "MolecularPerturbations.hpp"

std::pair<std::shared_ptr<SubgraphPerturbation>, 
          std::shared_ptr<SubgraphPerturbation>> Crossover(
  const RDKit::ROMol& molecule1,
  const RDKit::ROMol& molecule2,
  const SubgraphRequests& requests,
  std::mt19937& prng) {
  SubgraphRequests requests1 (requests);
  // We disable this behaviour because we want to do it manually.
  requests1.pad_attachment_points_if_insufficient = false;
  requests1.shuffle_attachment_points = false;

  boost::dynamic_bitset<> subgraph1 = RandomSubgraph(molecule1, requests1, prng);
  boost::dynamic_bitset<> inverted_subgraph1 = ~subgraph1;
  std::vector<AttachmentPoint> attachment_points1 = 
    AttachmentPoints(molecule1, subgraph1, requests1, prng);
  std::vector<AttachmentPoint> inverted_attachment_points1 = 
    AttachmentPoints(molecule1, inverted_subgraph1, requests1, prng);

  SubgraphRequests requests2 (requests1);
  requests2.n_attachment_points = inverted_attachment_points1.size();

  boost::dynamic_bitset<> subgraph2 = RandomSubgraph(molecule2, requests2, prng);
  boost::dynamic_bitset<> inverted_subgraph2 = ~subgraph2;
  std::vector<AttachmentPoint> attachment_points2 = 
    AttachmentPoints(molecule2, subgraph2, requests2, prng);
  std::vector<AttachmentPoint> inverted_attachment_points2 = 
    AttachmentPoints(molecule2, inverted_subgraph2, requests2, prng);

  requests1.n_attachment_points = inverted_attachment_points2.size();
  requests2.n_attachment_points = inverted_attachment_points1.size();
  PadAttachmentPoints(attachment_points1, molecule1, subgraph1, requests1, prng);
  PadAttachmentPoints(attachment_points2, molecule2, subgraph2, requests2, prng);
  if (requests.shuffle_attachment_points) {
    std::shuffle(attachment_points1.begin(), attachment_points1.end(), prng);
    std::shuffle(attachment_points2.begin(), attachment_points2.end(), prng);
  };

  assert(attachment_points1.size() >= inverted_attachment_points2.size());
  assert(attachment_points2.size() >= inverted_attachment_points1.size());

  std::shared_ptr<SubgraphPerturbation> perturbation1 (new SubgraphPerturbation(
    molecule1, molecule2, subgraph2, subgraph1));
  const SubgraphConstruction& subgraph_construction1 = 
    perturbation1->GetSubgraphConstruction();
  const std::unordered_map<AtomIdx, AtomIdx>& atom_mapping1 = 
    subgraph_construction1.GetAtomMapping(); 
  Tag max_bond_tag1 = GetMaxBondTag(molecule1) + 
    subgraph_construction1.GetBondConstructions().size();
  for (std::size_t i = 0; i < inverted_attachment_points1.size(); ++i) {
    perturbation1->AddBondConstruction(
      ++max_bond_tag1,
      inverted_attachment_points1[i].atom_idx,
      atom_mapping1.at(attachment_points2[i].atom_idx));
  };

  std::shared_ptr<SubgraphPerturbation> perturbation2 (new SubgraphPerturbation(
    molecule2, molecule1, subgraph1, subgraph2));
  const SubgraphConstruction& subgraph_construction2 = 
    perturbation2->GetSubgraphConstruction();
  const std::unordered_map<AtomIdx, AtomIdx>& atom_mapping2 = 
    subgraph_construction2.GetAtomMapping(); 
  Tag max_bond_tag2 = GetMaxBondTag(molecule2) + 
    subgraph_construction2.GetBondConstructions().size();
  for (std::size_t i = 0; i < inverted_attachment_points2.size(); ++i) {
    perturbation2->AddBondConstruction(
      ++max_bond_tag2,
      inverted_attachment_points2[i].atom_idx,
      atom_mapping2.at(attachment_points1[i].atom_idx));
  };

  return {std::move(perturbation1), std::move(perturbation2)};
};

#endif // _CROSSOVER_HPP_
