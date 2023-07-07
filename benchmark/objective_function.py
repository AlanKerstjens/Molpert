import molpert as mpt
from rdkit import Chem
from guacamol.scoring_function import \
  MoleculewiseScoringFunction, ScoringFunctionBasedOnRdkitMol

class ObjectiveFunction:
    def __init__(self, 
        objective: MoleculewiseScoringFunction):
        self.objective = objective
        self.n_calls = 0
        self.n_corrupt_scores = 0
        self.smiles_based = True
        # Certain GuacaMol scoring functions support scoring RDKit molecules
        # directly while others can only score SMILES.
        if isinstance(objective, ScoringFunctionBasedOnRdkitMol):
            self.smiles_based = False
    def __call__(self, molecule: Chem.Mol) -> float:
        sanitized_molecule = Chem.Mol(molecule)
        mpt.PartialSanitization(sanitized_molecule)
        if self.smiles_based:
            score = self.objective.score(Chem.MolToSmiles(sanitized_molecule))
        else:
            # score_mol calculates a "raw score" (e.g. the value of a descriptor).
            # modify_score puts "raw_score" in the [0,1] range based on how good
            # the raw score is (e.g. how close the descriptor is to a target value).
            score = self.objective.modify_score(
                self.objective.score_mol(sanitized_molecule))
        if score < 0.0:
            score = 0.0
            self.n_corrupt_scores += 1
        self.n_calls += 1
        return score