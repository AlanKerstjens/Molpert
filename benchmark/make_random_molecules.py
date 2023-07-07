import pickle
import argparse
import molpert as mpt
from rdkit import Chem
from molecule_manipulation import MoleculePerturber
from molecular_constraints import ValenceConstraintGenerator

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Make random molecules using specific constraints",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("output", type=str,
        help="Path to output SMILES file containing designed molecules.")
    parser.add_argument("-vc", "--valence_constraints", action="store_true",
        help="Flag to use valence atom constraints.")
    parser.add_argument("-acg", "--atom_constraint_generator", type=str,
        help="Path to pickled atom constraint generator.")
    parser.add_argument("-bcg", "--bond_constraint_generator", type=str,
        help="Path to pickled bond constraint generator.")
    parser.add_argument("-n", "--n_molecules", type=int, default=1,
        help="Number of molecules to generate.")
    parser.add_argument("-na", "--n_atoms", type=int, default=29,
        help="Number of atoms per molecule.")
    parser.add_argument("-nb", "--n_bonds", type=int, default=32,
        help="Number of bonds per molecule.")
    return parser.parse_args()

def Main():
    args = ParseArgs()

    atom_constraint_generator = None
    bond_constraint_generator = None
    if args.valence_constraints:
        atom_constraint_generator = ValenceConstraintGenerator
    elif args.atom_constraint_generator:
        with open(args.atom_constraint_generator, "rb") as file:
            atom_constraint_generator = pickle.load(file)
    if args.bond_constraint_generator:
        with open(args.bond_constraint_generator, "rb") as file:
            bond_constraint_generator = pickle.load(file)

    constraints = None
    if atom_constraint_generator or bond_constraint_generator:
        constraints = mpt.MolecularConstraints(
            atom_constraint_generator=atom_constraint_generator,
            bond_constraint_generator=bond_constraint_generator)

    perturber = MoleculePerturber(constraints)

    with open(args.output, "w") as file:
        for i in range(args.n_molecules):
            print(i)
            molecule = perturber.RandomMolecule(
                Chem.MolFromSmiles("C"), args.n_atoms, args.n_bonds)
            mpt.PartialSanitization(molecule)
            file.write(f"{Chem.MolToSmiles(molecule)}\n")

if __name__ == "__main__":
    Main()
