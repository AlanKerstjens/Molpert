import pickle
import argparse
from rdkit import Chem
from functools import partial
from molecular_constraints import LocalAtomKey, RingAwareAtomKey, HashAtomKey, \
                                  EnvironmentKeyDictionary

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Make an environment key dictionary/constraint generator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str,
        help="Path to input SMILES file.")
    parser.add_argument("output", type=str,
        help="Path to output atom key dictionary.")
    parser.add_argument("-r", "--environment_radius", type=int, default=1,
        help="Circular atomic environment radius.")
    parser.add_argument("-t", "--atom_key_type", type=str, 
        choices=["local", "ring_aware"], default="local",
        help="Atom key type.")
    return parser.parse_args()

def Main():
    args = ParseArgs()

    if args.atom_key_type == "local":
        atom_key_generator = LocalAtomKey
    elif args.atom_key_type == "ring_aware":
        atom_key_generator = RingAwareAtomKey

    dictionary = EnvironmentKeyDictionary(
        atom_hasher=partial(HashAtomKey, atom_key_generator=atom_key_generator),
        env_radius=args.environment_radius)

    supplier = Chem.SmilesMolSupplier(args.input, titleLine=False, nameColumn=-1)
    for molecule in supplier:
        if molecule:
            dictionary.AddMolecule(molecule)

    with open(args.output, "wb") as file:
        pickle.dump(dictionary, file)

if __name__ == "__main__":
    Main()