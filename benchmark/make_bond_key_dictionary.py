import pickle
import argparse
from rdkit import Chem
from molecular_constraints import LocalBondKey, RingAwareBondKey, \
                                  BondKeyDictionary

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Make a bond key dictionary/constraint generator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str,
        help="Path to input SMILES file.")
    parser.add_argument("output", type=str,
        help="Path to output bond key dictionary.")
    parser.add_argument("-t", "--bond_key_type", type=str, 
        choices=["local", "ring_aware"], default="local",
        help="Bond key type.")
    return parser.parse_args()

def Main():
    args = ParseArgs()

    if args.bond_key_type == "local":
        bond_key_generator = LocalBondKey
    elif args.bond_key_type == "ring_aware":
        bond_key_generator = RingAwareBondKey

    dictionary = BondKeyDictionary(bond_key_generator)

    supplier = Chem.SmilesMolSupplier(args.input, titleLine=False, nameColumn=-1)
    for molecule in supplier:
        if molecule:
            dictionary.AddMolecule(molecule)

    with open(args.output, "wb") as file:
        pickle.dump(dictionary, file)

if __name__ == "__main__":
    Main()