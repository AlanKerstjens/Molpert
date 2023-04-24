#include "ChemicalDictionary.hpp"
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>

int main(int argc, char* argv[]) {
	if (argc < 3) {
		return 1;
	};

	bool kekulize = false;
  unsigned environment_radius = 2;
	std::string input_file_path (argv[1]);
	std::string output_file_path (argv[2]);
  if (argc > 3) {
    environment_radius = std::stoi(argv[3]);
  };

	RDKit::SmilesMolSupplier supplier(input_file_path, " \t", 0, 1, false, true);

  ChemicalDictionary dictionary (environment_radius);

  std::cout << "Creating ChemicalDictionary with circular atomic environments "
            << "of radius " << environment_radius << std::endl;

	RDKit::ROMOL_SPTR molecule;
	while (!supplier.atEnd()) {
		molecule.reset(supplier.next());
    if (!molecule) {
      continue;
    };
    if (kekulize) {
      RDKit::RWMol kekulized_molecule (*molecule);
      if (!RDKit::MolOps::KekulizeIfPossible(kekulized_molecule)) {
        continue;
      };
      dictionary.AddMolecule(kekulized_molecule);
    } else {
      dictionary.AddMolecule(*molecule);
    };
	};
  dictionary.BuildPartialKeyDictionaries();

  dictionary.Save(output_file_path);

	return 0;
};