import os
import json
import pickle
import argparse
import molpert as mpt
from rdkit import Chem
from evolutionary_algorithm import *
from objective_function import ObjectiveFunction
from molecule_manipulation import MoleculeSimilarity, MoleculePerturber, Mimicry
from molecular_constraints import ValenceConstraintGenerator
from guacamol.benchmark_suites import goal_directed_benchmark_suite

def ParseArgs():
    parser = argparse.ArgumentParser(description="Evaluates the performance of \
        a Molpert-based evolutionary algorithm with the GuacaMol benchmarks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("output", type=str,
        help="Path to output file or directory.")
    parser.add_argument("-vc", "--valence_constraints", action="store_true",
        help="Flag to use valence atom constraints.")
    parser.add_argument("-acg", "--atom_constraint_generator", type=str,
        help="Path to pickled atom constraint generator.")
    parser.add_argument("-bcg", "--bond_constraint_generator", type=str,
        help="Path to pickled bond constraint generator.")
    parser.add_argument("-s", "--suite", default="v2",
        choices=["trivial", "v1", "v2"],
        help="Name of the goal-directed GuacaMol benchmark suite.")
    parser.add_argument("-b", "--benchmarks", type=int, nargs="+", default=range(20),
        help="Indices of the benchmarks to run.")
    parser.add_argument("-r", "--replicas", type=int, default=1,
        help="Number of replicas to run for each benchmark.")
    args = parser.parse_args()
    if os.path.isfile(args.output) and \
        (len(args.benchmarks) > 1 or args.replicas > 1):
        raise parser.error(
            "output is a file, but multiple benchmarks or replicas are requested.")
    return args

def SizeBiasedTruncationSelection(
    population: PopulationType,
    scores: ScoresType,
    n: int) -> SelectedType:
    # Some scoring functions don't provide meaningful guidance to the algorithm
    # unless the molecule's size is substantial. This selection operation 
    # normalizes the scores by molecular size up until a certain size threshold.
    size_threshold = 6.0
    sizes = np.clip(np.fromiter(
        (m.GetNumAtoms() for m in population), np.single), 0.0, size_threshold)
    return np.argsort(scores * sizes / size_threshold)[-n:]

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

    population_size = 100
    perturber = MoleculePerturber(
        constraints=constraints, start_molecule_id=population_size)

    mimicry = Mimicry(
        constraints=constraints,
        min_fraction_atoms=0.1,
        max_fraction_atoms=0.5,
        start_molecule_id=2**31//2)

    similarity = MoleculeSimilarity()

    evolution = EvolutionaryAlgorithm(
        reproduction=Reproduction(
            mutation=perturber,
            crossover=mimicry,
            parent_selection=FitnessProportionateSelection,
            crossover_probability=0.25,
            similarity_criterion=similarity,
            similarity_threshold=0.90),
        survivor_selection=SizeBiasedTruncationSelection,
        n_children=population_size,
        n_survivors=population_size)

    benchmark_suite = goal_directed_benchmark_suite(args.suite)
    benchmark_ids = sorted([[bidx, 0] for bidx in args.benchmarks
                            if bidx < len(benchmark_suite)
                            for _ in range(args.replicas)])
    if len(benchmark_ids) > 1:
        for i in range(1, len(benchmark_ids)):
            if benchmark_ids[i-1][0] == benchmark_ids[i][0]:
                benchmark_ids[i][1] = benchmark_ids[i-1][1] + 1

    for benchmark_idx, replica_idx in benchmark_ids:
        benchmark = benchmark_suite[benchmark_idx]
        print(f"\nRunning benchmark {benchmark.name} (replica {replica_idx})")

        objective_function = ObjectiveFunction(benchmark.objective)

        perturber.n_molecules = 0
        perturber.design_time = 0
        mimicry.n_molecules = 0
        mimicry.design_time = 0
        similarity.clear()

        if benchmark.starting_population:
            population = [Chem.MolFromSmiles(s) for s in benchmark.starting_population]
        else:
            population = [Chem.Mol() for _ in range(population_size)]
        for i in range(len(population)):
            population[i].SetIntProp("ID", i)

        population, scores, n_generations = evolution(
            population=population,
            objective_function=objective_function,
            n_generations=10000,
            termination_criterion=StuckOrScoreThresholdTermination(
                max_n_generations_stuck=1000,
                score_threshold=1.0))

        for molecule in population:
            mpt.PartialSanitization(molecule)

        print(f"# generations: {n_generations}")
        print(f"Max score: {np.max(scores)}, Average score: {np.average(scores)}")

        output_file_path = args.output
        if os.path.isdir(args.output):
            output_file_path = os.path.join(
                args.output, f"{benchmark.name.replace(' ', '_')}")
            if args.replicas > 1:
                output_file_path += f"_{replica_idx}"
            output_file_path += ".json"
        with open(output_file_path, "w", encoding="utf-8") as file:
            data = {
                "population": [Chem.MolToSmiles(m) for m in population],
                "scores": [s for s in scores],
                "best_score": np.max(scores),
                "average_score": np.average(scores),
                "n_designed_molecules": perturber.n_molecules + mimicry.n_molecules,
                "design_time": perturber.design_time + mimicry.design_time,
                "n_objective_calls": objective_function.n_calls,
                "n_failed_objective_calls": objective_function.n_corrupt_scores}
            json.dump(data, file)

if __name__ == "__main__":
    Main()