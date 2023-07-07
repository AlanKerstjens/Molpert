import numpy as np
from typing import TypeVar, Union, Callable
from numpy.typing import NDArray

T = TypeVar("T") # State type
PopulationType = list[T]
ScoresType = NDArray[np.single]
SelectedType = NDArray[np.uint]
ObjectiveFunctionType = Callable[[T], float]
SelectionType = Callable[[PopulationType, ScoresType, int], SelectedType]
MutationType = Callable[[T], T]
CrossoverType = Callable[[T, T], T]
ReproductionType = Callable[[PopulationType, ScoresType, int],
                            tuple[PopulationType, ScoresType]]
SimilarityCriterionType = Union[Callable[[T, T], float], None]
TerminationCriterionType = Union[Callable[[PopulationType, ScoresType], bool], None]

class Reproduction:
    def __init__(self,
        mutation: MutationType,
        crossover: CrossoverType,
        parent_selection: SelectionType,
        crossover_probability: float = 0.1,
        similarity_criterion: SimilarityCriterionType = None,
        similarity_threshold: float = 1.0):
        self.mutation = mutation
        self.crossover = crossover
        self.parent_selection = parent_selection
        self.crossover_probability = crossover_probability
        self.similarity_criterion = similarity_criterion
        self.similarity_threshold = similarity_threshold

    def AcceptChild(self,
        population: PopulationType,
        child: T) -> bool:
        if not self.similarity_criterion or self.similarity_threshold >= 1.0:
            return True
        for individual in population:
            similarity = self.similarity_criterion(individual, child)
            if similarity > self.similarity_threshold:
                return False
        return True

    def __call__(self,
        population: PopulationType,
        scores: ScoresType,
        objective_function: ObjectiveFunctionType,
        n_children: int) -> tuple[PopulationType, ScoresType]:
        if n_children <= 0 or not len(population):
            return
        n_crossovers = 0
        if len(population) > 1:
            n_crossovers = n_children * self.crossover_probability
            if n_crossovers > 0.0 and n_crossovers < 1.0:
                if n_crossovers > np.random.rand():
                    n_crossovers = 1
            n_crossovers = int(n_crossovers)
        n_mutations = n_children - n_crossovers
        n_parents = n_mutations + n_crossovers * 2
        parent_indices = self.parent_selection(population, scores, n_parents)
        children_scores = []
        for i in range(n_mutations):
            parent = population[parent_indices[i]]
            child = self.mutation(parent)
            if self.AcceptChild(population, child):
                population.append(child)
                children_scores.append(objective_function(child))
        for i in range(n_mutations, n_parents, 2):
            parent1 = population[parent_indices[i]]
            parent2 = population[parent_indices[i+1]]
            child = self.crossover(parent1, parent2)
            if self.AcceptChild(population, child):
                population.append(child)
                children_scores.append(objective_function(child))
        return population, np.append(scores, children_scores)

class EvolutionaryAlgorithm:
    def __init__(self,
        reproduction: ReproductionType,
        survivor_selection: SelectionType,
        n_children: int,
        n_survivors: int):
        self.reproduction = reproduction
        self.survivor_selection = survivor_selection
        self.n_children = n_children
        self.n_survivors = n_survivors
    
    def __call__(self,
        population: PopulationType,
        objective_function: ObjectiveFunctionType,
        n_generations: int,
        termination_criterion: TerminationCriterionType = None):
        scores = np.fromiter(
            (objective_function(i) for i in population), np.single)
        for generation in range(n_generations):
            population, scores = self.reproduction(population, scores,
                objective_function, self.n_children)
            survivor_indices = self.survivor_selection(
                population, scores, self.n_survivors)
            population = [population[i] for i in survivor_indices]
            scores = scores[survivor_indices]
            if termination_criterion and termination_criterion(population, scores):
                break
        return population, scores, generation + 1

def FitnessProportionateSelection(
    population: PopulationType,
    scores: ScoresType, 
    n: int, 
    bias: float = 1.0,
    replace: bool = True) -> SelectedType:
    biased_scores = scores
    if bias != 1.0:
        biased_scores = np.power(scores, bias)
    scores_sum = np.sum(biased_scores)
    sampling_probabilities = None
    if scores_sum > 0.0:
        sampling_probabilities = biased_scores / scores_sum
    return np.random.choice(biased_scores.size, size=n,
        replace=replace, p=sampling_probabilities)

def TruncationSelection(
    population: PopulationType,
    scores: ScoresType,
    n: int) -> SelectedType:
    return np.argpartition(scores, -n)[-n:]

def BreakTiesTruncationSelection(
    population: PopulationType,
    scores: ScoresType,
    n: int) -> SelectedType:
    # Many individuals may have the same score. In those cases sorting 
    # algorithms retain the original order, and truncation always selects the 
    # same individuals. If you want to pick random individuals in case of score 
    # ties use this selection operator.
    sort_order = np.argsort(scores)
    selected = sort_order[-n:]
    selected_scores = scores[selected]
    worst_selected_score = selected_scores[0]
    n_worst_selection = np.searchsorted(
        selected_scores, worst_selected_score, side="right")
    if n_worst_selection == 1:
        return selected
    sorted_scores = scores[sort_order]
    first_worst_idx = np.searchsorted(
        sorted_scores, worst_selected_score, side="left")
    last_worst_idx = np.searchsorted(
        sorted_scores, worst_selected_score, side="right")
    worst_selection = np.random.choice(
       sort_order[first_worst_idx:last_worst_idx],
       size=n_worst_selection, replace=False)
    return np.concatenate((worst_selection, selected[n_worst_selection:]))

class StuckOrScoreThresholdTermination:
    def __init__(self, 
        max_n_generations_stuck: int,
        score_threshold: float):
        self.best_score = 0.0
        self.score_threshold = score_threshold
        self.n_generations = 0
        self.n_generations_stuck = 0
        self.max_n_generations_stuck = max_n_generations_stuck

    def __call__(self, 
        population: PopulationType,
        scores: ScoresType) -> bool:
        max_score = np.max(scores)
        self.n_generations += 1
        if max_score > self.best_score:
            self.best_score = max_score
            self.n_generations_stuck = 0
        else:
            self.n_generations_stuck += 1
        if self.best_score >= self.score_threshold:
            return True
        if self.n_generations_stuck > self.max_n_generations_stuck:
            return True
        if not self.n_generations % 100:
            print(f"Generation: {self.n_generations}, Score: {self.best_score}")
        return False