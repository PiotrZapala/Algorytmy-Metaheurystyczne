from random import choice, randint, uniform
import time
from pandas import array
import tsp_initialize as init
import tsp_functions as func

class selection_population:

    def __init__(self, type_of_selection_population = '', population = [], best_solutions = []):
        self.type_of_selection_population = type_of_selection_population
        self.total_probability = 0
        self.total_fitness = 0
        self.field = 0
        self.roulette = []
        self.fitness = []
        self.phenotypes = []
        self.probability = []
        self.new_population = []
        self.population = population
        self.best_solutions = best_solutions
        self.selection_population_algorithm()

    def selection_population_algorithm(self):
        if self.type_of_selection_population == 'roulette':
            self.roulette_selection_population()
        elif self.type_of_selection_population == 'tournament':
            self.tournament_selection_population()
        else:
            print("There is no \"", self.type_of_selection_population, "\" type of selection population")

    def roulette_selection_population(self):
        for i in range(len(self.population)):
            self.total_fitness = self.total_fitness + self.population[i].fitness
        for n in range(len(self.population)):
            self.probability.append(self.population[n].fitness/self.total_fitness)
            self.field = self.field + self.probability[n]
            self.roulette.append(self.field)
        for _ in range(int(len(self.population)/2)):
            r = uniform(0, 1)
            for p in range(len(self.roulette)):
                if r <= self.roulette[p]:
                    self.new_population.append(self.population[p])
                    break
                else:
                    continue
        for i in range(len(self.new_population)):
            self.phenotypes.append(self.new_population[i].phenotype)
        self.best_solutions.append(min(self.phenotypes))
        return self.new_population

    def tournament_selection_population(self):
        pass
