from random import choice, randint, uniform
import time
from pandas import array
import tsp_initialize as init
import tsp_functions as func

class mutation:

    def __init__(self, size_of_population, probability_of_mutation, type_of_mutation = '', set_of_individuals = []):
        self.size_of_population = size_of_population
        self.probability_of_mutation = probability_of_mutation
        self.type_of_mutation = type_of_mutation
        self.set_of_individuals = set_of_individuals 
        self.mutation_algorithm()

    def mutation_algorithm(self):
        if self.type_of_mutation == 'swap':
            self.swap_mutation()
        elif self.type_of_mutation == 'invert':
            self.invert_mutation()
        else:
            print("There is no \"", self.type_of_mutation, "\" type of mutation")

    def swap_mutation(self):
        r = uniform(0, 1)
        if r <= self.probability_of_mutation:
            for i in range(20):
                rand = randint(0, self.size_of_population-1)
                rand1 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                rand2 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                func.swap_positions(self.set_of_individuals[rand].genotype, rand1, rand2)
        return self.set_of_individuals

    def invert_mutation(self):
        r = uniform(0, 1)
        if r <= self.probability_of_mutation:
            for _ in range(5):
                rand = randint(0, self.size_of_population-1)
                rand1 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                rand2 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                func.invert(self.set_of_individuals[rand].genotype, rand1, rand2)
        return self.set_of_individuals
