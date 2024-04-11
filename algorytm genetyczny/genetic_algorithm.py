from random import choice, randint, uniform
import time
from pandas import array
import tsp_initialize as init
import tsp_functions as func
import selection 
import crossing
import selection_population
import mutation
import numpy as np

class genetic_algorithm:
    
    def __init__(self, type_of_selection, type_of_selection_population, type_of_crossing, type_of_mutation, probability_of_mutation, size_of_population, parents = [], best_solutions = [], distance_matrix = []):
        self.type_of_selection = type_of_selection
        self.type_of_selection_population = type_of_selection_population
        self.type_of_crossing = type_of_crossing
        self.type_of_mutation = type_of_mutation
        self.probability_of_mutation = probability_of_mutation
        self.size_of_population = size_of_population
        #self.type_of_stop_condition = type_of_stop_condition
        #self.number_for_stop = number_for_stop
        self.parents = parents
        self.best_solutions = best_solutions
        self.distance_matrix = distance_matrix
        self.children = []
        self.algorithm()

    def algorithm(self):
        #stop = stop_condition(self.type_of_stop_condition)
        #while(stop.stop_condition == False):
        for i in range(20000):
            sel = selection(self.type_of_selection, self.parents)
            fathers = sel.selected_fathers
            mothers = sel.selected_mothers
            cro = crossing(self.type_of_crossing, fathers, mothers, self.distance_matrix)
            self.children = cro.children
            self.population = np.concatenate((self.children, self.parents))
            sel_pop = selection_population(self.type_of_selection_population, self.population, self.best_solutions)
            self.parents = sel_pop.new_population
            mut = mutation(self.size_of_population, self.probability_of_mutation, self.type_of_mutation, self.parents)
            self.parents = mut.set_of_individuals
            #stop.update_stop_condition()
        return self.parents, self.best_solutions