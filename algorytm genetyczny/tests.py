from audioop import avg
from random import choice, randint, random, uniform
import time
import numpy as np
import tsp_initialize as init
import tsp_functions as func
from genetic import individual  
from genetic import selection 
from genetic import selection_population 
from genetic import crossing
from genetic import mutation 

def initi(file_name):
    number_of_cities, distance_matrix = init.read_file(file_name, type='tsp')
    return number_of_cities, distance_matrix

class population:

    def __init__(self, size_of_population, file_name):
        self.size_of_population = size_of_population
        self.size_of_individual = 0
        self.distance_matrix = []
        self.phenotype = []
        self.genotype = []
        self.set_of_individuals = []
        self.best_solution = []
        self.file_name = file_name
        self.initialize_population()

    def initialize_2opt_genotype(self):
        for _ in range(self.size_of_population):
            self.genotype.append(func.opt22(self.size_of_individual, self.distance_matrix))

    def initialize_genotype(self):
        for _ in range(self.size_of_population):
            self.genotype.append(func.random_solution(self.size_of_individual))

    def initialize_phenotype(self):
        for i in range(self.size_of_population):
            self.phenotype.append(func.get_weight(self.genotype[i], self.distance_matrix))

    def initialize_population(self):
        self.size_of_individual, self.distance_matrix = initi(self.file_name)
        self.initialize_2opt_genotype()
        self.initialize_phenotype()
        for i in range(len(self.genotype)):
            self.set_of_individuals.append(individual(self.genotype[i], self.phenotype[i], (1/self.phenotype[i])))
        index = self.phenotype.index(min(self.phenotype))
        self.best_solution.append(self.set_of_individuals[index])

        return self.set_of_individuals, self.best_solution, self.size_of_population, self.distance_matrix

class genetic_algorithm:
    
    def __init__(self, type_of_selection, type_of_selection_population, type_of_crossing, type_of_mutation, probability_of_mutation, size_of_population, iterations, best_known, parents = [], best_solutions = [], distance_matrix = []):
        self.type_of_selection = type_of_selection
        self.type_of_selection_population = type_of_selection_population
        self.type_of_crossing = type_of_crossing
        self.type_of_mutation = type_of_mutation
        self.probability_of_mutation = probability_of_mutation
        self.size_of_population = size_of_population
        self.iterations = iterations
        self.best_known = best_known
        self.parents = parents
        self.best_solutions = best_solutions
        self.distance_matrix = distance_matrix
        self.children = []
        self.current_best = []
        self.number_of_iterations_without_update = 0
        self.algorithm()


    def algorithm(self):
        i = 0
        self.current_best.append(self.best_solutions[0])
        while i < 4000: 
            sel = selection(self.type_of_selection, self.parents)
            fathers = sel.selected_fathers
            mothers = sel.selected_mothers
            cro = crossing(self.type_of_crossing, fathers, mothers, self.distance_matrix)
            self.children = cro.children
            if self.number_of_iterations_without_update > self.iterations:
                self.number_of_iterations_without_update = 0
                for _ in range(int(self.size_of_population/10)):
                    rand = randint(0, self.size_of_population-1)
                    child = func.opt23(self.size_of_population, self.distance_matrix, self.children[rand].genotype)
                    self.children.remove(self.children[rand])
                    self.children.insert(rand, individual(child, func.get_weight(child, self.distance_matrix), 1/func.get_weight(child, self.distance_matrix)))
            self.population = np.concatenate((self.children, self.parents))
            sel_pop = selection_population(self.type_of_selection_population, self.population, self.best_solutions)
            self.parents = sel_pop.new_population
            self.best_solutions = sel_pop.best_solutions
            if self.best_solutions[-1].phenotype < self.current_best[-1].phenotype:
                self.number_of_iterations_without_update += 1
                rand = randint(0, self.size_of_population-1)
                self.parents.remove(self.parents[rand])
                self.parents.insert(rand, self.current_best[-1])
                self.best_solutions.append(self.current_best[-1])
            else:
                self.current_best.append(self.best_solutions[-1])
            mut = mutation(self.distance_matrix, self.size_of_population, self.probability_of_mutation, self.type_of_mutation, self.parents)
            self.parents = mut.set_of_individuals
            i += 1
        return self.best_solutions, self.parents

def test(file_name, best_known_solution, size_of_population, type_of_crossing1, type_of_mutation1, type_of_selection1, type_of_selection_population1):
    probability_of_mutation = [15]
    number_of_iterations_without_update = [50]
    for it in range(len(size_of_population)):
        file = open('tests/'+file_name+'-size'+str(size_of_population[it])+'+'+type_of_crossing1+'+'+type_of_mutation1+'.txt', 'w')
        file.write('1.(size of population)' + ',' + '2.(probability of mutation)' + ',' + '3.(number of iterations without update)' + ',' + '4.(time)' + ',' + '5.(best solution)' + ',' + '6.(prd)' + '\n')
        size_of_population1 = size_of_population[it]
        pop1 = population(size_of_population1, file_name)
        parents1 = pop1.set_of_individuals
        best_solutions1 = pop1.best_solution
        size_of_population1 = pop1.size_of_population
        distance_matrix1 = pop1.distance_matrix
        for jt in range(len(probability_of_mutation)):
            probability_of_mutation1 = probability_of_mutation[jt]
            for ka in range(len(number_of_iterations_without_update)):
                number_of_iterations_without_update1 = number_of_iterations_without_update[ka]
                best_individual1 = []
                start = time.time()         
                gen = genetic_algorithm(str(type_of_selection1), str(type_of_selection_population1), str(type_of_crossing1), str(type_of_mutation1), probability_of_mutation1, size_of_population1, number_of_iterations_without_update1, best_known_solution, parents1, best_solutions1, distance_matrix1)
                end = time.time()
                my_time = end - start
                for j in range(len(gen.best_solutions)):
                    best_individual1.append(gen.best_solutions[j].phenotype)
                minn = min(best_individual1)
                prd = 100 * ((minn - best_known_solution) / best_known_solution)
                file.write(str(probability_of_mutation[jt]) + ',' + str(number_of_iterations_without_update[ka]) + ',' + str(my_time) + ',' + str(prd) + '\n')

def test2(file_name, best_known_solution, size_of_population, type_of_crossing, type_of_mutation, type_of_selection, type_of_selection_population):

        probability_of_mutation = 15
        number_of_iterations_without_update = 50
        file = open('test4/gr24'+'-'+str(type_of_crossing)+'-'+str(type_of_mutation)+'-'+str(type_of_selection)+'-'+str(type_of_selection_population)+'.txt', 'w')
        pop = population(size_of_population, file_name)
        parents = pop.set_of_individuals
        best_solutions = pop.best_solution
        size_of_population = pop.size_of_population
        distance_matrix1 = pop.distance_matrix
        result_time = []
        result_prd = []
        for _ in range(10):
            best_individual1 = []
            start = time.time()         
            gen = genetic_algorithm(str(type_of_selection), str(type_of_selection_population), str(type_of_crossing), str(type_of_mutation), probability_of_mutation, size_of_population, number_of_iterations_without_update, best_known_solution, parents, best_solutions, distance_matrix1)
            end = time.time()
            my_time = end - start
            for j in range(len(gen.best_solutions)):
                    best_individual1.append(gen.best_solutions[j].phenotype)
            minn = min(best_individual1)
            prd = 100 * ((minn - best_known_solution) / best_known_solution)
            result_time.append(my_time)
            result_prd.append(prd)
        avgtime = sum(result_time)/len(result_time)
        avgprd = sum(result_prd)/len(result_prd)
        file.write(str(avgtime)+','+str(avgprd) +'\n')

if __name__ == '__main__':
    file_name = 'gr24.tsp'
    size = 96
    best_known_solution = 1272
    ox='OX'
    pmx = 'PMX'
    swap='swap'
    invert = 'invert'
    random ='random'
    roulette = 'roulette'
    tournament = 'tournament'

    test2(file_name, best_known_solution, size, ox, swap, random, random)
    test2(file_name, best_known_solution, size, ox, invert, random, random)
    test2(file_name, best_known_solution, size, ox, swap, roulette, roulette)
    test2(file_name, best_known_solution, size, ox, invert, roulette, roulette)
    test2(file_name, best_known_solution, size, ox, swap, tournament, tournament)
    test2(file_name, best_known_solution, size, ox, invert, tournament, tournament)
    test2(file_name, best_known_solution, size, pmx, swap, random, random)
    test2(file_name, best_known_solution, size, pmx, invert, random, random)
    test2(file_name, best_known_solution, size, pmx, swap, roulette, roulette)
    test2(file_name, best_known_solution, size, pmx, invert, roulette, roulette)
    test2(file_name, best_known_solution, size, pmx, swap, tournament, tournament)
    test2(file_name, best_known_solution, size, pmx, invert, tournament, tournament)
  
