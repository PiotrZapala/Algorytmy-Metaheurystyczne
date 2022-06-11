from random import choice, randint, uniform
import time
import numpy as np
import tsp_initialize as init
import tsp_functions as func
import concurrent.futures
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

def migration(parents1, parents2, size):
    for _ in range(int(size/5)):
        rand1 = randint(0, size-1)
        rand2 = randint(0, size-1)
        ind1 = parents1[rand1]
        ind2 = parents2[rand2]
        parents1.append(ind2)
        parents1.remove(ind1)
        parents2.append(ind1)
        parents2.remove(ind2)
    return parents1, parents2

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
        #self.type_of_stop_condition = type_of_stop_condition
        #self.number_for_stop = number_for_stop
        self.parents = parents
        self.best_solutions = best_solutions
        self.distance_matrix = distance_matrix
        self.children = []
        self.current_best = []
        self.number_of_iterations_without_update = 0
        self.algorithm()


    def algorithm(self):
        #stop = stop_condition(self.type_of_stop_condition)
        i = 0
        self.current_best.append(self.best_solutions[0])
        while i < 4000: #self.global_iterations_without_update
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
            #stop.update_stop_condition()
            i += 1
        return self.best_solutions, self.parents

def test():
    size_of_population = [40, 80, 120, 160, 200, 240]
    probability_of_mutation = [5, 10, 15, 20]
    number_of_iterations_without_update = [50, 100, 150, 200]
    best_known_solution = 7544
    file_name = 'berlin52.tsp'
    type_of_selection1 = 'roulette'
    type_of_selection2 = 'roulette' 
    type_of_selection_population1 = 'roulette'
    type_of_selection_population2 = 'roulette'
    type_of_crossing1 = 'PMX'
    type_of_crossing2 = 'PMX'
    type_of_mutation1 = 'swap'
    type_of_mutation2 = 'swap'
    for it in range(len(size_of_population)):
        file = open('tests/'+file_name+'-size'+str(size_of_population[it])+'.txt', 'w')
        file.write('1.(size of population)' + ',' + '2.(probability of mutation)' + ',' + '3.(number of iterations without update)' + ',' + '4.(time)' + ',' + '5.(best solution)' + ',' + '6.(prd)' + '\n')
        size_of_population1 = size_of_population[it]
        size_of_population2 = size_of_population[it]
        pop1 = population(size_of_population1, file_name)
        pop2 = population(size_of_population1, file_name)
        best_known = best_known_solution
        parents1 = pop1.set_of_individuals
        parents2 = pop2.set_of_individuals
        best_solutions1 = pop1.best_solution
        best_solutions2 = pop2.best_solution
        size_of_population1 = pop1.size_of_population
        size_of_population2 = pop2.size_of_population
        distance_matrix1 = pop1.distance_matrix
        distance_matrix2 = pop2.distance_matrix
        for jt in range(len(probability_of_mutation)):
            probability_of_mutation1 = probability_of_mutation[jt]
            probability_of_mutation2 = probability_of_mutation[jt]
            for ka in range(len(number_of_iterations_without_update)):
                number_of_iterations_without_update1 = number_of_iterations_without_update[ka]
                number_of_iterations_without_update2 = number_of_iterations_without_update[ka]
                best_individual1 = []
                best_individual2 = []
                minn = []
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    start = time.time()         
                    gen1 = executor.submit(genetic_algorithm, str(type_of_selection1), str(type_of_selection_population1), str(type_of_crossing1), str(type_of_mutation1), probability_of_mutation1, size_of_population1, number_of_iterations_without_update1, best_known, parents1, best_solutions1, distance_matrix1)
                    gen2 = executor.submit(genetic_algorithm, str(type_of_selection2), str(type_of_selection_population2), str(type_of_crossing2), str(type_of_mutation2), probability_of_mutation2, size_of_population2, number_of_iterations_without_update2, best_known, parents2, best_solutions2, distance_matrix2)
                    my_time = end - start
                    end = time.time()
                    results1 = gen1.result()
                    results2 = gen2.result()
                length1 = len(results1.best_solutions)
                length2 = len(results2.best_solutions)
                for j in range(length1):
                    best_individual1.append(results1.best_solutions[j].phenotype)
                for k in range(length2):
                    best_individual2.append(results2.best_solutions[k].phenotype)
                minn.append(min(best_individual1))
                minn.append(min(best_individual2))
                beesstt = min(minn)
                prd = 100 * ((beesstt - best_known) / best_known)
                file.write(str(size_of_population[it]) + ',' + str(probability_of_mutation[jt]) + ',' + str(number_of_iterations_without_update[ka]) + ',' + str(my_time) + ',' + str(beesstt) +  ',' + str(prd) + '\n')

if __name__ == '__main__':
    test()
