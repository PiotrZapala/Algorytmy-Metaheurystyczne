from random import choice, randint, uniform
import time
import numpy as np
import tsp_initialize as init
import tsp_functions as func
import concurrent.futures

class individual:

    def __init__(self, genotype, phenotype, fitness):
        self.genotype = genotype
        self.phenotype = phenotype
        self.fitness = fitness

class population:

    def __init__(self, size_of_population):
        self.size_of_population = size_of_population
        self.size_of_individual = 0
        self.distance_matrix = []
        self.phenotype = []
        self.genotype = []
        self.set_of_individuals = []
        self.best_solution = []
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
        self.size_of_individual, self.distance_matrix = func.initialize_problem()
        self.initialize_genotype()
        self.initialize_phenotype()
        for i in range(len(self.genotype)):
            self.set_of_individuals.append(individual(self.genotype[i], self.phenotype[i], (1/self.phenotype[i])))
        index = self.phenotype.index(min(self.phenotype))
        self.best_solution.append(self.set_of_individuals[index])

        return self.set_of_individuals, self.best_solution, self.size_of_population, self.distance_matrix

class selection:

    def __init__(self, type_of_selection = '', set_of_individuals = []):
        self.type_of_selection = type_of_selection
        self.set_of_individuals = set_of_individuals
        self.total_fitness = 0
        self.field = 0
        self.fitness = []
        self.roulette = []
        self.selected_fathers = []
        self.selected_mothers = []
        self.probability = []
        self.selection_algorithm()

    def selection_algorithm(self):
        if self.type_of_selection == 'random':
            self.random_selection()
        elif self.type_of_selection == 'roulette':
            self.roulette_selection()
        elif self.type_of_selection == 'tournament':
            self.tournament_selection()
        else:
            print("There is no \"", self.type_of_selection, "\" type of selection")

    def random_selection(self):
        size_of_population = len(self.set_of_individuals)
        values = list(range(size_of_population))
        for _ in range(int(size_of_population/2)):
            rand1 = choice(values)
            self.selected_fathers.append(self.set_of_individuals[rand1])
            values.remove(rand1)
            rand2 = choice(values)
            self.selected_mothers.append(self.set_of_individuals[rand2])
            values.remove(rand2)
        return self.selected_fathers, self.selected_mothers

    def roulette_selection(self):
        for i in range(len(self.set_of_individuals)):
            self.total_fitness = self.total_fitness + self.set_of_individuals[i].fitness
        for n in range(len(self.set_of_individuals)):
            self.probability.append(self.set_of_individuals[n].fitness/self.total_fitness)
            self.field = self.field + self.probability[n]
            self.roulette.append(self.field)
        for _ in range(int(len(self.set_of_individuals)/2)):
            r1 = uniform(0,1)
            for m in range(len(self.roulette)):
                if r1 <= self.roulette[m]:
                    self.selected_fathers.append(self.set_of_individuals[m])
                    break
                else:
                    continue  
            r2 = uniform(0,1)
            for n in range(len(self.roulette)):
                if r2 <= self.roulette[n]:
                    self.selected_mothers.append(self.set_of_individuals[n])
                    break
                else:
                    continue  

        return self.selected_fathers, self.selected_mothers

    def tournament_selection(self):
        competition1 = []
        competition2 = []
        individuals1 = []
        individuals2 = []
        size_of_population = len(self.set_of_individuals)
        values = list(range(size_of_population))
        for _ in range(int(size_of_population/2)):
            for _ in range(int(size_of_population/20)):
                rand_father = choice(values)
                competition1.append(self.set_of_individuals[rand_father].fitness)
                individuals1.append(self.set_of_individuals[rand_father])
                rand_mother = choice(values)
                competition2.append(self.set_of_individuals[rand_mother].fitness)
                individuals2.append(self.set_of_individuals[rand_mother])
            best_father = max(competition1)
            index1 = competition1.index(best_father)
            self.selected_fathers.append(individuals1[index1])
            best_mother = max(competition2)
            index2 = competition2.index(best_mother)
            self.selected_mothers.append(individuals1[index2])

        return self.selected_fathers, self.selected_mothers

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
        if self.type_of_selection_population == 'random':
            self.random_selection_population()
        elif self.type_of_selection_population == 'roulette':
            self.roulette_selection_population()
        elif self.type_of_selection_population == 'tournament':
            self.tournament_selection_population()
        else:
            print("There is no \"", self.type_of_selection_population, "\" type of selection population")

    def random_selection_population(self):
        size_of_population = len(self.population)
        values = list(range(size_of_population))
        for _ in range(int(size_of_population/2)):
            rand = choice(values)
            self.new_population.append(self.population[rand])
            values.remove(rand)

        return self.new_population

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
        index = self.phenotypes.index(min(self.phenotypes))
        self.best_solutions.append(self.new_population[index])

        return self.new_population

    def tournament_selection_population(self):
        competition = []
        individuals = []
        size_of_population = len(self.population)
        values = list(range(size_of_population))
        for _ in range(int(size_of_population/2)):
            for _ in range(int(size_of_population/20)):
                rand = choice(values)
                competition.append(self.population[rand].fitness)
                individuals.append(self.population[rand])
            best_individual = max(competition)
            index = competition.index(best_individual)
            self.new_population.append(individuals[index])

        return self.new_population

class crossing:
    def __init__(self, type_of_crossing = '', selected_fathers = [], selected_mothers = [], distance_matrix = []):
        self.type_of_crossing = type_of_crossing
        self.selected_fathers = selected_fathers
        self.selected_mothers = selected_mothers
        self.distance_matrix = distance_matrix
        self.children = []
        self.crossing_algorithm()

    def crossing_algorithm(self):
        if self.type_of_crossing == 'HX':
            self.half_crossover()
        elif self.type_of_crossing == 'OX':
            self.order_crossover()
        elif self.type_of_crossing == 'PMX':
            self.partially_mapped_crossover()
        else:
            print("There is no \"", self.type_of_crossing, "\" type of crossing")

    def half_crossover(self):
        pass

    def order_crossover(self):
        def crossover(parent1, parent2, start, stop):
            child = [None]*len(parent1)
            child[start:stop] = child[start:stop]
            index1 = stop
            index2 = stop
            length = len(parent1)
            while None in child:
                if parent2[index1%length] not in child:
                    child[index2%length] = parent2[index1%length]
                    index2 += 1
                index1 += 1
            return child
        
        for i in range(len(self.selected_fathers)):       
            half = len(self.selected_fathers[i].genotype) // 2
            start = randint(0, len(self.selected_fathers[i].genotype)-half)
            stop = start + half
            child1_genotype = crossover(self.selected_fathers[i].genotype, self.selected_mothers[i].genotype, start, stop)
            child2_genotype = crossover(self.selected_mothers[i].genotype, self.selected_fathers[i].genotype, start, stop)
            child1_phenotype = func.get_weight(child1_genotype, self.distance_matrix)
            child2_phenotype = func.get_weight(child2_genotype, self.distance_matrix)
            self.children.append(individual(child1_genotype, child1_phenotype, 1/child1_phenotype))
            self.children.append(individual(child2_genotype, child2_phenotype, 1/child2_phenotype))

        return self.children

    def partially_mapped_crossover(self):
        def __pmx(parent1, parent2, start, stop):
            child = [None]*len(parent1)
            child[start:stop] = parent1[start:stop]
            for index1, x in enumerate(parent2[start:stop]):
                index1 += start
                if x not in child:
                    while child[index1] != None:
                        index1 = parent2.index(parent1[index1])
                    child[index1] = x
            for index1, x in enumerate(child):
                if x == None:
                    child[index1] = parent2[index1]
            return child

        for i in range(len(self.selected_fathers)):
            half = len(self.selected_fathers[i].genotype) // 2
            start = randint(0, len(self.selected_fathers[i].genotype)-half)
            stop = start + half
            child1_genotype = __pmx(self.selected_fathers[i].genotype, self.selected_mothers[i].genotype, start, stop)
            child2_genotype = __pmx(self.selected_mothers[i].genotype, self.selected_fathers[i].genotype, start, stop)
            child1_phenotype = func.get_weight(child1_genotype, self.distance_matrix)
            child2_phenotype = func.get_weight(child2_genotype, self.distance_matrix)
            self.children.append(individual(child1_genotype, child1_phenotype, 1/child1_phenotype))
            self.children.append(individual(child2_genotype, child2_phenotype, 1/child2_phenotype))

        return self.children

class mutation:

    def __init__(self, distance_matrix, size_of_population, probability_of_mutation, type_of_mutation = '', set_of_individuals = []):
        self.distance_matrix = distance_matrix
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
        r = randint(1, 100)
        if r <= self.probability_of_mutation:
            for _ in range(int(self.size_of_population/10)):
                rand = randint(0, self.size_of_population-1)
                rand1 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                rand2 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                tour = func.swap_positions(self.set_of_individuals[rand].genotype, rand1, rand2)
                self.set_of_individuals.remove(self.set_of_individuals[rand])
                self.set_of_individuals.insert(rand,individual(tour, func.get_weight(tour, self.distance_matrix), 1/func.get_weight(tour, self.distance_matrix)))

        return self.set_of_individuals

    def invert_mutation(self):
        r = randint(1, 100)
        if r <= self.probability_of_mutation:
            for _ in range(int(self.size_of_population/10)):
                rand = randint(0, self.size_of_population-1)
                rand1 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                rand2 = randint(0, len(self.set_of_individuals[rand].genotype)-1)
                tour = func.invert(self.set_of_individuals[rand].genotype, rand1, rand2)
                self.set_of_individuals.remove(self.set_of_individuals[rand])
                self.set_of_individuals.insert(rand, individual(tour, func.get_weight(tour, self.distance_matrix), 1/func.get_weight(tour, self.distance_matrix)))

        return self.set_of_individuals

class stop_condition:

    def __init__(self, type_of_stop_condition, number_for_stop):
        self.type_of_stop_condition = type_of_stop_condition
        self.number_for_stop = number_for_stop
        self.stop_condition = False
        self.start_time = time.time()
        self.iterations = 0
        self.iterations_without_update = 0

    def selection_stop_condition(self):
        if self.type_of_stop_condition == '0':
            self.stop_condition_0()
        elif self.type_of_stop_condition == '1':
            self.stop_condition_1()
        elif self.type_of_stop_condition == '2':
            self.stop_condition_2()
        else:
            print("There is no \"", self.type_of_stop_condition, "\" type of stop condition")

    def update_stop_condition(self):
        self.selection_stop_condition()

    def stop_condition_0(self):
        if self.iterations_without_update > self.number_for_stop:
            self.stop_condition = True
            self.iterations_without_update += 1
            print('stop_condition_0 works!')
        
    def stop_condition_1(self):
        self.iterations += 1
        if self.iterations >= self.number_for_stop:
            self.stop_condition = True
        print('stop_condition_1 works!')

    def stop_condition_2(self):
        end = time.time()
        if end - self.start_time > self.number_for_stop:
            self.stop_condition = True
        print('stop_condition_2 works!')

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
        while i < 2000: #self.global_iterations_without_update
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

def basic():
    condition = True
    print()
    print("Genetic Algorithm for 'Travelling salesman problem'\n")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
    while condition:
        print()
        print("To continue press '1'\nTo exit press '0'")
        option = int(input().strip())
        if option == 1:
            print("Enter size of population:")
            size_of_population = int(input().strip())
            print("Enter the type of selection individuals for the crossing:")
            print("1.random\n2.roulette\n3.tournament")
            type_of_selection = str(input().strip())
            print("Enter the type of selection the new generation:")
            print("1.random\n2.roulette\n3.tournament")
            type_of_selection_population = str(input().strip())
            print("Enter the type of the crossing:")
            print("1.HX\n2.OX\n3.PMX")
            type_of_crossing = str(input().strip())
            print("Enter the type of the mutation:")
            print("1.swap\n2.invert")
            type_of_mutation = str(input().strip())
            print("Enter the probability of the mutation:")
            probability_of_mutation = int(input().strip())
            print("Enter number of iterations without update:")
            number_of_iterations_without_update = int(input().strip())
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
            best = []
            pop = population(size_of_population)
            parents = pop.set_of_individuals
            best_solutions = pop.best_solution
            size_of_population = pop.size_of_population
            distance_matrix = pop.distance_matrix
            gen = genetic_algorithm(str(type_of_selection), str(type_of_selection_population), str(type_of_crossing), str(type_of_mutation), probability_of_mutation, size_of_population, number_of_iterations_without_update, parents, best_solutions, distance_matrix)
            for i in range(len(gen.best_solutions)):
                best.append(gen.best_solutions[i].phenotype)
            print(min(best))
        else:
            condition = False

def islands():
    condition = True
    print()
    print("Genetic Algorithm for 'Travelling salesman problem'\n")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
    while condition:
        print()
        print("To continue press '1'\nTo exit press '0")
        option = int(input().strip())
        if option == 1:
            print("Enter size of populations:")
            size_of_population1 = int(input().strip())
            size_of_population2 = int(input().strip())
            pop1 = population(size_of_population1)
            pop2 = population(size_of_population1)
            print("Enter the types of selection individuals for the crossing:")
            print("1.random\n2.roulette\n3.tournament")
            type_of_selection1 = str(input().strip())
            type_of_selection2 = str(input().strip())
            print("Enter the types of selection the new generation:")
            print("1.random\n2.roulette\n3.tournament")
            type_of_selection_population1 = str(input().strip())
            type_of_selection_population2 = str(input().strip())
            print("Enter the types of the crossing:")
            print("1.HX\n2.OX\n3.PMX")
            type_of_crossing1 = str(input().strip())
            type_of_crossing2 = str(input().strip())
            print("Enter the types of the mutation:")
            print("1.swap\n2.invert")
            type_of_mutation1 = str(input().strip())
            type_of_mutation2 = str(input().strip())
            print("Enter the probabilities of the mutation:")
            probability_of_mutation1 = int(input().strip())
            probability_of_mutation2 = int(input().strip())
            print("Enter number of iterations without update:")
            number_of_iterations_without_update1 = int(input().strip())
            number_of_iterations_without_update2 = int(input().strip())   
            print("Enter best known solution:")
            best_known = int(input().strip())
            print("Enter number of migrations:")
            number_of_migrations = int(input().strip())        
            parents1 = pop1.set_of_individuals
            parents2 = pop2.set_of_individuals
            best_solutions1 = pop1.best_solution
            best_solutions2 = pop2.best_solution
            size_of_population1 = pop1.size_of_population
            size_of_population2 = pop2.size_of_population
            distance_matrix1 = pop1.distance_matrix
            distance_matrix2 = pop2.distance_matrix
            best_individual1 = []
            best_individual2 = []
            with concurrent.futures.ThreadPoolExecutor() as executor:   
                for _ in range(number_of_migrations):         
                    gen1 = executor.submit(genetic_algorithm, str(type_of_selection1), str(type_of_selection_population1), str(type_of_crossing1), str(type_of_mutation1), probability_of_mutation1, size_of_population1, number_of_iterations_without_update1, best_known, parents1, best_solutions1, distance_matrix1)
                    gen2 = executor.submit(genetic_algorithm, str(type_of_selection2), str(type_of_selection_population2), str(type_of_crossing2), str(type_of_mutation2), probability_of_mutation2, size_of_population2, number_of_iterations_without_update2, best_known, parents2, best_solutions2, distance_matrix2)
                    results1 = gen1.result()
                    results2 = gen2.result()
                    parents1 = results1.parents
                    parents2 = results2.parents
                    parents1, parents2 = migration(parents1, parents2, size_of_population1)
            length1 = len(results1.best_solutions)
            length2 = len(results2.best_solutions)
            for j in range(length1):
                best_individual1.append(results1.best_solutions[j].phenotype)
            for k in range(length2):
                best_individual2.append(results2.best_solutions[k].phenotype)
            print(min(best_individual1))
            print(min(best_individual2))
        else:
            condition = False

def main():
    print()
    print("Genetic Algorithm for 'Travelling salesman problem'\n")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
    print()
    print("Enter the type of algorithm\n1.basic\n2.island")
    type_of_algorithm = int(input().strip())
    if type_of_algorithm == 1:
        basic()
    elif type_of_algorithm == 2:
        islands()
    else:
        main()

if __name__ == '__main__':
    main()