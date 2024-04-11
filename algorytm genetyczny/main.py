from random import choice, randint, uniform
import time
from pandas import array
import tsp_initialize as init
import tsp_functions as func
import population
import genetic_algorithm

def main():
    condition = True
    print()
    print("Genetic Algorithm for 'Travelling salesman problem'\n")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
    while condition:
        print()
        print("To continue press '1'\nTo exit press '0")
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
            probability_of_mutation = float(input().strip())
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
            best = []
            pop = population(size_of_population)
            parents = pop.set_of_individuals
            best_solutions = pop.best_solution
            size_of_population = pop.size_of_population
            distance_matrix = pop.distance_matrix
            gen = genetic_algorithm(str(type_of_selection), str(type_of_selection_population), str(type_of_crossing), str(type_of_mutation), probability_of_mutation, size_of_population, parents, best_solutions, distance_matrix)
            new_population = gen.parents
            best_individual = min(gen.best_solutions)
            print(best_individual)
            for i in range(len(parents)):
                best.append(new_population[i].phenotype)
            print(min(best))
        else:
            condition = False

if __name__ == '__main__':
    main()