from random import choice, randint, uniform
import time
from pandas import array
import tsp_initialize as init
import tsp_functions as func
import population as po
import individual 

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
        def __ocover(parent1, parent2, start, stop):
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
            child1_genotype = __ocover(self.selected_fathers[i].genotype, self.selected_mothers[i].genotype, start, stop)
            child2_genotype = __ocover(self.selected_mothers[i].genotype, self.selected_fathers[i].genotype, start, stop)
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