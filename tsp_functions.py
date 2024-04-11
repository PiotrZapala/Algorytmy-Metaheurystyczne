import math
import random
import tsplib95
from tsplib95 import distances
import numpy as np
import matplotlib.pyplot as plt
import tsp_initialize as init

'''
Initialize TSP problem with set arguments
return nesesary variables for further calculation
'''
def initialize_problem():
    x, y = init.initialization()
    return x, y

def initialize_problem2():
    number_of_cities = 7
    min_distance = 10
    max_distance = 100
    seed = 0
    np.random.seed(seed)
    rand_matrix = np.random.randint(min_distance, max_distance + 1, size=(number_of_cities, number_of_cities))
    for i in range(number_of_cities):
        rand_matrix[i][i] = 0
        for j in range(number_of_cities):
            rand_matrix[j][i] = rand_matrix[i][j]
    return number_of_cities, rand_matrix

'''
Initialize random solution for TSP problem
return order of visited citis
'''
def random_solution(n):
    tour = list(range(n))
    random.shuffle(tour)
    return tour

'''
Initialize nearest neighbour solution for TSP problem
return order of visited citis
'''

def nearest_neighbor_extended_solution(n, distance_matrix):
    best_tour = nearest_neighbor_solution(n, 0, distance_matrix)
    for i in range(1,n):
        current_tour = nearest_neighbor_solution(n, i, distance_matrix)
        if get_weight(current_tour, distance_matrix) < get_weight(best_tour, distance_matrix):
            best_tour = current_tour
    return best_tour

def nearest_neighbor_solution(n, i, distance_matrix):
    unvisited = list(range(n))
    unvisited.remove(i)
    last = i
    tour = [i]
    while unvisited:
        nexxt = nearest(last, unvisited, distance_matrix)
        tour.append(nexxt)
        unvisited.remove(nexxt)
        last = nexxt
    return tour

def nearest(last, unvisited, distance_matrix):
    near = unvisited[0]
    min_distance = distance_matrix[last][near]
    for i in unvisited[1:]:
        if distance_matrix[last][i] < min_distance:
            near = i
            min_distance = distance_matrix[last][near]
    return near


def get_weight(tour, distance_matrix):
    weight = distance_matrix[tour[0]][tour[-1]]
    for i in range(1, len(tour)):
        weight += distance_matrix[tour[i]][tour[i - 1]]
    return weight

def swap_positions(tour, i, j):
    tour[i], tour[j] = tour[j], tour[i]
    return tour

def invert(tour, i, j):
    length_sub_list = j - i
    current_tour = tour.copy()
    for k in range(int(length_sub_list / 2 + 1)):
        current_tour = swap_positions(current_tour, i, j - k)
        i += 1
    return current_tour

def opt2(n, distance_matrix):
    best_solution = [None] * n
    current_solution = [None] * n
    results = []
    solution = random_solution(n)
    best_solution = solution
    min_weight = get_weight(best_solution, distance_matrix)
    improved = True
    while improved:
        i = 0
        while i < n - 1:
            j = i + 1
            while j < n:
                current_solution = best_solution.copy()
                current_solution = invert(current_solution, i, j)
                current_distance = get_weight(current_solution, distance_matrix)
                if current_distance < min_weight:
                    best_solution = current_solution
                    min_weight = current_distance
                    i = 0
                    j = 0
                results.append(min_weight)
                j += 1
            i += 1
        improved = False
    return best_solution

def opt22(n, distance_matrix):
    best_solution = [None] * n
    current_solution = [None] * n
    results = []
    solution = random_solution(n)
    best_solution = solution
    min_weight = get_weight(best_solution, distance_matrix)
    improved = True
    while improved:
        i = 0
        while i < (n - 1)/3:
            j = i + 1
            while j < n/3:
                current_solution = best_solution.copy()
                current_solution = invert(current_solution, i, j)
                current_distance = get_weight(current_solution, distance_matrix)
                if current_distance < min_weight:
                    best_solution = current_solution
                    min_weight = current_distance
                    i = 0
                    j = 0
                results.append(min_weight)
                j += 1
            i += 1
        improved = False
    return best_solution

def opt23(n, distance_matrix, child):
    best_solution = [None] * n
    current_solution = [None] * n
    results = []
    solution = child
    best_solution = solution
    min_weight = get_weight(best_solution, distance_matrix)
    improved = True
    while improved:
        i = 0
        while i < (n - 1)/2:
            j = i + 1
            while j < n/4:
                current_solution = best_solution.copy()
                current_solution = invert(current_solution, i, j)
                current_distance = get_weight(current_solution, distance_matrix)
                if current_distance < min_weight:
                    best_solution = current_solution
                    min_weight = current_distance
                    i = 0
                    j = 0
                results.append(min_weight)
                j += 1
            i += 1
        improved = False
    return best_solution
