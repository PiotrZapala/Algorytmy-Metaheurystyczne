import math
import random
import tsplib95
from tsplib95 import distances
import numpy as np

def matrix(coord, distance):
    n = len(coord)
    distance_matrix = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            distance_matrix[i][j] = int(distance((coord[i][0], coord[i][1]), (coord[j][0], coord[j][1])))
    return n, distance_matrix


def symetric_random_instance(number_of_cities, min_distance, max_distance, seed):
    np.random.seed(seed)
    rand_matrix = np.random.randint(min_distance, max_distance + 1, size=(number_of_cities, number_of_cities))
    for i in range(number_of_cities):
        rand_matrix[i][i] = 0
        for j in range(number_of_cities):
            rand_matrix[j][i] = rand_matrix[i][j]
    return number_of_cities, rand_matrix


def asymetric_random_instance(number_of_cities, min_distance, max_distance, seed):
    np.random.seed(seed)
    rand_matrix = np.random.randint(min_distance, max_distance + 1, size=(number_of_cities, number_of_cities))
    for i in range(number_of_cities):
        rand_matrix[i][i] = 0
    return number_of_cities, rand_matrix


def euclidean_random_instance(number_of_cities, max_x, max_y, seed):
    distance = distances.euclidean
    np.random.seed(seed)
    coord = np.array([(np.random.randint(0, max_x, dtype=int), np.random.randint(0, max_y, dtype=int)) for _ in range(number_of_cities)])
    distance_matrix = np.zeros((number_of_cities, number_of_cities))
    for i in range(0, number_of_cities):
        for j in range(0, number_of_cities):
            distance_matrix[i][j] = int(distance((coord[i][0], coord[i][1]), (coord[j][0], coord[j][1])))
    return number_of_cities, distance_matrix, coord


def explicit_instance(filename, type):
    if type == 'tsp':
        instance = tsplib95.load('instance/tsp/' +filename)
    elif type =='atsp':
        instance = tsplib95.load('instance/atsp/' +filename) 
    edges = list(instance.get_edges())
    nodes = len(list(instance.get_nodes()))
    length = len(edges)
    weight = [None]*length
    matrix = np.zeros((nodes, nodes))
    for i in range(0, length):
        weight[i] = instance.get_weight(edges[i][0], edges[i][1])
        matrix[edges[i][0]-1][edges[i][1]-1] = weight[i]

    return nodes, matrix


def euclidean_instance(filename):
    instance = tsplib95.load('instance/tsp/' +filename)
    nodes =  len(list(instance.get_nodes()))
    edges =  list(instance.get_edges())
    length = len(edges)
    weight = [None]*length
    matrix = np.zeros((nodes, nodes))
    for i in range(0, length):
        weight[i] = instance.get_weight(edges[i][0], edges[i][1])
        matrix[edges[i][0]-1][edges[i][1]-1] = weight[i]

    return nodes, matrix


def read_file(filename, type):

    if type == 'tsp':
        f = open('instance/tsp/' +filename, 'r')
    elif type == 'atsp':
        f = open('instance/atsp/' +filename, 'r')

    name = f.readline().strip().split()[1]
    filetype = f.readline().strip().split()[1]
    comment = f.readline().strip().split()[1]
    dimension = f.readline().strip().split()[1] 

    line = f.readline()

    while line.find("EDGE_WEIGHT_TYPE") == -1:
        line = f.readline()

    if line.find("EUC_2D") != -1:
        number_of_cities, distance_matrix = euclidean_instance(filename)
    elif line.find("MAN_2D") != -1:
        distance = distances.manhattan
    elif line.find("MAX_2D") != -1:
        distance = distances.maximum
    elif line.find("GEO") != -1:
        distance = distances.geographical
    elif line.find("CEIL_2D") != -1:
        distance = distances.functools.partial(distances.euclidean, round=math.ceil)
    elif line.find("EXPLICIT") != -1:
        number_of_cities, distance_matrix = explicit_instance(filename, type)
    else:
        raise Exception

    return number_of_cities, distance_matrix

def initialization():
    while True:
        # chose instance
        print('1. import instance \n2. random instance')
        option = int(input().strip())
        if option == 1:
            # import instance from file - set optimum solution
            print("1. symetric \n2. asymetric")
            option = int(input().strip())
            if option == 1:
                print("Enter name of file: ")
                filename = input().strip()
                number_of_cities, distance_matrix = read_file(filename, type='tsp')
                #print(distance_matrix)
                                  
            elif option == 2:
                print("Enter name of file:")
                filename = input().strip()
                number_of_cities, distance_matrix = read_file(filename, type='atsp')
                print(distance_matrix)

            print("To continue press '0'")
            x = int(input().strip())
            if x == 0:
                break
            
        elif option == 2:
            # random instance (chose seed) - set optimum solution
            print("1. random euclidean 2D \n2. random symetric \n3. random asymetric")
            option = int(input().strip())
            if option == 1:
                print("Enter the number of cities:")
                n = int(input().strip())
                print("Enter max x:")
                max_x = int(input().strip())
                print("Enter max y:")
                max_y = int(input().strip())
                print("Enter seed:")
                seed = int(input().strip())
                number_of_cities, distance_matrix = euclidean_random_instance(n, max_x, max_y, seed)
                print(distance_matrix)

            elif option == 2:
                print("Enter the number of cities:")
                n = int(input().strip())
                print("Enter min distance:")
                min_distance = int(input().strip())
                print("Enter max distance:")
                max_distance = int(input().strip())
                print("Enter seed:")
                seed = int(input().strip())
                number_of_cities, distance_matrix = symetric_random_instance(n, min_distance, max_distance, seed)
                print(distance_matrix)

            elif option == 3:
                print("Enter the number of cities:")
                n = int(input().strip())
                print("Enter min distance:")
                min_distance = int(input().strip())
                print("Enter max distance:")
                max_distance = int(input().strip())
                print("Enter seed:")
                seed = int(input().strip())
                number_of_cities, distance_matrix = asymetric_random_instance(n, min_distance, max_distance, seed)
                print(distance_matrix)
            
            else:
                print("No choise - plis do it again.")

            print("To continue press '0'")
            x = int(input().strip())
            if x == 0:
                break
        else:
            print('wrong choise \nDo you want to exit? press 0 otherwise press 1')
            option = int(input().strip())
            if option == 0:
                break
            else:
                instance_exist = False

    return number_of_cities, distance_matrix