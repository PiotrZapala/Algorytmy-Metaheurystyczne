from random import choice, randint, uniform
import time
from pandas import array
import tsp_initialize as init
import tsp_functions as func

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