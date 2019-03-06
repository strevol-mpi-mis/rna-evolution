'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA 
import random 
from Individual import Individual
import matplotlib.pyplot as plt
import pandas
import os
import RNAEvolution


class Initializer(object) : 
    
    def __init__(self, landscape, population_size): 
        
        self.population_size = population_size 
        self.landscape = landscape


    def init(self) : 
        pos = self.get_bp_position()
        nucluotides = ["A", "U", "G", "C"]
        base_paire = ["AU","UA","GU","GC","UG","CG"]
        arn = numpy.random.choice(nucluotides,len(self.landscape.target_structure))
        pop = []
        i = 0
        while i < self.population_size : 
            
            for bp_cord in pos : 
                bp = numpy.random.choice(base_paire,1)
                arn[bp_cord[0]] = bp[0][0]
                arn[bp_cord[1]] = bp[0][1]
            pop.append(''.join(arn))
            i = len(set(pop))
            
        return numpy.array(list(set(pop)))
    
    

    def initialize(self) : 
        n = 0 
        population = []
        nucluotides = ["A", "G", "U", "C"]
        init_depth = len(self.landscape.target_structure)
        for i in range(self.population_size):
            if i < 4 : 
                arn = numpy.random.choice(nucluotides[i:i+1],init_depth)
                seq = ''.join(arn)
            else : 
                arn = numpy.random.choice(nucluotides,init_depth)
                seq = ''.join(arn)
            (strc, mfe) = RNA.fold(seq)

            ind = Individual(seq,strc,self.landscape.fitness(strc), mfe)
            population.append(ind)

        return population
