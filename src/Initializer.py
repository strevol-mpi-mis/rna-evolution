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
    
    def __init__(self, target_structure, population_size): 
        
        self.target_structure = target_structure 
        self.population_size = population_size 


    def init(self) : 
        pos = self.get_bp_position()
        nucluotides = ["A", "U", "G", "C"]
        base_paire = ["AU","UA","GU","GC","UG","CG"]
        arn = np.random.choice(nucluotides,len(self.target_structure))
        pop = []
        i = 0
        while i < self.population_size : 
            
            for bp_cord in pos : 
                bp = np.random.choice(base_paire,1)
                arn[bp_cord[0]] = bp[0][0]
                arn[bp_cord[1]] = bp[0][1]
            pop.append(''.join(arn))
            i = len(set(pop))
            
        return np.array(list(set(pop)))
    
    #Compute the fitness of an RNA Structure
    def fitness(self, structure) : 
        ref_xstrc = RNA.expand_Full(self.target_structure)
        xstrc = RNA.expand_Full(structure)
        
        return 1./(1.+RNA.tree_edit_distance(RNA.make_tree(ref_xstrc), RNA.make_tree(xstrc)))

    def initialize(self) : 
        n = 0 
        population = []
        nucluotides = ["A", "G", "U", "C"]
        init_depth = len(self.target_structure)
        for i in range(self.population_size):
            arn = numpy.random.choice(nucluotides,init_depth)
            seq = ''.join(arn)
            (strc, mfe) = RNA.fold(seq)

            ind = Individual(seq,strc,self.fitness(strc), mfe)
            population.append(ind)

        return population
