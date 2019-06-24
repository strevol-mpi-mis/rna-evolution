'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy as np
import RNA as rna 
import random 
import matplotlib.pyplot as plt
import RNAEvolution
import Landscape
import pandas


#Main function 
def main() : 
    
   
    population_size = 100
    number_of_generation = 200
    mut_rate = 0.4
    mut_bp = 0.5
    lamda = -1
    k= 15
    type_ = "ham"
    method = "MDE"
    
    data = pandas.read_csv("defect-data.csv", sep=";")
    print len(data.values), len(list(set(data.values[:,0]))), data.values[:,0]
    for line in (data.values): 
        print "================================================================================"
        print line[0], line[1]
        print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        
        target_structure = line[1]
        init_deph = len(target_structure)
        print "Length ==", init_deph 
        mut_prob = 1./(init_deph) 
        #mut_prob = 0.001
        if init_deph > 300 : 
            number_of_generation = 150
        mut_probs = np.array(rna.ptable(target_structure)[1:])
        mut_probs = mut_probs + mut_prob
        mut_probs[mut_probs>mut_prob] = 0

        landscape = Landscape.Landscape(lamda, k, type_, target_structure)
        rna_evolution = RNAEvolution.RNAEvolution(population_size,0,None,landscape, method)

        result = rna_evolution.run(number_of_generation,mut_probs,mut_bp,4, str(line[0])+"/")

        list_fitness = [] 
        good_result = []
        for ind in result : 
            #print ind.RNA_structure, ind.fitness
            if ind.fitness !=0 : 
                list_fitness.append(ind.fitness)
                good_result.append(ind)
        list_fitness = sorted(list_fitness, reverse=True) 
        #print list_fitness
        sorted_pop = [ ] 
        for fitness in list_fitness : 
            for ind in good_result : 
                if ind.fitness == fitness : 
                    sorted_pop.append(ind)
        for i in range(10) : 
            print sorted_pop[i].RNA_seq, sorted_pop[i].mfe, sorted_pop[i].fitness
        
        print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    
if __name__== "__main__" : 

    main()
