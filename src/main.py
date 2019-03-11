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


#Main function 
def main() : 
    
   
    population_size = 100
    number_of_generation = 300
    mut_rate = 0.4
    mut_bp = 0.5
    lamda = 1
    k= 15
    type_ = "ham"
    target_structure = "............((((......)))).....(((((((((........)).))))))).......((.(((((..(((....)))..)))))))......"
    init_deph = len(target_structure)
    print "Length ==", init_deph
    mut_prob = 1./(init_deph)
    mut_prob = 0.001

    mut_probs = np.array(rna.ptable(target_structure)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0

    landscape = Landscape.Landscape(lamda, k, type_, target_structure)
    rna_evolution = RNAEvolution.RNAEvolution(population_size,0,None,landscape)

    result = rna_evolution.run(number_of_generation,mut_probs,mut_bp,1)

    list_fitness = [] 
    good_result = []
    for ind in result : 
        print ind
        if ind.fitness !=0 : 
            list_fitness.append(ind.fitness)
            good_result.append(ind)
    list_fitness = sorted(list_fitness, reverse=True) 

    sorted_pop = [ ] 
    for fitness in list_fitness : 
        for ind in good_result : 
            if ind.fitness == fitness : 
                sorted_pop.append(ind)
    for i in range(10) : 
        print sorted_pop[i].RNA_seq, sorted_pop[i].mfe, sorted_pop[i].fitness
if __name__== "__main__" : 

    main()
