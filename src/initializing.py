'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy as np
import RNA as rna 
import random 
from Individual import Individual
import matplotlib.pyplot as plt
import pandas
import os
import RNAEvolution2
import RNAEvolution


def get_bp_position(ref_ss) : 
    position = rna.ptable(ref_ss)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
        if position[i] != 0 :
            base_paire_pos.append((i,position[i]-1))

    return base_paire_pos


def init(ref_ss, pop_size) : 
    pos = get_bp_position(ref_ss)
    nucluotides = ["A", "U", "G", "C"]
    base_paire = ["AU","UA","GU","GC","UG","CG"]
    arn = np.random.choice(nucluotides,len(ref_ss))
    pop = []
    i = 0
    while i < pop_size : 
        
        for bp_cord in pos : 
            bp = np.random.choice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]
        pop.append(''.join(arn))
        i = len(set(pop))
        
    return np.array(list(set(pop)))



#Main function 
def main() : 
    
    population_size = 500
    number_of_generation = 200
    init_depth = 40
    #mut_prob = 1./init_deph
    mut_prob = 0.01
    mut_rate = 0.5
    nucluotides = ["A", "G", "U", "C"]
    arn = np.random.choice(nucluotides,init_depth)
    rna_str = ''.join(arn)
    print ("rn", rna_str)
    
    # compute minimum free energy (MFE) and corresponding  referent structure
    (ref_ss, mfe) = rna.fold(rna_str)
    ref_ind = Individual(rna_str,ref_ss,0)
    
    data = [[ref_ind.RNA_seq, ref_ind.RNA_structure, RNAEvolution.fitness(ref_ind.RNA_structure,ref_ind.RNA_structure)]] 
    

    
    init_pop = init(ref_ind.RNA_structure,population_size)
    
    for seq in init_pop : 
        data.append([seq, rna.fold(seq)[0], RNAEvolution.fitness(rna.fold(seq)[0],ref_ind.RNA_structure)])


    dataFrame = pandas.DataFrame(data)

    print (dataFrame)

    dataFrame.to_csv("../Logs/init_data/data1000_new"+str(init_depth)+".csv")
    
      
    
   
if __name__== "__main__" : 

    main()