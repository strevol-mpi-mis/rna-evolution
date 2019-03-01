'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA 
import random 
import Individual 
import matplotlib.pyplot as plt
import pandas
import os
import Logger

#Get the base paire position 
def get_bp_position(ref_ss) : 
    position = RNA.ptable(ref_ss)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
        if position[i] != 0 :
            base_paire_pos.append((i,position[i]-1))

    return base_paire_pos


#Mutation function 
def mutateOne(ref_strc, ind, mut_p) :  
    
    nucluotides = ["A", "G", "U", "C"]
    RNA_seq = []
    for i in range(len(ind.RNA_seq)) :
        r = random.uniform(0,1)
    
        if r < mut_p : 
            selct = numpy.random.choice(nucluotides,size=1)
            RNA_seq.append(selct[0])
        else : 
            RNA_seq.append(ind.RNA_seq[i])
    
    (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
    return Individual.Individual(''.join(RNA_seq), RNA_strc,fitness(ref_strc,RNA_strc))

#Advance mutation function 
def adv_mutateOne(ref_strc, ind, mut_p, mut_bp) :  
    base_paire = ["AU","UA","GU","GC","UG","CG"]
    nucluotides = ["A", "G", "U", "C"]

    RNA_seq = []
    for i in range(len(ind.RNA_seq)) :
        r = random.uniform(0,1)
    
        if r < mut_p[i] : 
            selct = numpy.random.choice(nucluotides,size=1)
            RNA_seq.append(selct[0])
        else : 
            RNA_seq.append(ind.RNA_seq[i])
    pos = get_bp_position(ref_strc)

    for bp_cord in pos : 
        r = random.uniform(0,1)
        if r < mut_bp : 
            bp = numpy.random.choice(base_paire,1)
            RNA_seq[bp_cord[0]] = bp[0][0]
            RNA_seq[bp_cord[1]] = bp[0][1]
    

    (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
    return Individual.Individual(''.join(RNA_seq), RNA_strc,fitness(ref_strc,RNA_strc))



def mutateAll(ref_strc, population, mut_prob ) : 
    mutated_pop = [] 
    for individual in population : 
        mutated_pop.append(mutateOne(ref_strc,individual,mut_prob))

    return mutated_pop

def adv_mutateAll(ref_strc, population, mut_prob, mut_bp ) : 
    mutated_pop = [] 
    for individual in population : 
        mutated_pop.append(adv_mutateOne(ref_strc,individual,mut_prob,mut_bp))

    return mutated_pop


#Compute the fitness of an RNA Structure
def fitness(ref_strc, structure) : 
    ref_xstrc = RNA.expand_Full(ref_strc)
    xstrc = RNA.expand_Full(structure)
    
    return 1./(1.+RNA.tree_edit_distance(RNA.make_tree(ref_xstrc), RNA.make_tree(xstrc)))


# Natural selection based on fitness proportionate method
def select(population,size) : 
    
    sum_fitness = 0 
    for ind in population : 
        sum_fitness += ind.fitness 

    proportion_prob = []
    for ind in population : 
        proportion_prob.append(ind.fitness/sum_fitness)
    
    choices = numpy.random.choice(population,size=size,p=proportion_prob)
    return choices


#Generate the initial population of RNA of size N 
'''
    This method is used to genrate a N-random RNA sequence as an initial population.
        INPUT
        =====
                
        OUTPUT
        ======
                The function returns a list of RNA individuals.
'''
def initialize(ref_strc,population_size, init_depth=50) : 
    n = 0 
    population = []
    nucluotides = ["A", "G", "U", "C"]
   
    for i in range(population_size):
        arn = numpy.random.choice(nucluotides,init_depth)
        seq = ''.join(arn)
        (strc, mfe) = RNA.fold(seq)

        ind = Individual.Individual(seq,strc,fitness(ref_strc,strc))
        population.append(ind)

    return population

#Transform a list of RNA sequence to a list of RNA structure 
def toRNAstructure(list_of_RNA_sequences) : 
    result = [] 
    for seq in list_of_RNA_sequences : 
        result.append(RNA.fold(seq)[0])
    return result 

'''
 This function is implementing the simple genetic algorithm
        INPUT
        =====
                population_size(Type: int) : the number of RNA sequences to generate.
                number_of_generation(Type: int) : the number of generation.
                mut_prob(Type: float between [0,1]): the mutation probability.
        OUTPUT
        ======
                The function returns a list of individual of the new generation.
'''

def evolution(ref_ind, init_pop,log_folder, population_size, number_of_generation, mut_prob) : 
   

    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    
    n = number_of_generation
    i = 1
    histFitness = []

    histFitness.append(sorted(getFitnesses(prev_population)))
    logger = Logger.Logger(str(log_folder),str("Exp3/")+str(mut_prob))
    logger.save_population(init_pop,0) 

    n = number_of_generation

    while n > 0 : 
    
        print ('Generation '+str(number_of_generation - n))
        
        newgeneration = []

        selected_ind = select(prev_population,population_size)
        newgeneration = mutateAll(ref_ind.RNA_structure,selected_ind,mut_prob)
        prev_population = numpy.copy(newgeneration)

        n -=1
        logger.save_population(newgeneration,number_of_generation-n)
    return newgeneration


def getFitnesses(population) : 
    result = [] 
    for ind in population : 
        result.append(ind.fitness)
    return result

