'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA 
import random 
import Individual 
import RNAEvolution
import Logger


#Create neutral subsets for a given population 
def clustering(pop, nb_cluster) : 
    result = []
    population = list(numpy.copy(pop))
    for i in range(nb_cluster-1) : 
        cluster = [population[0]]
        population.remove(population[0])
        i = 0
        while i < len(population)  : 

            if cluster[0].fitness == population[i].fitness : 
                cluster.append(population[i])
                population.remove(population[i])


            else : 
                i +=1
        result.append(cluster)
    result.append(population)
    return result

#Crossover function 
def crossOver(ref_strc, ind1, ind2) : 
    child1 = []
    child2 = []
    
    for i in range(len(ind1.RNA_seq)) : 
        child1.append(ind1.RNA_seq[i])
        child2.append(ind2.RNA_seq[i])
    single_point = random.randint(1,len(ind1.RNA_seq)-1)
    
    cross_part = child1[single_point:]
    child1[single_point:] = child2[single_point:]
    child2[single_point:] = cross_part

    (child1_strc, child1_mef) = RNA.fold(''.join(child1))
    (child2_strc, child2_mef) = RNA.fold(''.join(child2))

    new_ind1 = Individual.Individual(''.join(child1), child1_strc,RNAEvolution.fitness(ref_strc,child1_strc))
    new_ind2 = Individual.Individual(''.join(child2), child2_strc,RNAEvolution.fitness(ref_strc,child2_strc))
    
    return new_ind1, new_ind2

#Crossover function with cross probability
def crossOver2(ref_strc, ind1, ind2, cross_prob) : 
    child1 = []
    child2 = []
    
    for i in range(len(ind1.RNA_seq)) : 
        child1.append(ind1.RNA_seq[i])
        child2.append(ind2.RNA_seq[i])
    
    if random.uniform(0,1) < cross_prob :
        single_point = random.randint(1,len(ind1.RNA_seq)-1)
        cross_part = child1[single_point:]
        child1[single_point:] = child2[single_point:]
        child2[single_point:] = cross_part

    (child1_strc, child1_mef) = RNA.fold(''.join(child1))
    (child2_strc, child2_mef) = RNA.fold(''.join(child2))

    new_ind1 = Individual(''.join(child1), child1_strc,RNAEvolution.fitness(ref_strc,child1_strc))
    new_ind2 = Individual(''.join(child2), child2_strc,RNAEvolution.fitness(ref_strc,child2_strc))
    
    return new_ind1, new_ind2


def evolution(ref_ind, init_pop, log_folder, population_size, number_of_generation, mut_prob, mut_rate) : 
   

    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    
    n = number_of_generation
    i = 1
    histFitness = []

    #histFitness.append(sorted(getFitnesses(prev_population)))
    logger = Logger.Logger(str(log_folder),str("Exp4/")+str(mut_prob))
    logger.save_population(init_pop,0) 

    n = number_of_generation
      
    while n > 0 : 
        print ('Generation '+str(number_of_generation - n))
         

        newgeneration = []

        selected_ind = RNAEvolution.select(prev_population,int(mut_rate*population_size))
        newgeneration = RNAEvolution.mutateAll(ref_ind.RNA_structure,selected_ind,mut_prob)
        
        cross_ind_num = population_size-int(mut_rate*population_size)
        while cross_ind_num > 0 : 
            selected = RNAEvolution.select(prev_population,2)
            new_ind1, new_ind2 = crossOver(ref_ind.RNA_structure,selected[0],selected[1])
            
            newgeneration.append(new_ind1)
            newgeneration.append(new_ind2)
            cross_ind_num = cross_ind_num - 2


        prev_population = numpy.copy(newgeneration)


        n -=1
        logger.save_population(newgeneration,number_of_generation-n)
    return newgeneration

def evolution2(ref_ind, init_pop, number_of_generation, mut_prob, mut_rate) : 
   

    prev_population = init_pop
    population_size = len(init_pop)
    n = number_of_generation
    i = 1
    histFitness = []

    histFitness.append(sorted(RNAEvolution.getFitnesses(prev_population)))
      
    while number_of_generation > 0 : 
    
        print ('Generation '+str(i))
        i += 1 

        newgeneration = []

        selected_ind = RNAEvolution.select(prev_population,int(mut_rate*population_size))
        newgeneration = RNAEvolution.mutateAll(ref_ind.RNA_structure,selected_ind,mut_prob)
        
        cross_ind_num = population_size-int(mut_rate*population_size)
        while cross_ind_num > 0 : 
            selected = RNAEvolution.select(prev_population,2)
            new_ind1, new_ind2 = crossOver(ref_ind.RNA_structure,selected[0],selected[1])
            
            newgeneration.append(new_ind1)
            newgeneration.append(new_ind2)
            cross_ind_num = cross_ind_num - 2


        prev_population = numpy.copy(newgeneration)


        histFitness.append(sorted(RNAEvolution.getFitnesses(prev_population)))
        
        number_of_generation -=1

    return histFitness

def evolution_both_mut_cross(ref_ind, init_pop, number_of_generation, mut_prob, cross_prob) : 
   

    prev_population = init_pop
    population_size = len(init_pop)
    n = number_of_generation
    i = 1
    histFitness = []

    histFitness.append(sorted(RNAEvolution.getFitnesses(prev_population)))
      
    while number_of_generation > 0 : 
    
        print ('Generation '+str(i))
        i += 1 

        newgeneration = []

        selected_ind = RNAEvolution.select(prev_population,int(population_size))
        prev_population = RNAEvolution.mutateAll(ref_ind.RNA_structure,selected_ind,mut_prob)
        
        cross_ind_num = population_size
        while cross_ind_num > 0 : 
            selected = RNAEvolution.select(prev_population,2)
            new_ind1, new_ind2 = crossOver2(ref_ind.RNA_structure,selected[0],selected[1],cross_prob)
            
            newgeneration.append(new_ind1)
            newgeneration.append(new_ind2)
            cross_ind_num = cross_ind_num - 2


        prev_population = numpy.copy(newgeneration)


        histFitness.append(sorted(RNAEvolution.getFitnesses(prev_population)))
        
        number_of_generation -=1

    return histFitness
