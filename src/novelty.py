'''
    @author: Nono Saha Cyrille Merleau 
    @email: nonosaha@mis.mpg.de 

    This is a small distributed program that the goal is to run an RNA evolution using novelty search. 

    In this experiment we consider the novelty fitness as follow:

        N(x) = (1/k)sum(f(x,ui))  where f(x,ui) = 1/(1+tree_edit_dist(x,ui))

    Enjoy!!!!

'''

#Standard libraries
import sys
import os

#External libraries 
import numpy 
from Individual import Individual
import RNAEvolution2
import RNAEvolution
import Logger
import pp 
import RNA 
import pandas
import matplotlib.pyplot as plt
import random 
import Archive 



"""
    Selection function based on the novelty function defined as: 
         
            f_lambda = (1-lambda)*f(x,ref_struc) + lambda*N(x)

"""
# Natural selection based on Novelty metric
def select(population,size,lamda,k) : 
    
    sum_fitness = 0 
    sum_novelty = 0
    saved_novelty = []

    for ind in population : 
        sum_fitness += ind.fitness 
        n = novelty(ind,population,k)
        sum_novelty += n
        saved_novelty.append(n)

    proportion_prob = []
    for i in range(len(population)) : 
        proportion_prob.append((1-lamda)*(population[i].fitness/sum_fitness) + lamda*(saved_novelty[i]/sum_novelty))
    
    choices = numpy.random.choice(population,size=size,p=proportion_prob)
    return choices
 
def tournament_selection(population,tournament_size=10) :
        
    selected = random.sample(population,tournament_size)
    fitest = selected[0]

    for gen in selected : 
        if gen.fitness > fitest.fitness : 
            fitest = gen
    return fitest


#Novelty of all indisvidual 
def novelties(population, k, lamda, log_folder ) : 
    list_strc = []
    for ind in population : 
        list_strc.append(ind.RNA_structure)
            

    dataFrame= pandas.DataFrame(list_strc)
    dataFrame.to_csv(str(log_folder)+str(lamda)+"structures",sep=" ",index=False, header=False)
            
    cmd = "RNAdistance -Xm < "+str(log_folder)+str(lamda)+"structures >" +str(log_folder)+str(lamda)+".csv"
            
    os.system(cmd) 

    result = open(str(log_folder)+str(lamda)+".csv","r")
    j = 1
    lines = result.readlines()
    Lower = []
    while j < len(lines) : 
        line =numpy.array(lines[j].split(), dtype = float)  
        if len(line) >0 : 
            Lower.append(line)
        j+=1
    os.remove(str(log_folder)+str(lamda)+".csv")
    os.remove(str(log_folder)+str(lamda)+"structures")

    matrix_distance = buildSymetricMethod(Lower)
    
    novelties = []
    for line in matrix_distance :
        sorted_line = sorted(line) 
        #sorted_line = numpy.array(sorted_line)
        #sorted_line = sorted_line[sorted_line>10]
        novelties.append(sum(sorted_line[:k])/k)

    return novelties

# Natural selection based on Novelty metric
def optimal_select(population,size,lamda,k, log_folder)  : 
    
    sum_fitness = 0 
    saved_novelty = []

    for ind in population : 
        sum_fitness += ind.fitness 

    proportion_prob = []
    saved_novelty = novelties(population,k,lamda, log_folder)
    
    archive.archive(population, saved_novelty)

    sum_novelty = sum(saved_novelty)
    for i in range(len(population)) : 
        proportion_prob.append((1-lamda)*(population[i].fitness/sum_fitness) + lamda*(saved_novelty[i]/sum_novelty))
    
    choices = numpy.random.choice(population,size=size,p=proportion_prob)
    return choices

def optimal_select_with_archving(population,size,lamda,k, log_folder, archive)  : 
    
    sum_fitness = 0 
    saved_novelty = []

    for ind in population : 
        sum_fitness += ind.fitness 

    proportion_prob = []
    if lamda !=0 : 
        saved_novelty = novelties(population,k,lamda, log_folder)
        archive.archive(population, saved_novelty)
        print "Arhive size ==========", len(archive.archiving)
        sum_novelty = sum(saved_novelty)
        for i in range(len(population)) : 
            proportion_prob.append((1-lamda)*(population[i].fitness/sum_fitness) + lamda*(saved_novelty[i]/sum_novelty))
    else : 
        for i in range(len(population)) : 
            proportion_prob.append(population[i].fitness/sum_fitness)
    
    choices = numpy.random.choice(population,size=size,p=proportion_prob)
    return choices
#Novelty function 
def novelty(individual, population, k) : 
    list_metric = [] 

    for ind in population : 
        metric = RNA.tree_edit_distance(RNA.make_tree(RNA.expand_Full(individual.RNA_structure)),RNA.make_tree(RNA.expand_Full(ind.RNA_structure)))
        if metric == 0 : 
            continue
        list_metric.append(metric) 
    list_metric = sorted(list_metric)
    return sum(list_metric[:k])/k


#Reconstruct the symetric matrix distance from the lower triangular matrix 
def buildSymetricMethod(lower) : 
    result = []
    result.append(numpy.zeros(len(lower)+1))

    for m in lower : 
        c = numpy.zeros(len(lower)+1)
        c[:len(m)] = m
        result.append(c)
    result = numpy.array(result)
    return result + result.T

def reproduce(population, size) : 
    list_fitness = [] 
    for ind in population : 
        list_fitness.append(ind.fitness)
    list_fitness = sorted(list_fitness) 

    sorted_pop = [ ] 
    for fitness in list_fitness : 
        for ind in population : 
            if ind.fitness == fitness : 
                sorted_pop.append(ind)
    return sorted_pop[:size]

def adv_evolution_with_elitism(ref_ind, init_pop, number_of_generation, mut_probs, lamda, k, log_folder, mut_bp) : 
   
    print (" Starting of evolution ")
    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    population_size = len(init_pop)
    n = number_of_generation
    
    logger = Logger.Logger(str(log_folder),str(lamda))
    logger.save_population(init_pop,0)
    maxfitness = max([ind.fitness for ind in prev_population])  
    while (n > 0)and (maxfitness<1)  :
    
        print ('Generation '+str(number_of_generation - n))
        newgeneration = []
        newgeneration =  reproduce(prev_population,int(0.1*population_size))
        selected_ind = optimal_select(prev_population,population_size,lamda,k, log_folder)
        newgeneration =numpy.insert(newgeneration, len(newgeneration), RNAEvolution.adv_mutateAll(ref_ind.RNA_structure,selected_ind,mut_probs,mut_bp))
        
        prev_population = numpy.copy(newgeneration)
        maxfitness = max([ind.fitness for ind in prev_population])
        n -=1
        logger.save_population(newgeneration,number_of_generation-n)


    return newgeneration


def adv_evolution_with_archiving(ref_ind, init_pop, number_of_generation, mut_probs, lamda, k, log_folder, mut_bp, archive_size) : 
   
    print (" Starting of evolution ")
    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    population_size = len(init_pop)
    n = number_of_generation
    
    logger = Logger.Logger(str(log_folder),str(lamda))
    logger.save_population(init_pop,0)
    maxfitness = max([ind.fitness for ind in prev_population]) 

    archive = Archive.Archive(archive_size, "N")

    while (n > 0)and (maxfitness<1)  :
    
        print ('Generation '+str(number_of_generation - n))
        newgeneration = []
        newgeneration =  reproduce(prev_population,int(0.1*population_size))
        selected_ind = optimal_select_with_archving(numpy.insert(prev_population,len(prev_population),archive.archiving)
        ,population_size,lamda,k, log_folder, archive)
        newgeneration =numpy.insert(newgeneration, len(newgeneration), RNAEvolution.adv_mutateAll(ref_ind.RNA_structure,selected_ind,mut_probs,mut_bp))
        
        prev_population = numpy.copy(newgeneration)
        maxfitness = max([ind.fitness for ind in prev_population])
        n -=1
        print "Size of the archive =============== ",len(archive.archiving)
        logger.save_population(newgeneration,number_of_generation-n)


    return newgeneration




def evolution(ref_ind, init_pop, number_of_generation, mut_probs, lamda, k, log_folder, mut_bp) : 
   
    print (" Starting of evolution ")
    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    population_size = len(init_pop)
    n = number_of_generation
    
    logger = Logger.Logger(str(log_folder),str(lamda))
    logger.save_population(init_pop,0)
    maxfitness = max([ind.fitness for ind in prev_population])  
    while (n > 0)and (maxfitness<1)  :
    
        print ('Generation '+str(number_of_generation - n))

        newgeneration = []
        selected_ind = optimal_select(prev_population,population_size,lamda,k, log_folder)
        newgeneration = RNAEvolution.adv_mutateAll(ref_ind.RNA_structure,selected_ind,mut_probs,mut_bp)
        prev_population = numpy.copy(newgeneration)
        maxfitness = max([ind.fitness for ind in prev_population])
        n -=1
        logger.save_population(newgeneration,number_of_generation-n)


    return newgeneration

def evolution_with_crossover(ref_ind, init_pop, number_of_generation, mut_prob,lamda, k, mut_rate,log_folder) : 
   

    prev_population = numpy.copy(init_pop) #Initialize the population of RNA 
    population_size = len(init_pop)
    n = number_of_generation
    i = 1
    
    logger = Logger.Logger(str(log_folder),str(lamda))
    logger.save_population(init_pop,0)
    while n > 0 : 
    
        print ('Generation ' +str(i))
        i += 1 

        newgeneration = []

        selected_ind = select(prev_population,int(mut_rate*population_size),lamda,k)
        newgeneration = RNAEvolution.mutateAll(ref_ind.RNA_structure,selected_ind,mut_prob)
        
        cross_ind_num = population_size-int(mut_rate*population_size)
        while cross_ind_num > 0 : 
            selected = select(prev_population,2,lamda,k)
            new_ind1, new_ind2 = RNAEvolution2.crossOver(ref_ind.RNA_structure,selected[0],selected[1])
            
            newgeneration.append(new_ind1)
            newgeneration.append(new_ind2)
            cross_ind_num = cross_ind_num - 2


        prev_population = numpy.copy(newgeneration)
        n -=1
        logger.save_population(newgeneration,number_of_generation-n)

    return newgeneration


def combined_evolution(ref_ind, init_pop, number_of_generation, mut_prob, k, mut_rate,log_folder) : 
   

    n = int(number_of_generation/2)
  
    prev_population = evolution_with_crossover(ref_ind,init_pop,n,mut_prob,0,k,mut_rate,log_folder)

    newgeneration = evolution_with_crossover(ref_ind,prev_population,n,mut_prob,0.5,k,mut_rate,log_folder)

    return newgeneration


def main(): 
    #tuple of all parallel python servers to connect with
    ppservers = ()
    #ppservers = ("10.0.0.1",)
    
    if len(sys.argv) > 1:
        ncpus = int(sys.argv[1])
        folder_number = str(sys.argv[2])

       
        # Creates jobserver with ncpus workers
        job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers=ppservers)


    print "Starting pp with", job_server.get_ncpus(), "workers"
    print "Folder Number == ", folder_number
    population_size = 80
    number_of_generation = 60
    init_deph = 40
    mut_prob = 1./(init_deph)
    mut_rate = 0.4
    mut_bp = 1./(init_deph)
    lamdas = [0, 0.3]
    k= 10
    archive_size = 5

    dataFrame = pandas.read_csv("../Logs/init_data/data1000_40.csv",sep=",")
    
    
    ref_ind = Individual((dataFrame.values)[0,1:][0],(dataFrame.values)[0,1:][1],(dataFrame.values)[0,1:][2])
    print ("RNA_ref = " + ref_ind.RNA_seq)
    mut_probs = numpy.array(RNA.ptable(ref_ind.RNA_structure)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0

    init_pop = []
    for ind in (dataFrame.values)[1:population_size+1,1:] : 
        init_pop.append(Individual(ind[0],ind[1],ind[2]))
    
    print (population_size, len(init_pop))

    #print (novelties(init_pop,15,1))



    
    #Parallel evolution for every lamda value
    print "Start running jog"
    jobs = [(lamda, job_server.submit(adv_evolution_with_archiving, (ref_ind, init_pop, number_of_generation, mut_probs,lamda, k,folder_number,mut_bp, archive_size), (RNAEvolution.initialize, optimal_select, optimal_select_with_archving, RNAEvolution.fitness, RNAEvolution.mutateAll, RNAEvolution.mutateOne, RNAEvolution.getFitnesses, Individual, reproduce, novelties,buildSymetricMethod),
                                               ("numpy", "Individual", "RNAEvolution", "RNAEvolution2", "RNA", "random","Logger","pandas","os", "Archive"))) for lamda in lamdas]
    print "End of jobs"
    
    
    for lamda, job in jobs:
        job()
    """
    combined_evolution(ref_ind, init_pop, number_of_generation, mut_prob, k,mut_rate,folder_number,)
    """
if __name__== "__main__": 
    
    main()