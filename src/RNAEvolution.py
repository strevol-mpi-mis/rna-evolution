'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA 
import random 
import Individual 
import Initializer
import matplotlib.pyplot as plt
import pandas
import os
import Logger
import pp


class RNAEvolution(object) : 

    def __init__ (self, population_size, lamda, archive, target_structure,  init_depth=10) : 

        self.population_size = population_size 
        self.init_depth = init_depth 
        self.lamda = lamda
        self.archive = archive
        self.initializer = Initializer.Initializer(target_structure, population_size)
    
    
    #Compute the fitness of an RNA Structure
    def tree_edit_fitness(self, structure) : 
        ref_xstrc = RNA.expand_Full(self.initializer.target_structure)
        xstrc = RNA.expand_Full(structure)
        
        return 1./(1.+RNA.tree_edit_distance(RNA.make_tree(ref_xstrc), RNA.make_tree(xstrc)))

    #Compute the fitness of an RNA Structure
    def base_paire_fitness(self, structure) : 
        return 1./(1.+RNA.bp_distance(self.initializer.target_structure,structure))


    #Mutation function 
    def mutateOne(self, individual, mut_p, mut_bp) :  
        base_paire = ["AU","UA","GU","GC","UG","CG"]
        nucluotides = ["A", "G", "U", "C"]

        RNA_seq = []
        for i in range(len(individual.RNA_seq)) :
            r = random.uniform(0,1)
        
            if r < mut_p[i] : 
                selct = numpy.random.choice(nucluotides,size=1)
                RNA_seq.append(selct[0])
            else : 
                RNA_seq.append(individual.RNA_seq[i])
        pos = individual.get_bp_position(self.initializer.target_structure)

        for bp_cord in pos : 
            r = random.uniform(0,1)
            if r < mut_bp : 
                bp = numpy.random.choice(base_paire,1, p=[0.1, 0.2, 0.2, 0.2, 0.2, 0.1])
                RNA_seq[bp_cord[0]] = bp[0][0]
                RNA_seq[bp_cord[1]] = bp[0][1]
        

        (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
        return Individual.Individual(''.join(RNA_seq), RNA_strc,self.base_paire_fitness(RNA_strc), mef)



    def mutateAll(self, population, mut_prob, mut_bp ) : 
        mutated_pop = [] 
        for individual in population : 
            mutated_pop.append(self.mutateOne(individual,mut_prob,mut_bp))

        return mutated_pop


    # Natural selection based on fitness proportionate method
    def fitness_proportion_selection(self, population, size) : 
        
        sum_fitness = 0 
        for ind in population : 
            sum_fitness += ind.fitness 

        proportion_prob = []
        for ind in population : 
            proportion_prob.append(ind.fitness/sum_fitness)
         
        selected = numpy.random.choice(population,size=size,p=proportion_prob)
        return selected


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

    def simple_EA(self, number_of_generation, mut_probs, log_folder, mut_bp) : 
    
        print (" Starting of evolution ")
        prev_population = self.initializer.initialize() #Initialize the population of RNA
        population_size = self.initializer.population_size
        n = number_of_generation
        
        logger = Logger.Logger(str(log_folder),str(self.lamda))
        logger.save_population(prev_population,0)
        
        while (n > 0) :
        
            print ('Generation '+str(number_of_generation - n))
            newgeneration = []
            selected_ind = self.fitness_proportion_selection(prev_population,population_size)
            newgeneration = self.mutateAll(selected_ind,mut_probs,mut_bp)
            
            prev_population = numpy.copy(newgeneration)
            n -=1
            logger.save_population(newgeneration,number_of_generation-n)


        return newgeneration

    
    def run(self, number_of_generation,mut_probs, mut_bp, nbjobs) : 
        #tuple of all parallel python servers to connect with
        ppservers = ()
        job_server = pp.Server(nbjobs, ppservers=ppservers)
        
        tasks = range(nbjobs)
        result = numpy.array([Individual.Individual("","",0,0)])
        #Parallel evolution for every lamda value
        print "Start running jog"
        jobs = [(task, job_server.submit(self.simple_EA, (number_of_generation, mut_probs,task,mut_bp,), ( self.fitness_proportion_selection,self.tree_edit_fitness, self.mutateAll, self.mutateOne,),
                                               ("numpy", "Individual", "RNAEvolution2", "RNA", "random","Logger","pandas","os","Initializer"))) for task in tasks]
        
        for task, job in jobs : 
            gen = job()
            result = numpy.insert(result, 0, numpy.array(gen))
        
        print "End of jobs"

        return result



