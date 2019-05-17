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

    def __init__ (self, population_size, lamda, archive, landscape) : 

        self.population_size = population_size 
        self.lamda = lamda
        self.archive = archive
        self.initializer = Initializer.Initializer(landscape, population_size)
        self.landscape = landscape

    
    #Mutation function 
    def mutateOne(self, individual, mut_p, mut_bp) :  
        #base_paire = ["GC","CG"]
        #nucluotides = ["A", "G", "C"]

        base_paire = ["AU","UA","GU","GC","UG","CG"]
        nucleotides = ["A", "G", "U", "C"]

        #base_paire = ["GC","CG"]
        #nucleotides = [ "A","G"]
        RNA_seq = []
        for i in range(len(individual.RNA_seq)) :
            r = random.uniform(0,1)
        
            if r < mut_p[i] : 
                selct = numpy.random.choice(nucleotides,size=1)
                RNA_seq.append(selct[0])
            else : 
                RNA_seq.append(individual.RNA_seq[i])
        pos = individual.get_bp_position(self.landscape.target_structure)

        for bp_cord in pos : 
            r = random.uniform(0,1)
            if r < mut_bp : 
                bp = numpy.random.choice(base_paire,1, p=[0.1,0.1,0.2,0.2,0.2,0.2])
                RNA_seq[bp_cord[0]] = bp[0][0]
                RNA_seq[bp_cord[1]] = bp[0][1]
        
        (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
        return Individual.Individual(''.join(RNA_seq), RNA_strc, mef,self.landscape.fitness(RNA_strc))



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

    def ComputeEnsDiversity(self, seq): 
        
        cmd = "echo "+ str(seq)+"| RNAfold -p |tail -1|cut -d ' ' -f 11 > rest.txt" 
        os.system(cmd)
        dt = pandas.read_csv('rest.txt', header=None)  
        return dt.values[0][0]

    def ComputeEnsDiversities(self, listOfSeq, task): 
        df = pandas.DataFrame(listOfSeq)
        df.to_csv("seq"+str(task)+".csv", header=False, index=False)
        cmd = "RNAfold -p --noPS < seq"+str(task)+".csv | grep -nhr 'diversity' | cut -d ' ' -f 11 > ensDiv"+str(task)
        os.system(cmd)
        dt = pandas.read_csv("ensDiv"+str(task), header=None)

        os.remove("seq"+str(task)+".csv")
        os.remove("ensDiv"+str(task))
        os.remove("dot.ps")
        return dt.values[:, 0]
    
    def ensDiversity_proportion_selection(self,population, size, task) : 
        listOfSeq = [ind.RNA_seq for ind in population]
        ensDiv = self.ComputeEnsDiversities(listOfSeq, task)
        ensDiv = numpy.array(ensDiv)
        selected = numpy.random.choice(population,size=size,p=numpy.array(ensDiv)/sum(ensDiv))
        return selected


    #################
    # Natural selection based on fitness proportionate and novelty proportion 
    def novelty_selection(self,population,size) : 

        saved_fitness = []
        sum_novelty = 0
        saved_novelty = [] 
        lamda = self.landscape.lamda
        k = self.landscape.k
        for ind in population : 
            saved_fitness.append(ind.fitness) 
            n = self.landscape.novelty(ind.RNA_structure, population)
            sum_novelty += n 
            saved_novelty.append(n)
        #self.archiving.archive(population, saved_novelty)

        proportion_prob = []
        for i in range(len(population)) : 
            s = (1-lamda)*numpy.divide((population[i].fitness-min(saved_fitness)),(max(saved_fitness)-min(saved_fitness))) + lamda*numpy.divide(float(saved_novelty[i]-min(saved_novelty)),(max(saved_novelty)-min(saved_novelty)))
            proportion_prob.append(s)
        proportion_prob = numpy.array(proportion_prob)
        choices = numpy.random.choice(population,size=size,p=proportion_prob/sum(proportion_prob))
        return choices

    
    def reproduce(self, population, size) : 
        list_fitness = [] 
        for ind in population : 
            list_fitness.append(ind.fitness)
        list_fitness = sorted(list_fitness, reverse=True) 

        sorted_pop = [ ] 
        for fitness in list_fitness : 
            for ind in population : 
                if ind.fitness == fitness : 
                    sorted_pop.append(ind)
        return sorted_pop[:size]
 

    ##################
    '''
        This method is used to cross two given genotypes to form two new offspring for the next generation.
        
        For two given sequences of bit, the crossover consists of selecting n-bit in genotype1 and n-bit in genotype2 
        and exchange them to form two new genotype. 

            INPUT
            =====
                    genotype1(Type: string) : a sequence of bit of lenth N. that represent a parent
                    genotype2(Type: string) : a sequence of bit of lenth N. that represent a donor
                    number_of_bit(Type: int): the default value is 1. It is the number bit to cross(to exchange)
            OUTPUT
            ======
                    The function returns two new offsprings.
    '''
    def crossover(self, genotype1, genotype2, number_of_bit=1) : 
        vect1 = map(str, genotype1.RNA_seq)
        vect2 = map(str, genotype2.RNA_seq)

        r = numpy.random.randint(0, len(vect1))
        swap = vect1[: r]
        vect1[:r] = vect2[:r]
        vect2[:r] = swap 

        (child1_structure, mef_child1) = RNA.fold(''.join(vect1))
        (child2_structure, mef_child2) = RNA.fold(''.join(vect2))
        child1 = Individual.Individual(''.join(vect1), child1_structure,self.landscape.fitness(child1_structure), mef_child1)
        child2 = Individual.Individual(''.join(vect2), child2_structure,self.landscape.fitness(child2_structure), mef_child2)
        
        return  child1, child2





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
    
        print (" Start of evolution ")
        #prev_population = self.initializer.initialize_from_csv("../Logs/NoveltyExp/News/old-select/"+str(log_folder)+"/"+str(self.landscape.lamda)+"/gen500.csv") #Initialize the population of RNA
        #prev_population = self.initializer.initialize_from_csv("../Logs/init_data/data1000_ham_30.csv") #Initialize the population of RNA
        
        prev_population = self.initializer.initialize() #Initialize the population of RNA
        population_size = self.initializer.population_size
        n = number_of_generation
        
        logger = Logger.Logger(str(log_folder),str(self.landscape.lamda))
        logger.bt_save_population(prev_population,prev_population,0)
        max_fitness = max([ind.fitness for ind in prev_population])
        while (n > 0) and (max_fitness < 1):
        
            if (number_of_generation - n)%100 == 0 : 
                print ('Generation '+str(number_of_generation - n)), "Max fitness = ", max_fitness
            newgeneration = []
            newgeneration =  self.reproduce(prev_population,int(0.1*population_size))
            
            if self.landscape.lamda >0 : 
                selected_ind = self.novelty_selection(prev_population,population_size)
            elif self.landscape.lamda == 0:
                selected_ind = self.fitness_proportion_selection(prev_population,population_size)
            else : 
                selected_ind = self.ensDiversity_proportion_selection(prev_population, population_size, log_folder)
            
            mutated = self.mutateAll(selected_ind,mut_probs,mut_bp)
            #selected_ind = numpy.insert(newgeneration, len(newgeneration),selected_ind)
            
            #newgeneration = numpy.insert(newgeneration, len(newgeneration),mutated)
            prev_population = numpy.copy(mutated)
            n -=1
            logger.bt_save_population(selected_ind, prev_population,number_of_generation-n)
            #max_fitness = max([ind.fitness for ind in prev_population])
            new_max = max([ind.fitness for ind in prev_population])
            if new_max > max_fitness : 
                max_fitness = new_max
            n -=1

        return newgeneration

    
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

    def EA_with_crossover(self, number_of_generation, mut_probs, log_folder, mut_bp) : 
    
        print (" Start of evolution ")
        #prev_population = self.initializer.initialize_from_csv("../Logs/init_data/data1000_ham_ 98.csv") #Initialize the population of RNA
        prev_population = self.initializer.initialize() #Initialize the population of RNA
        population_size = self.initializer.population_size
        n = number_of_generation
        
        logger = Logger.Logger(str(log_folder),str(self.landscape.lamda))
        logger.save_population(prev_population,0)
        max_fitness = max([ind.fitness for ind in prev_population])
        while (n > 0) and (max_fitness <1):
        
            print ('Generation '+str(number_of_generation - n))
            newgeneration = []
            newgeneration =  self.reproduce(prev_population,int(0.1*population_size))
            
            if self.landscape.lamda !=0 : 
                selected_ind = self.novelty_selection(prev_population,int(0.4*population_size))
            else:
                selected_ind = self.fitness_proportion_selection(prev_population,int(0.4*population_size))
            
            newgeneration = numpy.insert(newgeneration, len(newgeneration),self.mutateAll(selected_ind,mut_probs,mut_bp))
            cross_ind_num= int(0.6*population_size)
            while cross_ind_num > 0 : 
                if self.landscape.lamda !=0 : 
                    selected  = self.novelty_selection(prev_population,2)
                else:
                    selected = self.fitness_proportion_selection(prev_population,2)
            
                new_ind1, new_ind2 = self.crossover(selected[0],selected[1])

                numpy.insert(newgeneration, len(newgeneration),[new_ind1, new_ind2])
                cross_ind_num = cross_ind_num - 2


            prev_population = numpy.copy(newgeneration)
            n -=1
            logger.save_population(newgeneration,number_of_generation-n)
            max_fitness = max([ind.fitness for ind in prev_population])

        return newgeneration



    def run(self, number_of_generation,mut_probs, mut_bp, nbjobs, log_fold) : 
        #tuple of all parallel python servers to connect with
        ppservers = ()
        job_server = pp.Server(nbjobs, ppservers=ppservers)
        
        tasks = range(0, nbjobs)
        result = numpy.array([Individual.Individual("","",0,0)])
        #Parallel evolution for every lamda value
        print "Start running job"
        jobs = [(task, job_server.submit(self.simple_EA, (number_of_generation, mut_probs,task,mut_bp,), ( self.fitness_proportion_selection,self.mutateAll, self.mutateOne,),
                                               ("numpy", "Individual", "RNA", "random","Logger","pandas","os","Initializer", "Landscape"))) for task in tasks]
        
        for task, job in jobs : 
            gen = job()
            result = numpy.insert(result, 0, numpy.array(gen)[:100])
        
        print "End of jobs"

        return result


    def run_with_crossover(self, number_of_generation,mut_probs, mut_bp, nbjobs) : 
        #tuple of all parallel python servers to connect with
        ppservers = ()
        job_server = pp.Server(nbjobs, ppservers=ppservers)
        
        tasks = range(nbjobs)
        result = numpy.array([Individual.Individual("","",0,0)])
        #Parallel evolution for every lamda value
        print "Start running job"
        jobs = [(task, job_server.submit(self.EA_with_crossover, (number_of_generation, mut_probs,task,mut_bp,), ( self.fitness_proportion_selection,self.mutateAll, self.mutateOne,self.crossover,),
                                               ("numpy", "Individual", "RNA", "random","Logger","pandas","os","Initializer", "Landscape"))) for task in tasks]
        
        for task, job in jobs : 
            gen = job()
            result = numpy.insert(result, 0, numpy.array(gen)[:10])
        
        print "End of jobs"

        return result



