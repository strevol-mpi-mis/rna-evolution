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
import subprocess

class RNAEvolution(object) : 

    def __init__ (self, population_size, lamda, archive, landscape, select_meth) : 

        self.population_size = population_size 
        self.lamda = lamda
        self.select_meth = select_meth
        self.archive = archive
        self.initializer = Initializer.Initializer(landscape, population_size)
        self.landscape = landscape

    
    #Mutation function 
    def mutateOne(self, individual, mut_probs, mut_bp) :  
        base_paire = ["GC","CG"]
        nucleotides = ["A", "G"]

        RNA_seq = numpy.array(list(individual.RNA_seq))
        r = numpy.random.rand(len(individual.RNA_seq))
        mut_pos =RNA_seq[r<mut_probs] 
        choices = numpy.random.choice(nucleotides, len(mut_pos),p=[0.2, 0.8])
        RNA_seq[r<mut_probs] = choices 
        pos = individual.get_bp_position(self.landscape.target_structure)

        for bp_cord in pos : 
            r = random.uniform(0,1)
            if r < mut_bp : 
                bp = numpy.random.choice(base_paire,1, p=[0.5, 0.5])
                RNA_seq[bp_cord[0]] = bp[0][0]
                RNA_seq[bp_cord[1]] = bp[0][1] 
        

        (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
        return Individual.Individual(''.join(RNA_seq), RNA_strc, mef,self.landscape.fitness(RNA_strc))



    def mutateAll(self, population, mut_prob, mut_bp ) : 
        mutated_pop = [self.mutateOne(ind,mut_prob,mut_bp) for ind in population] 
        return mutated_pop


    # Natural selection based on fitness proportionate method
    def fitness_proportion_selection(self, population, size) : 
        fitnesses = numpy.array([ind.fitness for ind in population])

        selected = numpy.random.choice(population,size=size,p=fitnesses/numpy.sum(fitnesses))
        return selected

    def ComputeEnsDiversity(self, seq): 
        
        cmd = "echo "+ str(seq)+"| RNAfold -p |tail -1|cut -d ' ' -f 11 > rest.txt" 
        os.system(cmd)
        dt = pandas.read_csv('rest.txt', header=None)  
        return dt.values[0][0]

    def ComputeEnsDiversities(self, listOfSeq, task): 
        df = pandas.DataFrame(listOfSeq)
        df.to_csv("seq"+str(task)+".csv", header=False, index=False)
        cmd = "RNAfold -p --noPS < seq"+str(task)+".csv | grep 'ensemble diversity' | cut -d ' ' -f 11 > ensDiv"+str(task)
        os.system(cmd)
        dt = pandas.read_csv("ensDiv"+str(task), header=None)

        os.remove("seq"+str(task)+".csv")
        os.remove("ensDiv"+str(task))
        #os.remove("dot.ps")
        return dt.values[:, 0]
    
    def ensDiversity_proportion_selection(self,population, size, task) : 
        listOfSeq = [ind.RNA_seq for ind in population]
        ensDiv = self.ComputeEnsDiversities(listOfSeq, task)
        ensDiv = numpy.array(ensDiv)
        selected = numpy.random.choice(population,size=size,p=numpy.array(ensDiv)/sum(ensDiv))
        return selected
    
    def ens_defect(self, sequence) : 
        p = subprocess.Popen("defect", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_input = sequence+"\n"+self.landscape.target_structure+"\n"
        #print cmd_input
        defect, error = p.communicate(cmd_input)
        defect = defect.split("\n")
        return 1/float(defect[-3])


    def ensDefect_proportion_selection(self, population, size) : 
        
        ensDefect = [self.ens_defect(ind.RNA_seq) for ind in population]
        selected = numpy.random.choice(population,size=size,p=numpy.array(ensDefect)/sum(ensDefect))
        return selected

    def min_ens_distance(self,sequence, min_energy=5.0) : 
        rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        rnasubopt_out, rnasubopt_err = rnasubopt.communicate(sequence) 

        cut_pipe = subprocess.Popen(['cut', '-f', '1', '-d', ' '], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cut_out, stderr = cut_pipe.communicate(rnasubopt_out)
        
        list_subopt =  cut_out.split()[1:] 

        return max([self.landscape.fitness(subopt) for subopt in list_subopt])


    def min_ens_distance_proportion_selection(self,population, size, target) : 
        
        ensDist = [self.min_ens_distance(ind.RNA_seq) for ind in population]
        selected = numpy.random.choice(population,size=size,p=numpy.array(ensDist)/sum(ensDist))
        return selected
    #################
    # Natural selection based on fitness proportionate and novelty proportion 
    def novelty_selection(self,population,size) : 

        saved_fitness = []
        saved_novelty = [] 
        lamda = self.landscape.lamda
        k = self.landscape.k
        for ind in population : 
            saved_fitness.append(ind.fitness) 
            n = self.landscape.novelty(ind.RNA_structure, population)
            saved_novelty.append(n)
        saved_novelty = numpy.array(saved_novelty) 
        saved_fitness = numpy.array(saved_fitness)
        proportion_prob =(1-lamda)*((saved_fitness-min(saved_fitness))/(max(saved_fitness)-min(saved_fitness))) + lamda*((saved_novelty-min(saved_novelty))/(max(saved_novelty)-min(saved_novelty)))
        
        choices = numpy.random.choice(population,size=size,p=proportion_prob/sum(proportion_prob))
        return choices

    
    def reproduce(self, population, size) : 
        list_fitness = []
        #ref = numpy.random.choice(a=[".",len(population[0].RNA_seq)])
        for ind in population : 
            list_fitness.append(ind.fitness)
        list_fitness = sorted(list_fitness, reverse=True) 

        sorted_pop = [ ] 
        for fitness in list_fitness : 
            for ind in population : 
                if ind.fitness == fitness : # or self.landscape.hamming_distance(ind.RNA_structure, ref) == 0 : 
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
        
        prev_population = self.initializer.init() #Initialize the population of RNA
        population_size = self.initializer.population_size
        n = number_of_generation
        
        logger = Logger.Logger(str(log_folder),str(self.select_meth))
        logger.save_population(prev_population,0)
        max_fitness = max([ind.fitness for ind in prev_population])
        newgeneration = numpy.copy(prev_population)
        while (n > 0) and (max_fitness < 1):
        
            if (number_of_generation - n)%100 == 0 : 
                print ('Generation '+str(number_of_generation - n)), "Max fitness = ", max_fitness
            newgeneration = []
            newgeneration =  self.reproduce(prev_population,int(0.1*population_size))
            
            if self.select_meth == "N" : 
                selected_ind = self.novelty_selection(prev_population,population_size)
            elif self.landscape.lamda == 0:
                selected_ind = self.fitness_proportion_selection(prev_population,population_size)
            elif self.select_meth == "MDE"  : 
                selected_ind = self.min_ens_distance(prev_population, population_size)
            
            mutated = self.mutateAll(selected_ind,mut_probs,mut_bp)
            #selected_ind = numpy.insert(newgeneration, len(newgeneration),selected_ind)
            
            newgeneration = numpy.insert(newgeneration, len(newgeneration),mutated)
            prev_population = numpy.copy(newgeneration)
            logger.save_population(prev_population,number_of_generation-n)
            #max_fitness = max([ind.fitness for ind in prev_population])
            new_max = max([ind.fitness for ind in prev_population])
            if new_max > max_fitness : 
                max_fitness = new_max
            n -=1

        return prev_population



    
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
        
        while (n > 0) :
        
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


        return newgeneration



    def run(self, number_of_generation,mut_probs, mut_bp, nbjobs, log_fold) : 
        #tuple of all parallel python servers to connect with
        ppservers = ()
        job_server = pp.Server(nbjobs, ppservers=ppservers)
        
        tasks = range(nbjobs)
        result = numpy.array([Individual.Individual("","",0,0)])
        #Parallel evolution for every lamda value
        print "Start running jog"
        jobs = [(task, job_server.submit(self.simple_EA, (number_of_generation, mut_probs,log_fold+str(task),mut_bp,), ( self.fitness_proportion_selection,self.mutateAll, self.mutateOne,),
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
        print "Start running jog"
        jobs = [(task, job_server.submit(self.EA_with_crossover, (number_of_generation, mut_probs,task,mut_bp,), ( self.fitness_proportion_selection,self.mutateAll, self.mutateOne,self.crossover,),
                                               ("numpy", "Individual", "RNA", "random","Logger","pandas","os","Initializer", "Landscape"))) for task in tasks]
        
        for task, job in jobs : 
            gen = job()
            result = numpy.insert(result, 0, numpy.array(gen)[:10])
        
        print "End of jobs"

        return result



