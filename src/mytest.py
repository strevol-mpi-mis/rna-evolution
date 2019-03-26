import pandas 
import Individual
import random
import numpy 
import RNA
import collections
import pp 
import os
import platform
import Logger


def hamming_distance(ind1, ind2) :
    distance = 0.
    for i in range(len(ind1)) : 
        if ind1[i] != ind2[i] :
            distance +=1
                
    return distance

def fitness (ind,target) : 
    return 1/(1.+hamming_distance(target,ind))

#Mutation function 
def mutateOne(individual, mut_p, target) :  
    nucluotides = ["G", "C"]
    RNA_seq = []
    for i in range(len(individual.RNA_seq)) :
        r = random.uniform(0,1)
        if r < mut_p: 
            selct = numpy.random.choice(nucluotides,size=1)
            RNA_seq.append(selct[0])
        else : 
            RNA_seq.append(individual.RNA_seq[i])
    
    (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
    return Individual.Individual(''.join(RNA_seq), RNA_strc, mef, fitness(RNA_strc, target))


def mutateAll(population, mut_prob, target) : 
    mutated_pop = [] 
    for individual in population : 
        mutated_pop.append(mutateOne(individual,mut_prob,target))
    return mutated_pop

def novelty(individual, population, k) : 
    list_metric = [] 

    for ind in population : 
        metric = hamming_distance(individual, ind.RNA_structure)
        
        list_metric.append(metric) 
    #list_metric = sorted(list_metric)
    return numpy.mean(list_metric)


def novelty_selection(population,size,k) : 

   
    saved_novelty = [] 
    for ind in population : 
        n = novelty(ind.RNA_structure, population,k)
        saved_novelty.append(n)
    
    proportion_prob = numpy.array(saved_novelty)/sum(saved_novelty)
    choices = numpy.random.choice(population,size=size, p= proportion_prob)
    return choices


def save_population(population,gen, root_path) : 
    
    data = []
    for ind in population : 
        data.append([ind.RNA_seq, ind.RNA_structure,ind.mfe, ind.fitness])

    dataFrame = pandas.DataFrame(data)
    dataFrame.to_csv(root_path+"/gen"+str(gen)+".csv")
        
        
def fitness_proportion_selection(population, size) : 
    
    fitnesses = []
    for ind in population : 
        fitnesses.append(ind.fitness)

    selected = numpy.random.choice(population,size=size,p=numpy.array(fitnesses)/sum(fitnesses))
    return selected
    
def frequency_proportion_selection(population, size) : 
    
    structures = [ind.RNA_structure for ind in population]
    dic = collections.Counter(structures)
    saved_frequencies = []
    for ind in population : 
        saved_frequencies.append(dic[ind.RNA_structure])
    
    weights = 1/(1.+numpy.array(saved_frequencies))

    selected = numpy.random.choice(population,size=size,p=numpy.array(weights)/float(sum(weights)))
    return selected

def stochastic_selection(population, size) : 

    selected = numpy.random.choice(population,size=size)
    return selected


def simple_EA(target, number_of_generation, mut_probs, init_pop,k, selection_method, log_folder) : 
    
    print (" Starting of evolution ")
    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    
    root_path = "../Logs/MyTest2/"+str(log_folder)+"/"+selection_method
    try:
        os.makedirs(root_path)
    except OSError :
        print (" Can not initialize the log folder ")
    
    population_size =  len(init_pop)
    n = number_of_generation
    save_population(prev_population,0,root_path)
    maxfitness = max([ind.fitness for ind in prev_population])
    while (n > 0) and (maxfitness<1)   :
        
        if (number_of_generation - n)%100 == 0 : 
            print ('Generation '+str(number_of_generation - n))
        newgeneration = []
        if selection_method == "N": 
            selected_ind = novelty_selection(prev_population,population_size, k)
        elif selection_method == "F": 
            selected_ind = fitness_proportion_selection(prev_population,population_size)
        elif selection_method == "FREQ": 
            selected_ind = frequency_proportion_selection(prev_population,population_size)
        else : 
            selected_ind = stochastic_selection(prev_population,population_size)
        newgeneration = mutateAll(selected_ind,mut_probs, target)
        
        prev_population = numpy.copy(newgeneration)
        
        maxfitness = max([ind.fitness for ind in prev_population])
        n -=1
        save_population(prev_population, number_of_generation - n, root_path)

    return newgeneration, number_of_generation - n


def run(number_of_generation, mut_probs, k, log_folder) : 
    nucluotides = ["G","C"]
    init_pop = []
    target = "((....)).((....))"
    init_depth =len(target)
    
    while len(init_pop)<100 :
        arn = numpy.random.choice(nucluotides,init_depth)
        seq = ''.join(arn)
        (strc, mfe) = RNA.fold(seq)
        init_pop.append(Individual.Individual(seq, strc, mfe, fitness(strc, target)))

    #ppservers = ()
    #job_server = pp.Server(4, ppservers=ppservers)
    methods = ["N", "F", "S", "FREQ"]

    for method in methods : 
        simple_EA(target,number_of_generation, mut_probs,init_pop, k, method,log_folder)
    
    

def main() : 
    import sys
    mut_prob = 1/10.
    number_of_generation = 500
    k = 100
    ppservers = ()

    
    number_of_run  = int(sys.argv[1])
    job_server = pp.Server(4, ppservers=ppservers)
    
    print "Start running jog", number_of_run
    
    jobs = [(task , job_server.submit(run, (number_of_generation, mut_prob, k, task,), ( fitness_proportion_selection,mutateAll, mutateOne, frequency_proportion_selection, stochastic_selection, save_population,novelty, novelty_selection, fitness, hamming_distance, simple_EA),
                                            ("numpy", "Individual", "RNA", "random","Logger","pandas","os","Initializer", "Landscape", "collections"))) for task in range(number_of_run)]
    
    for task, job in jobs : 
        job()


if __name__ == "__main__":
    main()