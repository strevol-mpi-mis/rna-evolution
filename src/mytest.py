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
import subprocess


def hamming_distance(ind1, ind2) :
    distance = 0.
    for i in range(len(ind1)) : 
        if ind1[i] != ind2[i] :
            distance +=1
                
    return distance

def fitness (ind,target) : 
    return 1/(1.+hamming_distance(target,ind))


#Mutation function 
def mutateOne(individual, mut_probs,target, mut_bp=0.5) :  
    base_paire = ["GC","CG"]
    nucleotides = ["A", "G"]

    RNA_seq = numpy.array(list(individual.RNA_seq))
    r = numpy.random.rand(len(individual.RNA_seq))
    mut_pos =RNA_seq[r<mut_probs] 
    choices = numpy.random.choice(nucleotides, len(mut_pos))
    RNA_seq[r<mut_probs] = choices 
    pos = individual.get_bp_position(target)

    for bp_cord in pos : 
        r = random.uniform(0,1)
        if r < mut_bp : 
            bp = numpy.random.choice(base_paire,1, p=[0.5, 0.5])
            RNA_seq[bp_cord[0]] = bp[0][0]
            RNA_seq[bp_cord[1]] = bp[0][1] 
    

    (RNA_strc, mef) = RNA.fold(''.join(RNA_seq))
    return Individual.Individual(''.join(RNA_seq), RNA_strc, mef,fitness(RNA_strc, target))


"""
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
"""

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


def bt_save_population(prev_pop, population,gen, root_path) : 
    data = []
    prev_data = []
    for i in range(len(population)) : 
        data.append([population[i].RNA_seq, population[i].RNA_structure,population[i].mfe, population[i].fitness])
        prev_data.append([prev_pop[i].RNA_seq, prev_pop[i].RNA_structure,prev_pop[i].mfe, prev_pop[i].fitness])

    
    dataFrame = pandas.DataFrame(data)
    prev_dataFrame = pandas.DataFrame(prev_data)
    prev_dataFrame.to_csv(root_path+"/prev_gen"+str(gen)+".csv")
    dataFrame.to_csv(root_path+"/gen"+str(gen)+".csv")


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

def ens_defect(sequence, structure) : 
    p = subprocess.Popen("defect", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmd_input = sequence+"\n"+structure+"\n"
    #print cmd_input
    defect, error = p.communicate(cmd_input)
    defect = defect.split("\n")
    return 1/float(defect[-3])


def ensDefect_proportion_selection(population, size, target) : 
    
    ensDefect = [ens_defect(ind.RNA_seq, target) for ind in population]
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDefect)/sum(ensDefect))
    return selected

def min_ens_distance(sequence, min_energy=1.0) : 
    rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rnasubopt_out, rnasubopt_err = rnasubopt.communicate(sequence) 

    cut_pipe = subprocess.Popen(['cut', '-f', '1', '-d', ' '], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cut_out, stderr = cut_pipe.communicate(rnasubopt_out)
    
    list_subopt =  cut_out.split()[1:] 

    return max([fitness(subopt) for subopt in list_subopt])


def min_ens_distance_proportion_selection(population, size, target) : 
    
    ensDist = [min_ens_distance(ind.RNA_seq) for ind in population]
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDist)/sum(ensDist))
    return selected

def my_ens_def(sequence, min_energy=1.0) : 
    rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rnasubopt_out, rnasubopt_err = rnasubopt.communicate(sequence) 
    result = numpy.array([[s.split()[0], s.split()[1]] for s in rnasubopt_out.split("\n")[1:-1]])
    d_s = numpy.array([fitness(sigma) for sigma in result[:,0]])
    kt = 0.612 
    z = sum(numpy.exp(-numpy.array(result[:,1], dtype=float)/kt))
    p = numpy.exp(-numpy.array(result[:,1], dtype=float)/kt)/z

    return sum(p*d_s)

def my_ens_def_proportion_selection(population, size, target) : 
    
    myED = [my_ens_def(ind.RNA_seq) for ind in population]
    selected = numpy.random.choice(population,size=size,p=numpy.array(myED)/sum(myED))
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

def ComputeEnsDiversity(seq): 
        
    cmd = "echo "+ str(seq)+"| RNAfold -p |tail -1|cut -d ' ' -f 11 > rest.txt" 
    os.system(cmd)
    dt = pandas.read_csv('rest.txt', header=None)  
    return dt.values[0][0]

def ComputeEnsDiversities(listOfSeq, task): 
    df = pandas.DataFrame(listOfSeq)
    df.to_csv("seq"+str(task)+".csv", header=False, index=False)
    cmd = "RNAfold -p --noPS < seq"+str(task)+".csv | grep 'ensemble diversity' | cut -d ' ' -f 11 > ensDiv"+str(task)
    #print cmd
    os.system(cmd)
    dt = pandas.read_csv("ensDiv"+str(task), header=None)

    os.remove("seq"+str(task)+".csv")
    os.remove("ensDiv"+str(task))
    #os.remove("dot.ps")
    return dt.values[:, 0]
    
def ensDiversity_proportion_selection(population, size, task) : 
    listOfSeq = [ind.RNA_seq for ind in population]
    ensDiv = ComputeEnsDiversities(listOfSeq, task)
    ensDiv = numpy.array(ensDiv)
    selected = numpy.random.choice(population,size=size,p=numpy.array(ensDiv)/sum(ensDiv))
    return selected



def simple_EA(target, number_of_generation, mut_probs, init_pop,k, selection_method, log_folder) : 
    
    print (" Starting of evolution ")
    prev_population = numpy.copy(init_pop) #Initialize the population of RNA
    
    root_path = "../Logs/Defect/66/"+str(log_folder)+"/"+selection_method
    try:
        os.makedirs(root_path)
    except OSError :
        print (" Can not initialize the log folder ")
    
    population_size =  len(init_pop)
    n = number_of_generation
    #save_population(prev_population,0,root_path)
    bt_save_population(prev_population,prev_population,0,root_path)

   #maxfitness = max([ind.fitness for ind in prev_population])
    while (n > 0) : #and (maxfitness<1)   :
        
        if (number_of_generation - n)%1 == 0 : 
            print ('Generation '+str(number_of_generation - n))
        newgeneration = []
        if selection_method == "N": 
            selected_ind = novelty_selection(prev_population,population_size, k)
        elif selection_method == "F": 
            selected_ind = fitness_proportion_selection(prev_population,population_size)
        elif selection_method == "MED": 
            selected_ind = min_ens_distance_proportion_selection(prev_population,population_size)
        elif selection_method == "ED": 
            selected_ind = ensDefect_proportion_selection(prev_population,population_size,target)
        elif selection_method == "EDV": 
            selected_ind = ensDiversity_proportion_selection(prev_population, population_size, log_folder)
        
        newgeneration = mutateAll(selected_ind,mut_probs, target)
        
        prev_population = numpy.copy(newgeneration)
        #maxfitness = max([ind.fitness for ind in prev_population])
        n -=1
        #save_population(prev_population, number_of_generation - n, root_path)
        bt_save_population(selected_ind, prev_population,number_of_generation-n,root_path)
       
    return newgeneration, number_of_generation - n


def run(number_of_generation, mut_probs, k, log_folder) : 
    
    target = "((....)).((....)).((....)).((....))"
    init_depth =len(target)
    pos = get_bp_position(target)
    nucleotides = ["A", "U", "G", "C"]
    base_paire = ["AU","UA","GU","GC","UG","CG"]
    pop = []
    i = 0
    while i < 100 : 
        
        if i < 4 : 
            arn = numpy.random.choice(nucleotides[i:i+1],init_depth)
        else : 
            arn = numpy.random.choice(nucleotides,len(target))

        for bp_cord in pos : 
            bp = numpy.random.choice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]
        pop.append(''.join(arn))
        i = len(set(pop))

    pop = numpy.array(list(set(pop)))
    import RNA as rna
    init_pop = []
    for seq in pop:
        strc,mfe = rna.fold(seq)
        init_pop.append(Individual.Individual(seq, strc, mfe, fitness(strc,target)))

    print len(init_pop)


    #ppservers = ()
    #job_server = pp.Server(4, ppservers=ppservers)
    methods = ["ED", "MED", "F"]

    for method in methods : 
        simple_EA(target,number_of_generation, mut_probs,init_pop, k, method,log_folder)
    
def get_bp_position(structure) : 
    position = RNA.ptable(structure)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
        if position[i] != 0 :
            if (position[i]-1,i) in base_paire_pos : 
                continue; 
            else : 
                base_paire_pos.append((i,position[i]-1))

    return base_paire_pos

def main() : 
    import sys

    target =  "((....)).((....)).((....)).((....))"
    init_depth =len(target)
    mut_prob = 1./init_depth
    number_of_generation = 100
    k = 100
    ppservers = ()

    mut_probs = numpy.array(RNA.ptable(target)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0


    
    number_of_run  = int(sys.argv[1])
    job_server = pp.Server(4, ppservers=ppservers)
    
    print "Start running jog", number_of_run
    
    jobs = [(task , job_server.submit(run, (number_of_generation, mut_probs, k, task,), (ComputeEnsDiversities,ensDiversity_proportion_selection, bt_save_population,fitness_proportion_selection,mutateAll,get_bp_position,ensDefect_proportion_selection, ens_defect,  mutateOne, frequency_proportion_selection, stochastic_selection, save_population,novelty, novelty_selection, fitness, hamming_distance, simple_EA),
                                            ("numpy", "Individual", "RNA", "random","Logger","pandas","os","Initializer", "Landscape", "collections","subprocess"))) for task in range(number_of_run)]
    
    for task, job in jobs : 
        job()


if __name__ == "__main__":
    main()
