


import pandas 
import numpy as np 
import matplotlib.pyplot as plt 
import Individual


#Extract the fitness from the log file
def getFitness(population_size,log_file) : 
    
    fitnesses = [] 
    dataFrame = pandas.read_csv(log_file, sep=",")
    for ind in (dataFrame.values)[:population_size+1,1:] : 
        fitnesses.append(ind[3])
    
    return fitnesses

#Extract the fitness from the log file
def getIndividuals(population_size,log_file) : 
    
    population = [] 
    dataFrame = pandas.read_csv(log_file, sep=",")
    for ind in (dataFrame.values)[:population_size+1,1:] : 
        population.append(Individual.Individual(ind[0],ind[1],ind[3],ind[2]))
    
    return population

def hamming_distance(ind1, ind2) :
    distance = 0.
    for i in range(len(ind1.RNA_structure)) : 
            if ind1.RNA_structure[i] != ind2.RNA_structure[i] : 
                    distance +=1
                    
    return distance

#Novelty function 
def novelty(individual, population, k) : 
    list_metric = [] 

    for ind in population : 
        metric = hamming_distance(individual, ind)
        if metric == 0 : 
            continue
        list_metric.append(metric) 
    list_metric = sorted(list_metric)
    return sum(list_metric[:k])/k



def main() : 

    """
    pop = getIndividuals(109, "../Logs/NoveltyExp/News/Selection/0/0/gen1000.csv")
    saved_fitness = []
    sum_novelty = 0 

    saved_novelty = []

    for ind in pop : 
        saved_fitness.append(ind.fitness) 
        n = novelty(ind,pop, k = 15)
        saved_novelty.append(n)
    fitness_dist = []
    novelty_dist = []


    for i in range(len(pop)) : 
        novelty_dist.append(saved_novelty[i]/np.mean(saved_novelty))
        fitness_dist.append(pop[i].fitness/np.mean(saved_fitness))
    print fitness_dist 


    df = pandas.DataFrame({
        "Fitness" : fitness_dist , 
        "Novelty" : novelty_dist
    })

    df.plot.kde()
    plt.show()

    """

    fold = range(4)
    population_size = 109
    lamdas = [0]#,0.1,0.2,0.5,1]
    histoData = []
    methods = ["ED", "EDV", "F"]
    for n in  fold: 
            print "Loading folder ", n, "S....." 
            
            maxmaxs = [] 
            for meth in methods: 
                maxs = []
                for g in range(101) : 
                    dF = getFitness(population_size,"../Logs/Defect/66/"+str(n)+"/"+meth+"/gen"+str(g)+".csv")
                    maxs.append(np.mean(dF))
                maxmaxs.append(maxs)
            histoData.append(maxmaxs)

    histoData = np.array(histoData) 
    means = np.mean(histoData, axis=0)
    print means.shape
    
    for i in range(len(methods)): 
        plt.plot(means[i],label = methods[i])
    plt.legend(loc="lower right")
    
    plt.ylabel("Mean fitness fitness", fontweight="bold")
    plt.xlabel("Generation", fontweight="bold")
    #plt.savefig("../Logs/NoveltyExp/max-fitness-histo.eps")
    #plt.savefig("../Logs/NoveltyExp/max-fitness-histo.png")
    plt.show()

if __name__ == "__main__" : 
    main()