


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

    fold = range(20)
    population_size = 109
    lamdas = [1]#,0.1,0.2,0.5,1]
    histoData = []
    for n in  fold: 
            print "Loading folder ", n, "S....." 
            
            maxmaxs = [] 
            for i in range(1) : 
                maxs = []
                for g in range(1000,1001) : 
                    dF = getFitness(population_size,"../Logs/NoveltyExp/News/old-select/"+str(n)+"/"+str(lamdas[i])+"/gen"+str(g)+".csv")
                    maxs.append(max(dF))
                maxmaxs.append(max(maxs))
            histoData.append(maxmaxs)

    histoData = np.array(histoData) 
    """
    success = [0,0,0,0,0]
    print histoData
    for data in histoData :
        if data[0] == 1 : 
            success[0] +=1 
        if data[1] == 1 : 
            success[1] +=1 
        if data[2] ==1 : 
            success[2] +=1

        if data[3] ==1 : 
            success[3] +=1

        if data[4] ==1 : 
            success[4] +=1

    print success
    """
    plt.boxplot(histoData, labels=[r"$1$"])
    
    plt.ylabel("Max fitness", fontweight="bold")
    plt.xlabel("Lambda", fontweight="bold")
    #plt.savefig("../Logs/NoveltyExp/max-fitness-histo.eps")
    #plt.savefig("../Logs/NoveltyExp/max-fitness-histo.png")
    plt.show()

if __name__ == "__main__" : 
    main()