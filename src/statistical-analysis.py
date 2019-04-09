'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy as np
import os
import RNA as rna 
import random 
import Logger
from numpy.polynomial.polynomial import polyfit
from Individual import Individual
import matplotlib.pyplot as plt
import pandas
import os
from scipy import interpolate
#import RNAEvolution2
from mpl_toolkits.mplot3d import Axes3D
#import novelty 






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
        population.append(Individual(ind[0],ind[1],ind[3],ind[2]))
    
    return population

def hamming_dist (ind1, ind2)  : 
   distance = 0.
   for i in range(len(ind2)): 
	if ind1[i] != ind2[i] :
	    distance +=1 
   return distance

#Novelty function 
def novelty(individual, population, k) : 
    list_metric = [] 

    for ind in population : 
       # metric = RNA.tree_edit_distance(RNA.make_tree(RNA.expand_Full(individual.RNA_structure)),RNA.make_tree(RNA.expand_Full(ind.RNA_structure)))
	metric = hamming_dist(individual.RNA_structure, ind.RNA_structure)

        if metric == 0 : 
            continue
        list_metric.append(metric) 
    list_metric = sorted(list_metric)
    return sum(list_metric[:k])/k

#
#Main function 
def main() : 
   
    population_size = 200
    lamdas = [0,1]
    meanHist = []
    mut_cross_probs = []

    means = []
    colors = {
        0: "g",
        0.1: "b",
        0.2: "r"    }
    
    i = 1

    """
    #Last Mean Fitness Distribution
    last_mean0 = []
    last_mean1 = []
    for n in [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20]: 
        for lamda in  lamdas : 
            df = pandas.read_csv("../Logs/NoveltyExp/l40_100g/"+str(n)+"S/Job"+str(lamda)+"/Means"+str(lamda)+".csv", sep=",")
            list_ = (df.values).T[1]
            if lamda == 0 :
                last_mean0.append(list_[-1])
            else : 
                last_mean1.append(list_[-1])
    #plt.plot(last_mean0)      
    #plt.plot(last_mean1)  
    
    #fig, ax = plt.subplots()
    
    dataFrame = pandas.DataFrame({
        'lambda = 0' : last_mean0,
        'lambda = 0.5' : last_mean1
    })

    dataFrame.plot.kde()
    
    print dataFrame.median()
    
    #plt.title(r"$(1-\lambda)f_0 + \lambda N(x)$")
    #plt.legend(loc='upper left', shadow=True, fontsize='12')  
    
    #plt.xlabel("Generation")
    #plt.ylabel("Mean Fitness")
    #plt.savefig("NoveltyVsFitness_Stat9.eps")
    #plt.savefig("NoveltyVsFitness_Stat9.png")
    plt.show()
    """

    """
    #Max Fitness plots of all experiments
    last_max0 = []
    last_max1 = []
    j = 1
    for n in range(1,21): 
        print "Loading folder ", n, "S....."
        for lamda in  lamdas : 
            maxs = []
            for i in range(0,201) : 
                fitnesses = getFitness(population_size,"../Logs/NoveltyExp/l80_200g/"+str(n)+"S/"+str(lamda)+"/gen"+str(i)+".csv")
                maxs.append(max(fitnesses))
            if lamda == 0 :
                last_max0.append(max(maxs))
            elif lamda == 0.5 : 
                last_max1.append(max(maxs))

    dataFrame = pandas.DataFrame({
        'lambda = 0' : last_max0,
        'lambda = 0.5' : last_max1
    })

    dataFrame.plot.kde()
    
    print dataFrame.median()
    plt.show()
    
    
    
    #Mean fitness Plots of 20 experiments
    j = 1
    for n in range(1,21):
        print ("Loading folder " + str(n)+ "S.....") 
        for lamda in  lamdas : 
            means= []
            for i in range(0,301) : 
                dF = getFitness(population_size,"../Logs/NoveltyExp/Cross_Mut80-9/"+str(n)+"S/"+str(lamda)+"/gen"+str(i)+".csv")
                means.append(np.mean(dF))
            if j > 2 : 
                plt.plot(means,color=colors[lamda])
            else :
                plt.plot(means,color=colors[lamda],label=r"$\lambda = $"+str(lamda))
        
            j = j +1
    plt.xlabel("Generation")
    plt.ylabel("Mean Fitness")
    plt.title(r"$(1-\lambda)f_0 + \lambda N(x)$")
    plt.legend(loc='upper left', shadow=True, fontsize='12')
    plt.savefig("../Logs/NoveltyExp/Cross_Mut80-9/Mean.eps")
    plt.savefig("../Logs/NoveltyExp/Cross_Mut80-9/Mean.png")
    plt.show()
    

    """
    # Plot on the jumps after each 10.000 generations  
    #Max fitness Plots of 20 experiments
    
    j = 1
    histoData = []
    successes = []
    fold = range(0,50)
    gen = [100,200,300]#,6000,7000,8000,9000,10000]
    #fold.remove(6);fold.remove(7);fold.remove(8);fold.remove(9);fold.remove(10)
    result = []
    for g  in gen:
        print "Loading folder ", g, "S....." 
        maxmax = []
        success = [0,0]
        for n in fold  : 
            maxs = []
            for lamda in lamdas : 
                dF = getFitness(population_size,"../Logs/NoveltyExp/News/old-select/"+str(n)+"/"+str(lamda)+"/gen"+str(g)+".csv")
                maxs.append(max(dF))
            #print maxs
            for i in range(len(success)): 
                if maxs[i]>0.33: 
                    success[i] = success[i] + 1
            maxmax.append(maxs)
        

            """
            if j > 3 : 
                plt.plot(maxs,colors[lamda])
            else :
                plt.plot(maxs,colors[lamda],label=r"$\lambda = $"+str(lamda))
            j = j +1
            """
        successes.append(success)
        print np.array(maxmax).shape
        result.append(np.median(np.array(maxmax), axis=0))
        #maxmax.append(max(maxs))
        histoData.append(maxmax)
   
    result = np.array(result)
    print result[:,-1]
    print "Number of success = ", successes 
    #histoData = np.array(histoData) 
    #print histoData.shape 
    #medians = np.median(histoData, axis=0)
    #print medians 
    #print medians.shape 
    #df = pandas.read_csv("maxs.csv", sep=",")
    #plt.show()
    
    

    #histoData = np.array(histoData)
    
    #print histoData.shape
    #print histoData , histoData[:,0] , histoData[:,1]
    #histoData = histoData[:,0] - histoData[:,1]
    #print histoData
    
    plt.plot(gen,result[:,0],'-bs',label=r"$\lambda = $"+str(lamdas[0]))
    #plt.plot(gen,result[:,1],'-gs',label=r"$\lambda = $"+str(lamdas[1]))
    #plt.plot(gen,result[:,2],'-r+',label=r"$\lambda = $"+str(lamdas[2]))
    #plt.plot(gen,result[:,3],'-b+',label=r"$\lambda = $"+str(lamdas[3]))
    plt.plot(gen,result[:,-1],'-rs',label=r"$\lambda = $"+str(lamdas[-1]))

    plt.legend(loc="upper left", fontsize='12')
    plt.show()
    """
    
    
    
    plt.hist(histoData,color='b')
    
    plt.ylabel("Density")
    plt.xlabel(r"$(Fitness_{max} - FitnesNovelty_{max})$")
    plt.savefig("../Logs/NoveltyExp/Cross_Mut80-10/stat1.eps")
    plt.savefig("../Logs/NoveltyExp/Cross_Mut80-10/stat1.png")
    histoData = np.array(histoData)
    nn = len(histoData[histoData>0])
    npp = len(histoData[histoData<0])
    print nn, npp

    plt.show()
    """
    """
    #plt.boxplot(histoData)
    plt.xlabel("generation")
    plt.ylabel("Max Fitness")
    plt.legend(loc='upper left', shadow=True, fontsize='12')
    plt.savefig("../Logs/NoveltyExp/News/15000-20.eps")
    plt.savefig("../Logs/NoveltyExp/News/15000-20.png")
    plt.show()
    """
    """
    #Last Mean fitness Plots of 20 experiments
    last_mean0 = []
    last_mean1 = []
    j = 1
    for n in range(1,21):
        print "Loading folder ", n, "S....." 
        for lamda in  lamdas : 
            means= []
            for i in range(0,201) : 
                dF = getFitness(population_size,"../Logs/NoveltyExp/l80_200g/"+str(n)+"S/"+str(lamda)+"/gen"+str(i)+".csv")
                means.append(np.mean(dF))
            if lamda == 0 :
                last_mean0.append(means[-1])
            else : 
                last_mean1.append(means[-1])
        
            j = j +1
    dataFrame = pandas.DataFrame({
        'lambda = 0' : last_mean0,
        'lambda = 0.5' : last_mean1
    })

    dataFrame.plot.kde()
    
    print dataFrame.median()
    #plt.savefig("NoveltyVsFitnessMax_Stat9.eps")
    #plt.savefig("NoveltyVsFitnessMax_Stat9.png")
    plt.show()
    """
    """
        
    i = 0 
    n = 1 
    pop  = getIndividuals(population_size,"../Logs/NoveltyExp/l80_200g/"+str(n)+"S/"+str(lamda)+"/gen"+str(i)+".csv")
    list_strc = []
    
    for ind in pop : 
        list_strc.append(ind.RNA_structure)
    

    dataFrame= pandas.DataFrame(list_strc)
    dataFrame.to_csv("structures",sep=" ",index=False, header=False)
    
    cmd = "RNAdistance -Xm < structures >result.csv"
    
    os.system(cmd) 

    result = open("result.csv","r")
    n = 1
    lines = result.readlines()
    sum1 = np.array([])
    while n < len(lines) : 
        line =np.array(lines[n].split(), dtype = int)  
        sum1 = np.concatenate([line,sum1])
        n +=1
    print (sum(sum1)*2 )/(200*199)

    """
    """
    #Max Novelty Plots of 20 experiments
    j = 1
    for n in range(1,2):
        print "Loading folder ", n, "S....." 
        maxs = []
        for i in range(0,201) : 
            list_strc = []
            print n
            pop  = getIndividuals(population_size,"../Logs/NoveltyExp/l80_200g/"+str(n)+"S/"+str(lamda)+"/gen"+str(i)+".csv")
            for ind in pop : 
                list_strc.append(ind.RNA_structure)
            

            dataFrame= pandas.DataFrame(list_strc)
            dataFrame.to_csv("structures",sep=" ",index=False, header=False)
            
            cmd = "RNAdistance -Xm < structures >result.csv"
            
            os.system(cmd) 

            result = open("result.csv","r")
            j = 1
            lines = result.readlines()
            sum1 = np.array([])
            while j < len(lines) : 

                line =np.array(lines[j].split(), dtype = int)  
                sum1 = np.concatenate([line,sum1])
                j+=1
            print (sum(sum1)*2 )/(200*199)
            maxs.append((sum(sum1)*2 )/(200*199))
            print "i == ", i
        
        plt.plot(maxs)
        
    plt.xlabel("Generation")
    plt.ylabel("Max Fitness")
    plt.title(r"$(1-\lambda)f_0 + \lambda N(x)$")
    #plt.legend(loc='upper left', shadow=True, fontsize='12')
    #plt.savefig("NoveltyVsFitnessMax_Stat9.eps")
    #plt.savefig("NoveltyVsFitnessMax_Stat9.png")
    plt.show()

    """

    """
    #Max Novelty Plots of 20 experiments
    j = 1
    lamdas = [0,0.5,1]
    for lamda in lamdas:
        print "Loading folder ", lamda, "S....." 
        maxs = []
        for i in range(0,101) : 
            list_strc = []
            pop  = getIndividuals(population_size,"../Logs/NoveltyExp/Cross_Mut80/"+str(lamda)+"/gen"+str(i)+".csv")
            for ind in pop : 
                list_strc.append(ind.RNA_structure)
            

            dataFrame= pandas.DataFrame(list_strc)
            dataFrame.to_csv("structures",sep=" ",index=False, header=False)
            
            cmd = "RNAdistance -Xm < structures >result.csv"
            
            os.system(cmd) 

            result = open("result.csv","r")
            j = 1
            lines = result.readlines()
            sum1 = np.array([])
            while j < len(lines) : 
                line =np.array(lines[j].split(), dtype = int)  
                sum1 = np.concatenate([line,sum1])
                j+=1
            print (sum(sum1)*2 )/(200*199)
            maxs.append((sum(sum1)*2 )/(200*199))
            print "i == ", i
        
        plt.plot(maxs,label=r"$\lambda = $"+str(lamda))
        
    plt.xlabel("Generation")
    plt.ylabel("Mean Novelty")
    plt.title(r"$(1-\lambda)f_0 + \lambda N(x)$")
    plt.legend(loc='upper left', shadow=True, fontsize='12')
    #plt.savefig("NoveltyVsFitnessMax_Stat9.eps")
    #plt.savefig("NoveltyVsFitnessMax_Stat9.png")
    plt.show()
    """
    """
    #List files 
    j = 1
    histoData = []
    
    for n in range(1,51):
        print "Loading folder ", n, "S....." 
        nb_gen = []
        for lamda in  lamdas : 
            maxs = []
            files = os.listdir("../Logs/NoveltyExp/Cross_Mut40-1/"+str(n)+"S/"+str(lamda)+"/")
            nb_gen.append(len(files)- 1) 
        histoData.append(nb_gen)

    histoData = np.array(histoData)
    list1 = []
    list2 = []
    list3 = []

    print histoData 
    for n in  range(1,51): 
        maxs = [] 
        for i in range(3) : 
            dF = getFitness(population_size,"../Logs/NoveltyExp/Cross_Mut40-1/"+str(n)+"S/"+str(lamdas[i])+"/gen"+str(histoData[n-1][i])+".csv")
            maxs.append(max(dF))
        list1.append(maxs[0])
        list2.append(maxs[1])
        list3.append(maxs[2])
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax = plt.subplot(121)
    plt.boxplot(histoData, labels=[r"$\lambda=0$", r"$\lambda=0.5$", r"$\lambda=1$"])
    plt.title("Maxinum number of generation", fontweight="bold")
    plt.ylabel("Generation", fontweight="bold")
    
    plt.savefig("../Logs/NoveltyExp/Cross_Mut40-1/cross-mut40-box.eps")
    plt.savefig("../Logs/NoveltyExp/Cross_Mut40-1/cross-mut40-box.png")
    #df = pandas.read_csv("maxs.csv", sep=",")
    #histoData = np.array((df.values)[:,1:])
    plt.show()
    #print histoData.shape
    #print histoData , histoData[:,0] , histoData[:,1]
    #histoData = histoData[:,0] - histoData[:,1]
    #print histoData
    
    #for i in range(3) : 
    #    plt.plot(sorted(histoData[:,i]),color=colors[lamdas[i]],label=r"$\lambda = $"+str(lamdas[i]))
    
    #ax.set_aspect(1)
    
    """
    
    j = 1
    histoData = []
    """
    for n in range(1,51):
        print "Loading folder ", n, "S....." 
        nb_gen = []
        for lamda in  lamdas : 
            maxs = []
            files = os.listdir("../Logs/NoveltyExp/News/Mut80-9-1500/"+str(n)+"S/"+str(lamda)+"/")
            nb_gen.append(len(files)- 1) 
        histoData.append(nb_gen)

    histoData = np.array(histoData)
    """
    """
    list1 = []
    list2 = []
    list3 = []
    fold = range(1,51)
    #print fold
    #fold.remove(6);fold.remove(7);fold.remove(8);fold.remove(9)
    histoData = np.zeros((5, 1001))
    for n in  fold: 
        print "Loading folder ", n, "S....." 
        maxs = []
        maxmaxs = [] 
        for i in range(5) : 
            maxs = []
            for g in range(1001) : 
                dF = getFitness(population_size,"../Logs/NoveltyExp/07.02.19/1/"+str(n)+"S/"+str(lamdas[i])+"/gen"+str(g)+".csv")
                maxs.append(max(dF))
            
            maxmaxs.append(maxs)
        maxmaxs = np.array(maxmaxs) 
        print maxmaxs.shape 
        histoData = histoData + maxmaxs 

      
    #plt.subplot(122)

    histoData =histoData/len(fold)
    #plt.boxplot(histoData, labels=[r"$0$", r"$0.1$", r"$0.2$",r"$0.3$", r"$0.4$", r"$0.5$",r"$0.6$", r"$0.7$", r"$0.8$",r"$0.9$", r"$1$"])
    #plt.title("mutation prob 0.025 with elitism")
    #plt.ylabel("Generation")

    #plt.savefig("../Logs/NoveltyExp/News/Mut80-9-1500/box.eps")
    #plt.savefig("../Logs/NoveltyExp/News/Mut80-9-1500/box.png")
    #plt.show()
    """
    """
    histoData = np.array(list1) -  np.array(list2)
    plt.hist(histoData,color='b')
    
    plt.ylabel("Frequence", fontweight="bold")
    plt.xlabel(r"$(Fitness_{max} - FitnesNovelty_{max})$", fontweight="bold")
    plt.savefig("../Logs/NoveltyExp/Cross_Mut80-11/cross-mut80-11-stat1.eps")
    plt.savefig("../Logs/NoveltyExp/Cross_Mut80-11/cross-mut80-11-stat1.png")
    histoData = np.array(histoData)
    nn = len(histoData[histoData>0])
    npp = len(histoData[histoData<0])
    print nn, npp
    """
    """
    for i in range(5) : 
        plt.plot(histoData[i],label=r"$\lambda = $"+str(lamdas[i]))
    plt.legend(loc='upper left', shadow=True, fontsize='12')
    plt.ylabel("Max fitness", fontweight="bold")
    plt.xlabel("Generation", fontweight="bold")
    plt.savefig("../Logs/NoveltyExp/07.02.19/1/max-fitness-all.eps")
    plt.savefig("../Logs/NoveltyExp/07.02.19/1/max-fitness-all.png")
    plt.show()
    """

    """
    #List files 
    j = 1
    
    mut_probs = [0.0125,0.025, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    means = []
    for mut in mut_probs : 
        histoData = []
        print "Loading folder ", mut, "....." 
        for n in range(1,51):
            nb_gen = []
            for lamda in  lamdas : 
                maxs = []
                files = os.listdir("../Logs/NoveltyExp/Mut_rate_turning-1/"+str(mut)+"/"+str(n)+"S/"+str(lamda)+"/")
                nb_gen.append(len(files)- 1) 
            histoData.append(nb_gen)

        histoData = np.array(histoData)
        

        print histoData, np.mean(histoData, axis=0)
        means.append(np.mean(histoData, axis=0))
    means = np.array(means)
    print means[:,0]
    plt.plot(mut_probs,means[:,0], '-rs', label=r"$\lambda = 0$")
    plt.plot(mut_probs,means[:,1], '-gs', label=r"$\lambda = 0.5$")
    plt.plot(mut_probs,means[:,2], '-bs', label=r"$\lambda = 1$")
    plt.ylabel("Generation", fontweight="bold")
    plt.xlabel("Mutation prob",fontweight="bold")
    plt.legend(loc='lower right', shadow=True, fontsize='12')
    plt.title("Mean generation respect to mutation prob", fontweight="bold")
    plt.savefig("../Logs/NoveltyExp/Mut_rate_turning-1/summary-mean.eps")
    plt.savefig("../Logs/NoveltyExp/Mut_rate_turning-1/summary-mean.png")
    plt.show()
    
    """
    """
    histoData = []
    
    for n in range(1,51):
        print "Loading folder ", n, "S....." 
        maxmax = []
        for lamda in  lamdas : 
            maxs = []
             
            dF = getFitness(population_size,"../Logs/NoveltyExp/News/Mut80-5/"+str(n)+"S/"+str(lamda)+"/gen1000.csv")
            maxs.append(max(dF))
            
        maxmax.append(maxs)
    
    
    #df = pandas.read_csv("maxs.csv", sep=",")
    
    histoData = np.array(maxmax)
    
    print histoData.shape
    #print histoData , histoData[:,0] , histoData[:,1]
    #histoData = histoData[:,0] - histoData[:,1]
    #print histoData
    
    #for i in range(3) : 
    #    plt.plot(sorted(histoData[:,i]),color=colors[lamdas[i]],label=r"$\lambda = $"+str(lamdas[i]))
    
    plt.plot(histoData[:,0],'-rs', label=r"$\lambda = 0$")
    plt.plot(histoData[:,1],'-gs', label=r"$\lambda = 0.5$")
    plt.plot(histoData[:,2],'-bs', label=r"$\lambda = 1$")
    
    plt.ylabel("Last Max Fitnes")
    plt.xlabel("Generation")
    #plt.savefig("../Logs/NoveltyExp/Cross_Mut80-10/stat1.eps")
    #plt.savefig("../Logs/NoveltyExp/Cross_Mut80-10/stat1.png")
    
    plt.show()

    """
    """
    pop = getIndividuals(109,"../Logs/NoveltyExp/News/4SS/1/gen30000.csv")
    sum_fitness = 0
    sum_novelty = 0 
    save_novelty = []

    for ind in pop : 
         sum_fitness += ind.fitness 
         n = novelty(ind,pop,k=15)
         save_novelty.append(n)
    fitness_dist = [] 
    novelty_dist = []

    for i in range(len(pop)) :
        novelty_dist.append(save_novelty[i]/np.mean(save_novelty))
        fitness_dist.append(pop[i].fitness/sum_fitness)

    print fitness_dist

    df = pandas.DataFrame({
        "Fitness":fitness_dist,
        "Novelty":novelty_dist})
    df.plot.kde()

    #plt.hist(fitness_dist)
    plt.show()
    """

if __name__== "__main__" : 

    main()
