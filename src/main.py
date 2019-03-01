'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy as np
import RNA as rna 
import random 
from Individual import Individual
import matplotlib.pyplot as plt
import pandas
import os
import RNAEvolution
import RNAEvolution2


#Main function 
def main() : 
    
   
    population_size = 500
    number_of_generation = 100
    init_deph = 40
    mut_prob = 1./(init_deph)
    mut_rate = 0.4
    mut_bp = 0.5
    lamdas = [0,0.5,1]
    k= 5

    dataFrame = pandas.read_csv("../Logs/init_data/data1000_40.csv",sep=",")
    
    
    ref_ind = Individual((dataFrame.values)[0,1:][0],(dataFrame.values)[0,1:][1],(dataFrame.values)[0,1:][2])
    print ("RNA_ref = " + ref_ind.RNA_seq)

    init_pop = []
    for ind in (dataFrame.values)[1:population_size+1,1:] : 
        init_pop.append(Individual(ind[0],ind[1],ind[2]))
    
    print (population_size, len(init_pop))

    """
    for mut_prob in [0.025,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] : 
        print("Mut prob = ", mut_prob)
        RNAEvolution.evolution(ref_ind,init_pop,"Repport",population_size, number_of_generation, mut_prob,)
        print("=================================")
    """
    RNAEvolution2.evolution(ref_ind,init_pop,"Repport",population_size, number_of_generation, mut_prob,mut_rate)
        
    """
    meanHist = []
    for i in range(100) : 
        meanHist.append(pandas.read_csv("Logs/log_"+rna_str[:7]+"fitnessGen"+str(i+1)+".csv",sep=','))
    
    means = []
    for df in meanHist :
        means.append(((df.agg(['mean']).round(decimals=2)).values)[0,1])
    #meanDf = pandas.DataFrame(np.array(sorted(means))) 
    #fig, ax = plt.subplots()
    #meanDf.plot.kde(ax=ax,legend=False, title='Mean Evolution')
    plt.plot(means,title='Mean Evolution')
    plt.savefig("Logs/log_"+rna_str[:7]+".mean.eps")
    plt.show()
    """
if __name__== "__main__" : 

    main()