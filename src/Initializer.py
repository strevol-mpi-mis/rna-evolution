'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA 
import random 
from Individual import Individual
import matplotlib.pyplot as plt
import pandas
import os
import RNAEvolution
import Landscape


class Initializer(object) : 
    
    def __init__(self, landscape, population_size): 
        
        self.population_size = population_size 
        self.landscape = landscape


    def init(self) : 
        pos = self.get_bp_position()
        nucluotides = ["A", "U", "G", "C"]
        base_paire = ["AU","UA","GU","GC","UG","CG"]
        arn = numpy.random.choice(nucluotides,len(self.landscape.target_structure))
        pop = []
        i = 0
        while i < self.population_size : 
            
            for bp_cord in pos : 
                bp = numpy.random.choice(base_paire,1)
                arn[bp_cord[0]] = bp[0][0]
                arn[bp_cord[1]] = bp[0][1]
            pop.append(''.join(arn))
            i = len(set(pop))
            
        return numpy.array(list(set(pop)))
    
    

    def initialize(self) : 
        n = 0 
        population = []
        nucluotides = ["G", "C", "A", "U"]
        init_depth = len(self.landscape.target_structure)
        for i in range(self.population_size):
            if i < 4 : 
                arn = numpy.random.choice(nucluotides[i:i+1],init_depth)
                seq = ''.join(arn)
            else : 
                arn = numpy.random.choice(nucluotides,init_depth)
                seq = ''.join(arn)
            (strc, mfe) = RNA.fold(seq)

            ind = Individual(seq,strc, mfe,self.landscape.fitness(strc))
            population.append(ind)

        return population

    def initialize_from_csv(self, path) : 
        
        dataFrame = pandas.read_csv(path,sep=",")
    
        ref_ind = Individual((dataFrame.values)[0,1:][0],(dataFrame.values)[0,1:][1],(dataFrame.values)[0,1:][3],(dataFrame.values)[0,1:][2])
        print ("RNA_ref = " + ref_ind.RNA_structure)
        
        init_pop = []
        for ind in (dataFrame.values)[1:self.population_size+1,1:] : 
            init_pop.append(Individual(ind[0],ind[1],ind[3],ind[2]))

        return init_pop

    #Main function 
def main() : 
        
    population_size = 1000
    number_of_generation = 500
    mut_rate = 0.4
    mut_bp = 0.5
    lamda = 0.3
    k= 15
    type_ = "ham"
    target_structure = ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................."
    init_depth = len(target_structure)
    print "Length ==", init_depth
    mut_prob = 1./(init_depth)
    #mut_prob = 0.001

    mut_probs = numpy.array(RNA.ptable(target_structure)[1:])
    mut_probs = mut_probs + mut_prob
    mut_probs[mut_probs>mut_prob] = 0

    landscape = Landscape.Landscape(lamda, k, type_, target_structure)
        
    data = [["NA", ".....((...((....))...((....))...((....))...((....))...((....))...((....))...))...................."
    ,0, 1]] 
        
    initializer = Initializer(landscape,population_size)
        
    init_pop = initializer.initialize()
        
    for ind in init_pop : 
        data.append([ind.RNA_seq, ind.RNA_structure, ind.mfe,ind.fitness])


    dataFrame = pandas.DataFrame(data)

    print (dataFrame)

    dataFrame.to_csv("../Logs/init_data/data1000_"+str(type_)+"_ "+str(init_depth)+".csv")
    
if __name__== "__main__" : 

     main()