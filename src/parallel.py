'''
    @author: Nono Saha Cyrille Merleau 
    @email: nonosaha@mis.mpg.de 

    This is a small distributed program that the goal is to run an RNA evolution in parallel on a simple cluster by variating the mutation probability. 

    Enjoy!!!!

'''

#Standard libraries
import sys
import os

#External libraries 
import numpy as np 
from Individual import Individual
import RNAEvolution
import pp 
import RNA as rna 
import pandas
import matplotlib.pyplot as plt



def main(): 
    
   # tuple of all parallel python servers to connect with
    ppservers = ()
    #ppservers = ("10.0.0.1",)

    if len(sys.argv) > 1:
        ncpus = int(sys.argv[1])
        # Creates jobserver with ncpus workers
        job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers=ppservers)

    print "Starting pp with", job_server.get_ncpus(), "workers"



    population_size = 500
    number_of_generation = 100
    init_deph = 80
    #mut_probs = [1./(init_deph-i) for i in range(17,80,18)]
    mut_probs = [0.01,0.02,0.03,0.04]
    print mut_probs
    nucluotides = ["A", "G", "U", "C"]
    arn = np.random.choice(nucluotides,50)
    rna_str = ''.join(arn)
    print "rna_String", rna_str
     # compute minimum free energy (MFE) and corresponding  referent structure
    (ref_ss, mfe) = rna.fold(rna_str)
    ref_ind = Individual(rna_str,ref_ss,0)
   
    jobs = [(mut_prob, job_server.submit(RNAEvolution.evolution,(ref_ind, population_size, number_of_generation, mut_prob,), (RNAEvolution.initialize,RNAEvolution.select,RNAEvolution.fitness,RNAEvolution.mutateAll,RNAEvolution.mutateOne,RNAEvolution.getFitnesses,Individual),
     ("numpy","Individual","RNAEvolution","RNA","random"))) for mut_prob in mut_probs]
   
    
    for mut_prob, job in jobs:
        print "Job ,", mut_prob
        fig, ax = plt.subplots()
        os.mkdir("Job"+str(mut_prob))
        means = []
        i = 1; 
        for gen in job() : 
            data = pandas.DataFrame(sorted(gen))
            data.to_csv("Job"+str(mut_prob)+"/gen"+str(i)+".csv",sep=",")
            means.append(np.mean(gen))
            i = i +1 
        plt.plot(means)
    
    plt.savefig("History.mean.eps")
    plt.show()
           
        #dataFrame = pandas.DataFrame(job())
        #dataFrame.plot.kde(ax=ax,legend=False, title='Histogram:F'+str(mut_prob))
        #dataFrame.to_csv("Histogram:F"+str(mut_prob)+".csv",sep=",")
        #plt.savefig("Histogram:F"+str(mut_prob)+".eps")

    



if __name__== "__main__": 
    
    main()