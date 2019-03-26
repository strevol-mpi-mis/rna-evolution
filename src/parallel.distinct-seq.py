'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA as rna 
import random 
from Individual import Individual
import matplotlib.pyplot as plt
import pandas
import os
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import sys
import pp


def hamming(genotype1, genotype2) :
    distance = 0
    b1 = genotype1
    b2 = genotype2
    #print (len(b1), len(b2))
    #print (b1, b2)
    for i in range(len(b1)) : 
        if b1[i] != b2[i] : 
            distance +=1
                
    return distance

def nb_clusters(vector) : 
    cluster = [vector[0]]
    vector.remove(vector[0])
    i = 1
    while i <= len(vector) + 1  : 

        if hamming(cluster[0],vector[0]) == 0 : 
            cluster.append(vector[0])
            vector.remove(vector[0])


        else : 
            i +=1
            cluster = [vector[0]]
            vector.remove(vector[0])
    
    return i

def clusters(vector) : 
    vector = list(vector)
    result = []
    cluster = [vector[0]]
    vector.remove(vector[0])
    while len(vector) > 0  : 
        pointer = np.copy(vector)
        for i in range(len(pointer)) : 
            if hamming(cluster[0],pointer[i]) == 0 : 
                cluster.append(pointer[i])
                vector.remove(pointer[i])
        result.append(cluster[0])
        if len(vector)!=0 : 
            cluster = [vector[0]]
        
    return result
"""
def get_distinct_structure(interval, folder_num) : 
    lamdas = [1]
    result = []
    for lamda in lamdas : 
        archive_strc = numpy.array([], dtype="S30")
        data = []
        for i in range(interval[0],interval[1]) : 
            for lamda in  lamdas : 
                maxs = []
                files = os.listdir("../Logs/NoveltyExp/News/Mut80-9-1500/"+str(n)+"S/"+str(lamda)+"/")
                nb_gen.append(len(files)- 1) 
        histoData.append(nb_gen)

    histoData = np.array(histoData)
    print 'Job numer ',interval," = ",len(strc)    
    return strc
"""
def get_min_generation(interval) : 
    lamdas = [1]
    result = []
    methods = ["N", "F", "S", "FREQ"]
    
    for i in range(interval[0],interval[1]) : 
        min_gen  = []
        for method in  methods : 
            try : 
                files = os.listdir("../Logs/MyTest2/"+str(i)+"/"+str(method)+"/")
                min_gen.append(len(files)- 1) 
            except : 
                print "Folder ", i, "missed"
                continue
            
        if len(min_gen) == 4 : 
            result.append(min_gen)
    result = numpy.array(result)
    print 'Job numer ',interval," = ",  result.shape
    return numpy.array(result)


def distinct(folder_num) : 
    lamdas = [0,0.5,1]
    result = []
    for lamda in lamdas : 
        archive_strc = numpy.array([], dtype="S80")
        archive_seq = numpy.array([], dtype="S80")
        data = []
        for i in range(1,1001) : 
            df = pandas.read_csv("../Logs/NoveltyExp/Cross_Mut80-11/"+str(folder_num)+"S/"+str(lamda)+"/gen"+str(i)+".csv", sep=",")
            population  = numpy.array([ind[1] for ind in df.values])
            population = numpy.insert(population,len(population),archive_seq)
            #seq = clusters(population)
            seq = list(set(population))
            print i
            archive_seq = numpy.insert(archive_seq,len(archive_seq),numpy.array(seq))
            nb_seq = len(seq)
            population  = numpy.array( [ind[2] for ind in df.values])
            population = numpy.insert(population,len(population),archive_strc)
            #strc = clusters(population)
            strc = list(set(population))
            archive_strc = numpy.insert(archive_strc,len(archive_strc),numpy.array(strc))
               
            nb_strc = len(strc)
            data.append([nb_seq,nb_strc])
            
        print (len(data), data)  
        data = numpy.array(data)
        df = pandas.DataFrame(data) 
        df.to_csv("../Logs/NoveltyExp/Cross_Mut80-11/"+str(folder_num)+"S/"+str(lamda)+".csv", sep=",")



#Main function 
def main() : 
    
    #tuple of all parallel python servers to connect with
    ppservers = ()
    #ppservers = ("10.0.0.1",)
    
    if len(sys.argv) > 1:
        ncpus = int(sys.argv[1])
        folder_number = int(sys.argv[2])

        # Creates jobserver with ncpus workers
        job_server = pp.Server(ncpus, ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers=ppservers)

    print "Starting pp with", job_server.get_ncpus(), "workers"
    print "Folder Number == ", folder_number
    
    folders = range(1,  folder_number+25, 24)
    intervales = []
    for i in range(len(folders)-1) : 
        intervales.append((folders[i], folders[i+1]))
    print intervales
    #Parallel evolution for every lamda value
    print "Start running jog"
    jobs = [(interval, job_server.submit(get_min_generation, (interval,), modules=("numpy","pandas","os"))) for interval in intervales]
    
    #archive_strc = numpy.zeros((24,4))
    archive_strc = [] 
    for folder, job in jobs:
        res = job()
        archive_strc.append(res) 
    archive_strc = numpy.concatenate(archive_strc)

    print archive_strc.shape
    print "                                      ", ["|| N ||", "|| F ||", "|| S ||", "|| FREQ ||"]
    print "The average number of generations is ", numpy.median(archive_strc, axis=0)
    
          
   
    

    
if __name__== "__main__" : 

    main()
