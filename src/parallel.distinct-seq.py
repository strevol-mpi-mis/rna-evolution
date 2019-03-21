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

def get_distinct_structure(interval, folder_num) : 
    lamdas = [1]
    result = []
    for lamda in lamdas : 
        archive_strc = numpy.array([], dtype="S30")
        data = []
        for i in range(interval[0],interval[1]) : 
            df = pandas.read_csv("../Logs/NoveltyExp/News/old-select/old-select/"+str(folder_num)+"/"+str(lamda)+"/gen"+str(i)+".csv", sep=",")
            
            population  = numpy.array( [ind[2] for ind in df.values])
            #print population
            population = numpy.insert(population,len(population),archive_strc)
            strc = list(set(population))
            archive_strc = numpy.insert(archive_strc,len(archive_strc),numpy.array(strc))
            nb_strc = len(strc)
    print 'Job numer ',interval," = ",len(strc)    
    return strc




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
    
    folders = range(1,  folder_number+101, 100)
    intervales = []
    for i in range(len(folders)-1) : 
        intervales.append((folders[i], folders[i+1]))
    print intervales
    #Parallel evolution for every lamda value
    print "Start running jog"
    jobs = [(interval, job_server.submit(get_distinct_structure, (interval,1), modules=("numpy","pandas","os"))) for interval in intervales]
    
    archive_strc = numpy.array([], dtype="S30")

    for folder, job in jobs:
        res = job()
        archive_strc = numpy.insert(archive_strc,len(archive_strc),numpy.array(res))
    print "The total number of distinct structures is ", len(set(archive_strc))
    
          
   
    

    
if __name__== "__main__" : 

    main()