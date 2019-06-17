
'''
    @author: Nono Saha Cyrille Merleau  
    @email: nonosaha@mis.mpg.de

'''
import os
import pandas


class Logger(object) : 

    ROOT_LOG_PATH = "../Logs/BenshMark/Defect/"

    def __init__(self, logfolder, sublogfolder): 

        self.logfolder = logfolder
        self.root_path = self.ROOT_LOG_PATH+str(logfolder)+"/"+sublogfolder
        try:
            os.makedirs(self.root_path)
        except OSError :
            print (" Can not initialize the log folder ")

    def save_population(self, population,gen) : 
        data = []
        for ind in population : 
            data.append([ind.RNA_seq, ind.RNA_structure,ind.mfe, ind.fitness])
    
        dataFrame = pandas.DataFrame(data)
        dataFrame.to_csv(self.root_path+"/gen"+str(gen)+".csv")

    def bt_save_population(self, prev_pop, population,gen) : 
        data = []
        prev_data = []
        for i in range(len(population)) : 
            data.append([population[i].RNA_seq, population[i].RNA_structure,population[i].mfe, population[i].fitness])
            prev_data.append([prev_pop[i].RNA_seq, prev_pop[i].RNA_structure,prev_pop[i].mfe, prev_pop[i].fitness])

        
        dataFrame = pandas.DataFrame(data)
        prev_dataFrame = pandas.DataFrame(prev_data)
        prev_dataFrame.to_csv(self.root_path+"/prev_gen"+str(gen)+".csv")
        dataFrame.to_csv(self.root_path+"/gen"+str(gen)+".csv")
