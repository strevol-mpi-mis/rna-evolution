
'''
    @author: Nono Saha Cyrille Merleau  
    @email: nonosaha@mis.mpg.de

'''
import os
import pandas


class Logger(object) : 

    ROOT_LOG_PATH = "../Logs/NoveltyExp/News/"

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
            data.append([ind.RNA_seq, ind.RNA_structure, ind.fitness])
    
        dataFrame = pandas.DataFrame(data)
        dataFrame.to_csv(self.root_path+"/gen"+str(gen)+".csv")