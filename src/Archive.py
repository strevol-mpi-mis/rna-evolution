
import numpy 


class Archive : 


    RANDOM_METHOD = "R"
    NOVEL_METHOD = "N"

    archiving = []

    def __init__(self, max_size, method ) :  

        self.max_size = max_size
        self.method = method 
        self.archiving = []


    def archive(self, population, novelties) :
        to_archive = []
        if self.method == self.RANDOM_METHOD : 
            to_archive = self.random_archiving(population)
        elif self.method == self.NOVEL_METHOD : 
            to_archive = self.novel_archiving(population,novelties) 

        for ind in to_archive : 
            self.archiving.append(ind) 



    def novel_archiving(self, population, novelties) : 
        
        to_archive = []
        sorted_novelties = sorted(novelties)#,reverse = True) 
        #sorted_novelties = numpy.array(sorted_novelties)
        #sorted_novelties = sorted_novelties[sorted_novelties>10]
        lamda_novelties = sorted_novelties[:self.max_size]

        for novelty in lamda_novelties : 
            to_archive.append(population[novelties.index(novelty)])


        return  to_archive



    def random_archiving(self, population) : 

        return numpy.random.choice(population, self.max_size, replace=False) 



