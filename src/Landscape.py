'''
@author: Nono Saha Cyrille Merleau 
@email: csaha@aims.edu.gh/ nonosaha@mis.mpg.de

This program is an implementation of a landscape function. 

The landscape function is given by: 


'''


#Importing the libraries 
import numpy as np 
import random
import math
import RNA


class Landscape(object) : 

        TREE_TYPE = "tree"
        BP_TYPE = "bp" 
        HAM_TYPE = "ham"

        #Landscape function
        ''' This function generates one landscape function for a given parameters: 
                INPUT
                =====
                        
                OUTPUT
                ======
                        
        '''
        def __init__ ( self, lamda, k, type_, target_structure) :
                self.lamda = lamda 
                self.k = k
                self.type = type_ 
                self.target_structure = target_structure 



        ############################
        #Compute the fitness of an RNA Structure
        def tree_edit_fitness(self, structure) : 
                ref_xstrc = RNA.expand_Full(self.target_structure)
                xstrc = RNA.expand_Full(structure)
                
                return 1./(1.+RNA.tree_edit_distance(RNA.make_tree(ref_xstrc), RNA.make_tree(xstrc)))

        
        #Compute the fitness of an RNA Structure
        def base_paire_fitness(self,structure1, structure2) : 
                return 1./(1.+RNA.bp_distance(structure1,structure2))

        #Compute the fitness of an RNA Structure
        def hamming_distance(self,ref_struc, structure) :
                distance = 0.
                for i in range(len(structure)) : 
                        if ref_struc[i] != structure[i] : 
                                distance +=1
                                
                return distance
        
        def fitness(self, structure) : 
                if self.type == self.BP_TYPE : 
                        return self.base_paire_fitness(self.target_structure,structure)
                elif self.type == self.TREE_TYPE : 
                        return self.tree_edit_fitness(structure)
                else: 
                        return 1/(1. + self.hamming_distance(self.target_structure, structure))
        

        def novelty(self, structure, population) : 
                list_novelty_metrics = []

                for ind in population : 

                        #list_novelty_metrics.append(self.base_paire_fitness(structure,ind.RNA_structure))
                        list_novelty_metrics.append( 1/(1.+self.hamming_distance(structure,ind.RNA_structure)))

                
                return np.mean(list_novelty_metrics)
