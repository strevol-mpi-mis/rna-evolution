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

import subprocess 


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
                        metric = self.hamming_distance(structure,ind.RNA_structure)
                        list_novelty_metrics.append(1/(1.+metric))
                
                return np.mean(list_novelty_metrics)

        def min_ens_distance(self, sequence, min_energy) : 
                rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                rnasubopt_out, rnasubopt_err = rnasubopt.communicate(sequence) 

                cut_pipe = subprocess.Popen(['cut', '-f', '1', '-d', ' '], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                cut_out, stderr = cut_pipe.communicate(rnasubopt_out)
                
                list_subopt =  cut_out.split()[1:] 

                return max([self.fitness(subopt) for subopt in list_subopt])



        def ens_defect(self, sequence) : 
                p = subprocess.Popen("defect", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                cmd_input = sequence+"\n"+self.landscape.target_structure+"\n"
                
                defect, error = p.communicate(cmd_input)
                defect = defect.split("\n")
                return 1/float(defect[-3])
        
        def my_ens_def(self, sequence, min_energy) : 
                rnasubopt = subprocess.Popen(args=['RNAsubopt', '-e', str(min_energy)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                rnasubopt_out, rnasubopt_err = rnasubopt.communicate(sequence) 
                result = np.array([[s.split()[0], s.split()[1]] for s in rnasubopt_out.split("\n")[1:-1]])
                d_s = np.array([self.fitness(sigma) for sigma in result[:,0]])
                kt = 0.612 
                z = sum(np.exp(-np.array(result[:,1], dtype=float)/kt))
                p = np.exp(-np.array(result[:,1], dtype=float)/kt)/z

                return sum(p*d_s)