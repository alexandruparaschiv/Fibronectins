#!/usr/bin/env python3

''' This class is responsible for making the collective variables file'''

import random
import pickle
from fibronectin import Fibronectins
import argparse



parser = argparse.ArgumentParser()
parser.add_argument("percent_pulled", help="adjusts the fraction of fibronectins extended",
                    type=float,default=1.0)
args = parser.parse_args()
pulled = args.percent_pulled


class Colvars:

    def __init__(self):

        self.fibronectins = Fibronectins()

        self.colvar_filename = "config.colvars"
        self.rest_length = 250
        self.fraction = pulled 

        self.protein_ends = self.group_ends(list(filter(lambda x : x[2] in (6,7), self.fibronectins.atoms)))
        self.adjust_fraction_pulled()
        self.write_file()

        print("The collective variables file was created successfully!")

    def __str__(self):

        output = "Pulling " + str(int(self.fraction*100))+"% of the fibronectin molecules."
        
        return output

    def make_colvar(self,index):
        
        colvar = ["colvar {\n","\tname colvar"+str(index),"\n\twidth 0.2\n","\tdistance {\n","\t\tgroup1 {\n"]
        colvar += ["\t\t\tatomNumbers "+str(self.protein_ends[index][0])+" \n","\t\t}\n","\t\tgroup2 {\n"]
        colvar += ["\t\t\tatomNumbers "+str(self.protein_ends[index][1])+" \n","\t\t}\n"]
        colvar += ["\t}\n","}"]
        return ''.join(colvar)

    def group_colvars(self):

        no_colvars = len(self.protein_ends)
        header = ["harmonic {\n","\tname fibronectin_unfolding\n\tcolvars"]
        header += [" colvar"+str(i)+" " for i in range(no_colvars)]
        header += ["\n\tcenters"]
        header += [" "+ str(self.rest_length) for _ in range(no_colvars)]
        header += ["\n\tforceConstant 0.01\n}"]
        return ''.join(header)

    def group_ends(self,protein_ends):

        bonded_atoms = []
        for index,item in enumerate(protein_ends):
            if item[2] == 7:
                bonded_atoms.append((protein_ends[index-1][0],protein_ends[index][0]))
        return bonded_atoms

    def adjust_fraction_pulled(self): 

        no_protein_ends_pulled = self.fraction * len(self.protein_ends)
        while len(self.protein_ends) > no_protein_ends_pulled:
            elem = random.choice(self.protein_ends)
            self.protein_ends.remove(elem)


    def write_file(self):

        with open(self.colvar_filename,'w') as outfile:
            for index,item in enumerate(self.protein_ends):
                outfile.write("{}\n\n".format(self.make_colvar(index)))

        with open(self.colvar_filename,'a') as outfile:
            outfile.write("{}".format(self.group_colvars()))



if __name__ == "__main__":

    colvar = Colvars()
    print(colvar)
    pickle.dump(colvar, open( "fibronectin_system.pkl", "wb" ))


