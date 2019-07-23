#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pickle
import os
from colvars import Colvars
import sys
from fibrils import FibrilsNetwork
from periodic_kdtree.periodic_kdtree import PeriodicCKDTree
import itertools

try:
    from nematic import NematicCalculator
except ImportError:
    print(sys.exc_info())


class FibrilAnalyser:

    def __init__(self):

        self.fibro = pickle.load(open("fibronectin_system.pkl", 'rb'))
        xbound = 2*self.fibro.fibronectins.side_lengthx
        ybound = 2*self.fibro.fibronectins.side_lengthy
        zbound = 2*self.fibro.fibronectins.side_lengthz
        self.bounds = np.array([xbound, ybound, zbound])
        self.interacting_atom_types = [2, 3, 8, 9]

        self.column_names = ['type', 'x', 'y', 'z']
        self.fibronectin_ends = sorted(list(set(self.fibro.fibronectins.monomer_chunks)))
        self.fibronectin_length = self.fibronectin_ends[1] - self.fibronectin_ends[0]
        self.outputfile = "average_fibril_size.txt"

        print(os.getcwd())
        with open(self.outputfile,'w') as infile:
            pass

    def analyse_fibril_distribution(self,nematic_cutoff, radial_cutoff,
        types, x, y, z):

        arrays = [types,x,y,z]
        df = pd.DataFrame(dict(zip(self.column_names, arrays)), columns=self.column_names)

        """ Add columns corresponding to mapping between the atomic index and the monomer index. """
        fibronectin_indices = np.array([i//(self.fibronectin_length) for i in range(df.shape[0])])
        fibronectin_indices = pd.DataFrame({'fibronectin_indices': fibronectin_indices})
        no_monomers = df.shape[0]//self.fibronectin_length
        df = df.join(fibronectin_indices)

        fibronectins = df
        atomic_indices = [i for i in range(df.shape[0])]
        atomic_indices = pd.DataFrame({'atomic_indices': atomic_indices})
        df = df.join(atomic_indices)

        atom_to_monomer = { df['atomic_indices'][i]:df['fibronectin_indices'][i]
        for i in range(df.shape[0]) }

        """ Filters out only the relevant patches for the clustering algorithm. """
        df = df[df['type'].isin(self.interacting_atom_types)]

        """ Builds the periodic Cython-KDTree """
        data = np.array([df['x'].values, df['y'].values, df['z'].values]).T
        periodic_ckdt_tree = PeriodicCKDTree(self.bounds,data)

        """ Builds the fibronectin fibers network """
        fibrils = FibrilsNetwork(radial_cutoff=radial_cutoff, nematic_cutoff=nematic_cutoff, df=df,
            data=data, periodic_ckdt_tree=periodic_ckdt_tree, atom_to_monomer=atom_to_monomer, fibronectins=fibronectins)
        fibrils.generate_fibril_network(fibronectins)
 
        """ Cleans-up the fibrilar network."""
        k = sorted(fibrils.network)
        network = list(k for k,_ in itertools.groupby(k)) # eliminates duplicates in lists of lists
        isolated_fibrils = [len(fibril) for fibril in network]

        """ Outputs the analysis result (average fibril size, average nematic order parameter). """
        monomers = no_monomers - sum(isolated_fibrils)
        counts = np.zeros(no_monomers, dtype=np.int32)
        counts[0] += monomers
        for item in isolated_fibrils:
            counts[item-1] += 1

        """ Writes the average fibril size to a file. """    
        with open(self.outputfile,'a') as infile:
            infile.write("{}\n".format(1/np.mean(counts)))
