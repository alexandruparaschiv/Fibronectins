#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pickle
from colvars import Colvars
from tqdm import tqdm
import queue
import math
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()
from nematic import NematicCalculator
from periodic_kdtree.periodic_kdtree import PeriodicCKDTree



class MovieBuilder:
    def __init__(self,types,x,y,z):

        fibro = pickle.load(open("fibronectin_system.pkl",'rb'))

        xbound = 2*fibro.fibronectins.side_lengthx
        ybound = 2*fibro.fibronectins.side_lengthy
        zbound = 2*fibro.fibronectins.side_lengthz
        bounds = np.array([xbound, ybound, zbound])

        interacting_atom_types = [2,3,8,9]


        #df = pd.read_csv('last_frame.xyz',skiprows=2,delimiter=' ',names=['type','x','y','z'])
        names = ['type','x', 'y', 'z']
        arrays = [types,x,y,z]
        df = pd.DataFrame(dict(zip(names, arrays)), columns=names)
        coloured_df = df


        fibronectin_ends = sorted(list(set(fibro.fibronectins.monomer_chunks)))
        fibronectin_length = fibronectin_ends[1] - fibronectin_ends[0]
        fibronectin_indices = np.array([i//(fibronectin_length) for i in range(df.shape[0])])
        fibronectin_indices = pd.DataFrame({'fibronectin_indices': fibronectin_indices})
        df = df.join(fibronectin_indices)
        nematic = NematicCalculator()
        fibronectins = df
        atomic_indices = [ i for i in range(df.shape[0]) ]
        atomic_indices = pd.DataFrame({'atomic_indices': atomic_indices})
        df = df.join(atomic_indices)
        atom_to_monomer = { df['atomic_indices'][i]:df['fibronectin_indices'][i] for i in range(df.shape[0])}
        print(df)
        df = df[df['type'].isin(interacting_atom_types)]
        print(df)
        coordinates = np.array([df['x'].values,df['y'].values,df['z'].values])
        data = coordinates.T
        T = PeriodicCKDTree(bounds,data)


        class InteractionMatrixBuilder:
    

            def __init__(self):
        
                self.cutoff = 1.3*1.12
                no_monomers = df[-1:]['fibronectin_indices'].values[0]+1
                self.matrix = [[0]*no_monomers for _ in range(no_monomers)]    
        
            def get_interactions(self):
        
                for i in range(len(data)):
                    distances, neighbour_indices = T.query(data[i], k=5)
                    for j in range(len(neighbour_indices)):
                        if distances[j] < self.cutoff and atom_to_monomer[int(df.iloc[i].atomic_indices)] != atom_to_monomer[int(df.iloc[neighbour_indices[j]].atomic_indices)]:
                            u = atom_to_monomer[int(df.iloc[i].atomic_indices)]
                            v = atom_to_monomer[int(df.iloc[neighbour_indices[j]].atomic_indices)]
                            self.matrix[u][v],self.matrix[v][u] = 1,1  
                return self.matrix
        
            @staticmethod
            def clustering(interaction_matrix):
        
                visited = [ False for _ in range(len(interaction_matrix)) ]
                fibrils = []
                q = queue.Queue()
        
                for i in range(len(interaction_matrix)):
                    if visited[i]:
                        continue
                    else:
                        fibrils.append([i])
                        q.put(i)
                    while not q.empty():
                        to_visit = q.get()
                        if not visited[to_visit]:
                            visited[to_visit] = True
                            for j in range(len(interaction_matrix)):
                                if interaction_matrix[to_visit][j]:
                                    if not visited[j]:
                                        q.put(j)
                                        fibrils[-1].append(j)
                                        visited[j] = True
                return fibrils
    
            @staticmethod
            def get_all_order_params(fibrils):
        
                order_params = []
                for item in fibrils:
                    subunits = fibronectins[(fibronectins['fibronectin_indices'].isin(item)) & ((fibronectins['type']==2)|(fibronectins['type']==3))]
                    vectors = []
                    for i in range(subunits.shape[0]):
                        if  not i % fibro.fibronectins.no_rods:
                            dx = subunits.iloc[i+1].x - subunits.iloc[i].x 
                            dy = subunits.iloc[i+1].y - subunits.iloc[i].y 
                            dz = subunits.iloc[i+1].z - subunits.iloc[i].z 
                            dx,dy,dz =nematic.PBC(dx,dy,dz)
                            norm = math.sqrt(dx*dx+dy*dy+dz*dz)
                            vec = (dx/norm,dy/norm,dz/norm)
                            vectors.append(vec)
                    order_params.append(nematic.get_order_param(vectors))
                return order_params
    
            def get_alignment(fibril1,fibril2):

                subunits = fibronectins[(fibronectins['fibronectin_indices'].isin([fibril1,fibril2])) & ((fibronectins['type']==2)|(fibronectins['type']==3))]
                vectors = []
            
                for i in range(subunits.shape[0]):
                    if  not i % fibro.fibronectins.no_rods:
                    
                        dx = subunits.iloc[i+1].x - subunits.iloc[i].x 
                        dy = subunits.iloc[i+1].y - subunits.iloc[i].y 
                        dz = subunits.iloc[i+1].z - subunits.iloc[i].z 
                        dx,dy,dz = nematic.PBC(dx,dy,dz)
                        norm = math.sqrt(dx*dx+dy*dy+dz*dz)
                        vec = (dx/norm,dy/norm,dz/norm)
                        vectors.append(vec)
                return nematic.get_order_param(vectors)

        class Fibrils:
    
    
            def __init__(self):
        
                matrix = InteractionMatrixBuilder()
                interaction_matrix = matrix.get_interactions()
                self.clusters = InteractionMatrixBuilder.clustering(interaction_matrix)
                self.aggregates = list(filter(lambda x: len(x)>1, self.clusters))
                self.alignment_matrices = []
                self.network = []
                self.nematic_cutoff = 0.4
        
                self.generate_fibril_network()
        
            def generate_fibril_network(self):
        
                for aggregate in self.aggregates:       
                    alignment_matrix = [[0]*len(aggregate) for _ in range(len(aggregate))]    
                    for i in range(len(aggregate)):
                        for j in range(len(aggregate)):
                            alignment = InteractionMatrixBuilder.get_alignment(aggregate[i],aggregate[j])
                            if alignment > self.nematic_cutoff:
                                alignment_matrix[i][j],alignment_matrix[j][i] = 1,1
                    alignment_matrix = InteractionMatrixBuilder.clustering(alignment_matrix)
            
                    fibrils = []
            
                    for item in alignment_matrix:
                        fibril = []
                        for fibril_indices in item:
                            fibril.append(aggregate[fibril_indices])
                        fibrils.append(fibril[:])
                    self.network.extend(fibrils)

        fibrils = Fibrils()
        fibrils.generate_fibril_network()
        new_types = np.ones_like(coloured_df.type.values)
        for cluster in fibrils.network:
            for monomer in cluster:
                new_types[monomer*fibronectin_length:(monomer+1)*fibronectin_length] = len(cluster)
        coloured_df.type = new_types
        with open("coloured.xyz",'a') as infile:
            infile.write("50000\n")
            infile.write("Atoms. Timestep: 0\n")
            for i in range(coloured_df.type.values.shape[0]):
                infile.write("{} {} {} {}\n".format(coloured_df.type.values[i],coloured_df.x.values[i],coloured_df.y.values[i],coloured_df.z.values[i]))



		


