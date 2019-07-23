import queue
from nematic import NematicCalculator
import pickle
import math

class InteractionMatrixBuilder:

    def __init__(self,radial_cutoff,df,data,periodic_ckdt_tree):

        self.cutoff = radial_cutoff
        no_monomers = df[-1:]['fibronectin_indices'].values[0]+1
        self.matrix = [[0]*no_monomers for _ in range(no_monomers)]

    def get_interactions(self,data,periodic_ckdt_tree,atom_to_monomer,df):

        for i in range(len(data)):
            distances, neighbour_indices = periodic_ckdt_tree.query(data[i], k=5)
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



    def get_alignment(fibril1,fibril2,fibronectins):

        interacting_atom_types = [2,3,8,9]
        subunits = fibronectins[(fibronectins['fibronectin_indices'].isin([fibril1,fibril2])) & (fibronectins['type'].isin(interacting_atom_types))]
        vectors = []
        fibro = pickle.load(open("fibronectin_system.pkl",'rb'))
        nematic = NematicCalculator()

        for i in range(subunits.shape[0]):
            if  not i % fibro.fibronectins.no_rods:

                dx = subunits.iloc[i+1].x - subunits.iloc[i].x
                dy = subunits.iloc[i+1].y - subunits.iloc[i].y
                dz = subunits.iloc[i+1].z - subunits.iloc[i].z
                dx,dy,dz = nematic.periodic_boundaries(dx,dy,dz)
                norm = math.sqrt(dx*dx+dy*dy+dz*dz)
                vec = (dx/norm,dy/norm,dz/norm)
                vectors.append(vec)
        return nematic.get_order_param(vectors)
