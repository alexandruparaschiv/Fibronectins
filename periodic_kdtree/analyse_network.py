#!/usr/bin/env python3


''' This script analyses the clustering and the orientation of fibrils'''

import matplotlib.pyplot as plt
import math
import matplotlib.pyplot as plt
import seaborn as sns;sns.set()
import itertools
import pickle
import sys
import queue
from bfs import bfs

class Analysis:

    def __init__(self,filename):

        self.filename = filename
        self.right_ends = []
        self.left_ends = []

    def cluster_fibrils(self):

        with open(self.filename,'r') as infile:
            for item in infile:
                if len(item.split())>1  and 'Atoms' not in item:
                    if item[0] == '6':
                        at_type,x,y,z = item.split()
                        x = float(x);y=float(y);z=float(z)
                        self.right_ends.append((x,y,z))
                    if item[0] == '7':
                        at_type,x,y,z = item.split()
                        x = float(x);y=float(y);z=float(z)
                        self.left_ends.append((x,y,z))
    @staticmethod
    def get_distance(a1,a2):


        x1,y1,z1 = a1
        x2,y2,z2 = a2
        Lx = 100
        Ly = 300
        Lz = 300
        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        if dx>Lx/2:
            dx -= Lx
        elif dx<-Lx/2:
            dx += Lx
        if dy>Ly/2:
            dy -= Ly
        elif dy<-Ly/2:
            dy += Ly
        if dz>Lz/2:
            dz -= Lz
        elif dz<-Lz/2:
            dz += Lz
        return math.sqrt(dx**2+dy**2+dz**2)

    @staticmethod
    def get_vector(a1,a2):


        x1,y1,z1 = a1
        x2,y2,z2 = a2

        Lx = 100
        Ly = 300
        Lz = 300
        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        if dx>Lx/2:
            dx -= Lx
        elif dx<-Lx/2:
            dx += Lx
        if dy>Ly/2:
            dy -= Ly
        elif dy<-Ly/2:
            dy += Ly
        if dz>Lz/2:
            dz -= Lz
        elif dz<-Lz/2:
            dz += Lz
        return (dx,dy,dz)


    @staticmethod
    def dotproduct(a1,a2):
        x1,y1,z1 = a1
        x2,y2,z2 = a2
        dot_product = x1*x2+y1*y2+z1*z2
        norm1 = math.sqrt(x1**2+y1**2+z1**2)
        norm2 = math.sqrt(x2**2+y2**2+z2**2)
        dot_product /= (norm1*norm2)
        return dot_product

    @staticmethod
    def are_neighbouring_fibrils(a1,a2,b1,b2):

        cutoff = 50
        x1,y1,z1 = a1
        x2,y2,z2 = a2
        x3,y3,z3 = b1
        x4,y4,z4 = b2
        return Analysis.get_distance(a1,b1)<cutoff or Analysis.get_distance(a1,b2)<cutoff or Analysis.get_distance(a2,b1)<cutoff or Analysis.get_distance(a2,b2)<cutoff


if __name__ == "__main__":
    
    filename = sys.argv[1]
    analyser = Analysis(filename)
    analyser.cluster_fibrils()
    end_to_end_distances = []
    for i in range(len(analyser.right_ends)):
        a1 = (analyser.right_ends[i][0],analyser.right_ends[i][1],analyser.right_ends[i][2])
        a2 = (analyser.left_ends[i][0],analyser.left_ends[i][1],analyser.left_ends[i][2])
        dist = Analysis.get_distance(a1,a2)
        end_to_end_distances.append(dist)
    #sns.distplot(end_to_end_distances)
    #plt.show()

    orientations = []
    for i in range(len(analyser.right_ends)):
        u = (analyser.right_ends[i][0],analyser.right_ends[i][1],analyser.right_ends[i][2])
        v = (analyser.left_ends[i][0],analyser.left_ends[i][1],analyser.left_ends[i][2])
        vector = Analysis.get_vector(u,v)
        orientations.append(vector)

    dot_products = []
    for i in range(len(orientations)-1):
        for j in range(i+1,len(orientations)):
            dot_products.append(Analysis.dotproduct(orientations[i],orientations[j]))

    #print(dot_products)
    #sns.distplot(dot_products)
    #plt.show()

    # this part of the script calculates the orientations if fibrils are close
    neigh_orientations = []
    fibrils = []
    fibril_matrix = [ [0]*len(analyser.right_ends) for _ in range(len(analyser.right_ends))]
    for i in range(len(analyser.right_ends)-1):
        for j in range(i+1,len(analyser.right_ends)):
            if Analysis.are_neighbouring_fibrils(analyser.right_ends[i],analyser.left_ends[i],analyser.right_ends[j],analyser.left_ends[j]):
                u = Analysis.get_vector(analyser.right_ends[i],analyser.left_ends[i])
                v = Analysis.get_vector(analyser.right_ends[j],analyser.left_ends[j])
                #print(u,v)
                angle = math.degrees(math.acos(Analysis.dotproduct(u,v)))
                #print(angle)
                if angle < 20 and angle >-20:
                    fibril_matrix[i][j] = 1
                    fibril_matrix[j][i] = 1
    print(fibril_matrix)

    print(bfs(fibril_matrix))

    fibre_sizes = []
    fibrils = bfs(fibril_matrix)
    for item in fibrils:
        fibre_sizes.append(len(item))
    sns.kdeplot(fibre_sizes,shade=True)
    plt.show()
