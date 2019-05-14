#!/usr/bin/env python3

import math
import pickle
from colvars import Colvars

class NematicCalculator:

    def __init__(self):

        fibro = pickle.load(open("fibronectin_system.pkl",'rb'))
        self.xbound = xbound = 2*fibro.fibronectins.side_lengthx
        self.ybound = ybound = 2*fibro.fibronectins.side_lengthy
        self.zbound = zbound = 2*fibro.fibronectins.side_lengthz
        pass

    def PBC(self,dx,dy,dz):

        """applies periodic boundary conditions"""
        if dx > self.xbound/2:
            dx -= self.xbound
        elif dx <- self.xbound/2:
            dx += self.xbound
        if dy > self.ybound/2:
            dy -= self.ybound
        elif dy <- self.ybound/2:
            dy += self.ybound
        if dz > self.zbound/2:
            dz -= self.zbound
        elif dz <- self.zbound/2:
            dz += self.zbound
        return (dx,dy,dz)

    @staticmethod
    def get_cosine(u,v):

        dot_product = u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
        mod_u = math.sqrt(u[0]**2+u[1]**2+u[2]**2)
        mod_v = math.sqrt(v[0]**2+v[1]**2+v[2]**2)
        cosine = dot_product/(mod_u*mod_v)
        if cosine > 1:
            cosine = 1
        elif cosine < -1:
            cosine = -1
        return cosine

    def get_order_param(self,vectors):

        nematic_parameter = 0
        counter = 0
        for i in range(len(vectors)-1):
            for j in range(i+1,len(vectors)):
                counter += 1
                nematic_parameter += self.get_cosine(vectors[i],vectors[j])**2

        nematic_parameter /= counter
        return 1.5*nematic_parameter-0.5
