#!/usr/bin/python




def get_coords(alist):

        x = float(alist[1])
        y = float(alist[2])
        z = float(alist[3])
        return (x,y,z)

def get_distance( coord1, coord2):
        # coord1 and coord2 are tuples containing 
        # the coordinates 
        dx = coord1[0] - coord2[0]
        dy = coord1[1] - coord2[1]
        dz = coord1[2] - coord2[2]
        r = math.sqrt(dx*dx+dy*dy+dz*dz)
        r = 1
        return (dx/r,dy/r,dz/r)

def dot_product(vec1,vec2):

        dot = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]
        return dot




