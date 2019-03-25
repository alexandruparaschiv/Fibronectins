#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from utilities import get_coords, get_distance, dot_product
import os 


filename = "fibro.xyz"
no_atoms = 176
line_counter = 0
dummy_atom = [] 
smd_atom  = []
pulling_direction = (0,1,0)
v = 0.1
tau0 = 0.008
k = 10



settings_file = "in.fibronectins"

	
with open(filename,'r') as infile:
	for item in infile:
		line_counter += 1 
		if line_counter % (no_atoms+2) == 0:
			dummy_atom.append(get_coords(item.split()))			
		elif line_counter % (no_atoms+2) == no_atoms+1:
			smd_atom.append(get_coords(item.split()))
		


displacement = []
x0 = smd_atom[0][0]
y0 = smd_atom[0][1]
z0 = smd_atom[0][2]

for item in smd_atom:

	displacement.append(dot_product((item[0]-x0,item[1]-y0,item[2]-z0),pulling_direction))


vt = np.linspace(0,v*len(displacement)*100,len(displacement))
vt *= tau0
displacement = np.asarray(displacement)
force = k*(vt - displacement)


np.savetxt('force_displacement.txt', np.column_stack([force,displacement]))


