#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt

class Parser(object):

	def __init__ (self,filename):

		self.filename = filename
		self._count_atoms()

	def _count_atoms(self):

		with open(self.filename,'r') as infile:
			self._atoms_number = int(infile.readline())

	def go_through_frames(self):

		frame_index = 0
		atom_index = 0
		types = []
		coordinates = []

		with open(self.filename,'r') as infile:
			for item in infile:
				if len(item.split()) and types:

					frame_index +=1;
					types = np.asarray(types)
					coordinates = np.asarray(coordinates)
					types = [ ]
					coordinates = [ ]

				if len(item.split()) > 3:

					atom_type,x,y,z=item.split();
					atom_type = str(atom_type)
					x,y,z = float(x),float(y),float(z)
					types.append(atom_type)
					coordinates.append((x,y,z))
if __name__ == "__main__" :

	parser = Parser('fibro.xyz')
	parser._count_atoms()
	parser.go_through_frames()
