#!/usr/bin/env python3

import pickle
from colvars import Colvars
import pandas as pd

fibro = pickle.load(open("fibronectin_system.pkl",'rb'))
print(fibro.fibronectins.atoms)


interacting_atom_types  = [2,3,8,9]


atoms = fibro.fibronectins.atoms
patches = list(filter(lambda x: x[2] in interacting_atom_types, atoms))
print(patches)


df = pd.read_csv('fibro.xyz',skiprows=1)
print(df.size)

