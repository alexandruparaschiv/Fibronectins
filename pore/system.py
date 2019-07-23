#!/storage/users/aparaschiv/anaconda3/bin/python


""" This creates the initial system of fibronectin monomers.
The system can be initialised in either a random geometry or on a lattice. """

from itertools import product


class System:

    side_length = 5
    def __init__(self):
        pass

    @staticmethod
    def make_lattice(fibro):

        N = System.side_length
        [fibro.make_monomer(i*fibro.x_spacing,j*fibro.y_spacing,k*fibro.z_spacing) for [i, j, k] in product(list(range(N)), repeat=3)]
