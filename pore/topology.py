#!/usr/bin/env python3

import math
import os
import time

class Topology:

    def __init__(self):
        pass

    @staticmethod
    def write_topology_file(fibro):

        N = 150
        initial_atoms = str(fibro.atomic_index+N*N)
        header = "LAMMPS Description\n\n" + "\t" + str(fibro.atomic_index+N*N) + " atoms\n" + "\t" + \
        str(len(fibro.bonds)) + " bonds\n \t0 angles\n \t0 dihedrals\n \t0 impropers\n\n \t10 atom types\n \t1 bond types" + \
        "\n \t0 angle types\n \t0 dihedral types\n \t0 improper types\n\n \t" + str(-fibro.side_lengthx) + \
        " " + str(fibro.side_lengthx) + " xlo xhi\n \t" + str(-fibro.side_lengthy) + \
        " " + str(fibro.side_lengthy) + " ylo yhi\n \t" + str(-fibro.side_lengthz) + \
        " " + str(fibro.side_lengthz) + " zlo zhi\n \t"


        header += "\nMasses\n\n \t1 "+str(fibro.atom_mass)+"\n\t2 "+str(fibro.atom_mass)+\
        "\n\t3 "+str(fibro.atom_mass)+"\n\t4 "+str(fibro.atom_mass)+"\n\t5 "+str(fibro.atom_mass)+\
        "\n\t6 "+str(fibro.atom_mass)+"\n\t7 "+str(fibro.atom_mass)+"\n\t8 "+str(fibro.atom_mass)+\
        "\n\t9 "+str(fibro.atom_mass)+"\n"+"\t10 "+str(fibro.atom_mass)+"\n\n"

        header += "Atoms\n\n"
        z_offset = 100
        for i in range(len(fibro.atoms)):
            header += "\t"+" "+str(fibro.atoms[i][0])+" "+str(fibro.atoms[i][1])+" "+str(fibro.atoms[i][2])+" 0 "+\
            str(fibro.atoms[i][3])+" "+str(fibro.atoms[i][4])+" "+str(fibro.atoms[i][5]-z_offset)+" 0 0 0\n"

        initial_atomic = len(fibro.atoms) + 1
        initial_molecular = initial_atomic
        curr_index = -1
        xmid = ymid = 0.5*N/2
        xmid = ymid = 0
        r0 = 10
        for i in range(-N,N):
            for j in range(-N,N):
                for k in range(5):
                    #curr_index += 1
                    x = 1.0*i; y = 1.0*j; z = 50 +1*k
                    r = math.sqrt((xmid-x)*(xmid-x)+(ymid-y)*(ymid-y))
                    if r>r0:
                        curr_index += 1
                        header += "\t"+" "+str(initial_atomic+curr_index)+" "+str(initial_molecular+curr_index)+" 10 0 "+\
                        str(x)+" "+str(y)+" "+str(z)+" 0 0 0\n"
        #print(str(initial_atomic+curr_index)+" atoms ")
        final_atoms = str(initial_atomic+curr_index)

        header += "Bonds\n\n"

        for i in range(len(fibro.bonds)):
            header += "\t"+" "+str(i+1)+" 1 "+ str(fibro.bonds[i][0])+" "+\
                        str(fibro.bonds[i][1])+"\n"

        with open(fibro.topology_file,'w') as f:
            f.write("{}".format(header))

        with open(fibro.topology_file,'r') as infile:
            with open("safe_copy",'w') as outfile:
                for item in infile:
                    if " atoms" in item:
                        outfile.write("{} atoms\n".format(final_atoms))
                    else:
                        outfile.write("{}".format(item))
        os.system("mv safe_copy topology.dat")
