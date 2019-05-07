#!/usr/bin/env python3

class Topology:

    def __init__(self):
        pass

    @staticmethod
    def write_topology_file(fibro):

        header = "LAMMPS Description\n\n" + "\t" + str(fibro.atomic_index) + " atoms\n" + "\t" + \
        str(len(fibro.bonds)) + " bonds\n \t0 angles\n \t0 dihedrals\n \t0 impropers\n\n \t9 atom types\n \t1 bond types" + \
        "\n \t0 angle types\n \t0 dihedral types\n \t0 improper types\n\n \t" + str(-fibro.side_lengthx) + \
        " " + str(fibro.side_lengthx) + " xlo xhi\n \t" + str(-fibro.side_lengthy) + \
        " " + str(fibro.side_lengthy) + " ylo yhi\n \t" + str(-fibro.side_lengthz) + \
        " " + str(fibro.side_lengthz) + " zlo zhi\n \t"


        header += "\nMasses\n\n \t1 "+str(fibro.atom_mass)+"\n\t2 "+str(fibro.atom_mass)+\
        "\n\t3 "+str(fibro.atom_mass)+"\n\t4 "+str(fibro.atom_mass)+"\n\t5 "+str(fibro.atom_mass)+\
        "\n\t6 "+str(fibro.atom_mass)+"\n\t7 "+str(fibro.atom_mass)+"\n\t8 "+str(fibro.atom_mass)+\
        "\n\t9 "+str(fibro.atom_mass)+"\n\n"

        header += "Atoms\n\n"

        for i in range(len(fibro.atoms)):
            header += "\t"+" "+str(fibro.atoms[i][0])+" "+str(fibro.atoms[i][1])+" "+str(fibro.atoms[i][2])+" 0 "+\
            str(fibro.atoms[i][3])+" "+str(fibro.atoms[i][4])+" "+str(fibro.atoms[i][5])+" 0 0 0\n"

        header += "Bonds\n\n"

        for i in range(len(fibro.bonds)):
            header += "\t"+" "+str(i+1)+" 1 "+ str(fibro.bonds[i][0])+" "+\
                        str(fibro.bonds[i][1])+"\n"

        with open(fibro.topology_file,'w') as f:
            f.write("{}".format(header))
