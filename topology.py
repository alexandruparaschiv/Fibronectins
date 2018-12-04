#!/usr/bin/env python3

def write_topology_file(self):

    header = "LAMMPS Description\n\n"+"\t"+str(self.atomic_index)+" atoms\n"+"\t"+\
    str(len(self.bonds))+" bonds\n \t0 angles\n \t0 dihedrals\n \t0 impropers\n\n \t5 atom types\n \t1 bond types"+\
    "\n \t0 angle types\n \t0 dihedral types\n \t0 improper types\n\n \t"+str(-self.side_lengthx)+\
    " "+str(self.side_lengthx)+" xlo xhi\n \t"+str(-self.side_lengthy)+\
    " "+str(self.side_lengthy)+" ylo yhi\n \t"+str(-self.side_lengthz)+\
    " "+str(self.side_lengthz)+" zlo zhi\n \t"


    header += "\nMasses\n\n \t1 50.000\n\t2 50.000\n\t3 50.000\n\t4 50.000\n\t5 50.000\n\n"

    header += "Atoms\n\n"

    for i in range(len(self.atoms)):
        header += "\t"+" "+str(self.atoms[i][0])+" "+str(self.atoms[i][1])+" "+str(self.atoms[i][2])+" 0 "+\
        str(self.atoms[i][3])+" "+str(self.atoms[i][4])+" "+str(self.atoms[i][5])+" 0 0 0\n"

    header += "Bonds\n\n"

    for i in range(len(self.bonds)):
        header += "\t"+" "+str(i+1)+" 1 "+ str(self.bonds[i][0])+" "+\
                    str(self.bonds[i][1])+"\n"


    with open(self.topology_file,'w') as f:
        f.write("{}".format(header))
