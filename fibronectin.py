#!/usr/bin/env python3

import math


class Fibronectins:

    def __init__(self):

        self.side_lengthx = self.side_lengthy = self.side_lengthz = 200

        self.no_rods = 9
        self.beads_per_rod = 10
        self.beads_per_loop = 9
        self.bead_spacing = 1
        self.atoms = [];
        self.molecules = []
        self.atomic_index = 0
        self.molecular_index = 0
        self.bonded_atoms = []
        self.trajectory_file = "fibro.xyz"

        for i in range(1):
            for j in range(1):
                for k in range(1):
                    self.make_rigid_monomer(5*i,5*j,5*k)



        self.write_topology_file()
        self.write_lammps_file()

        #tester = Tests()
        #tester.test(self.bonded_atoms)

    def make_rigid_monomer(self,x,y,z):


        for i in range(self.no_rods):

            bead_type = 2 if i % 2 == 1 else 3
            self.molecular_index += 1

            for j in range(self.beads_per_rod):

                self.atomic_index += 1

                self.atoms.append((self.atomic_index,self.molecular_index,1,x,y+i*0.5*self.no_rods,z+j*self.bead_spacing))
                self.atoms.append((self.atomic_index+self.beads_per_rod,self.molecular_index,bead_type,x,y+0.15+i*0.5*self.no_rods-1,z+j*self.bead_spacing))

                if j != self.beads_per_rod - 1 :
                    self.bonded_atoms.append((self.atomic_index,self.atomic_index+1))
                self.bonded_atoms.append((self.atomic_index,self.atomic_index+self.beads_per_rod))



            self.atomic_index += self.beads_per_rod

            for k in range(self.beads_per_loop):

                self.atomic_index += 1
                orientation = 0
                if i % 2 == 0:
                    orientation = 1
                else:
                    orientation = -1
                theta = math.pi/(self.beads_per_loop-1)

                if orientation == 1:


                    self.atoms.append((self.atomic_index, self.molecular_index, 4, x, y+i*0.5*self.no_rods - 2*self.bead_spacing*(1-math.cos(k*theta)),-1+z-2*self.bead_spacing*math.sin(k*theta)))
                    if k == 0:
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-2*self.beads_per_rod))

                    elif k < self.beads_per_loop - 1:
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-1))

                    else:
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-1))
                        if i != 0:
                            print((self.atomic_index,self.atomic_index-4*self.beads_per_rod-self.beads_per_loop-k))
                            self.bonded_atoms.append((self.atomic_index,self.atomic_index-4*self.beads_per_rod-self.beads_per_loop-k))

                else:

                    self.atoms.append((self.atomic_index, self.molecular_index, 5, x, y+i*0.5*self.no_rods - 2*self.bead_spacing*(1-math.cos(k*theta)), 1+z+2*self.bead_spacing*math.sin(k*theta)+self.beads_per_rod*self.bead_spacing))

                    if k == 0:
                        print((self.atomic_index,self.atomic_index-self.beads_per_rod-1))
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-self.beads_per_rod-1))

                    elif k < self.beads_per_loop - 1:
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-1))

                    else:
                        print((self.atomic_index,self.atomic_index-1),(self.atomic_index,self.atomic_index-4*self.beads_per_rod))
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-1))
                        self.bonded_atoms.append((self.atomic_index,self.atomic_index-4*self.beads_per_rod-k))



    def write_topology_file(self):

        header = "LAMMPS Description\n\n"+"\t"+str(self.atomic_index)+" atoms\n"+"\t"+\
        str(len(self.bonded_atoms))+" bonds\n \t0 angles\n \t0 dihedrals\n \t0 impropers\n\n \t5 atom types\n \t1 bond types"+\
        "\n \t0 angle types\n \t0 dihedral types\n \t0 improper types\n\n \t"+str(-self.side_lengthx)+\
        " "+str(self.side_lengthx)+" xlo xhi\n \t"+str(-self.side_lengthy)+\
        " "+str(self.side_lengthy)+" ylo yhi\n \t"+str(-self.side_lengthz)+\
        " "+str(self.side_lengthz)+" zlo zhi\n \t"+"\nMasses\n\n \t1 50.000\n\t2 50.000\n\t3 50.000\n\t4 50.000\n\t5 50.000\n"+"Atoms\n\n"

        for i in range(len(self.atoms)):
            header += "\t"+" "+str(self.atoms[i][0])+" "+str(self.atoms[i][1])+" "+str(self.atoms[i][2])+" 0 "+\
            str(self.atoms[i][3])+" "+str(self.atoms[i][4])+" "+str(self.atoms[i][5])+" 0 0 0\n"

        header += "Bonds\n\n"

        for i in range(len(self.bonded_atoms)):
            header += "\t"+" "+str(i+1)+" 1 "+ str(self.bonded_atoms[i][0])+" "+\
            str(self.bonded_atoms[i][1])+"\n"

        with open("topology_fibronectins.dat",'w') as f:
            f.write("{}".format(header))


    def write_lammps_file(self):

        header = "units\tlj\n"+"boundary\tp p p\n"+ "atom_style\tfull\n"+ "read_data\t topology_fibronectins.dat\n"+ "neighbor\t0.3 bin\n"+\
        "neigh_modify\tevery 1 delay 1\n"+ "pair_style\tlj/cut 3.0\n"+ "pair_coeff\t* * 1 100 1.50 \npair_coeff\t2 3 10 1 1.5\n"+\
        "pair_modify\tshift yes\n"+"#neigh_modify exclude molecule all\n"+"bond_style\t harmonic\n"+"bond_coeff\t1 100 1.0"+"\n\n"+"group\tfibronectins type 1 2 3\n"+"velocity\tfibronectins create 1.0 1\n"+\
        "group loops type 4 5\n"+"dump            id all xyz 100 "+self.trajectory_file +"\n"+"comm_modify     cutoff 50\n"+\
        "dump_modify     id sort id\n"+"#restart         10000 restart.dat\n"+"thermo          100\n"+"thermo_modify   flush yes\n"+\
        "minimize 1.0e-4 1.0e-6 100 1000\n"+"fix             1 fibronectins rigid/nve molecule\n"+"fix             2 fibronectins langevin 1 1 10000 12345\n"+\
        "fix       3 loops nve\n"+"timestep        0.008\n"+"run             5000000\n"
        with open("in.fibronectins",'w') as f:
            f.write("{}".format(header))

if __name__ == "__main__":
    fibronectins = Fibronectins()
