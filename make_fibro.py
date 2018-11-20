#!/usr/bin/env python3

class Fibronectins(object):

    def __init__(self):

        self.side_lengthx = self.side_lengthy = 100
        self.side_lengthy = 100
        self.side_lengthz = self.side_lengthx
        self.domain_length = 10
        self.no_domains = 5
        self.domain_spacing = 2.5
        self.protein_beads = self.domain_length * self.no_domains
        self.fibronectin_separation = 15
        self.no_fibronectins = 1
        self.inner_bond_length = 1
        self.interaction_site_index = 0
        self.atoms = []
        self.molecules = []
        self.types = []
        self.bonds = []
        self.angles = []
        self.atomic_index = 0
        self.molecule_index = 0 
        self.bond_index = 0
        self.angle_index = 0
        self.bonded_atoms = []
        BEAD_TYPES = {"H":1,"He":2,"Li":3,"Be":4,"B":5}
        
        self.make_system()
        self.write_topology_file()
        self.write_lammps_file()

    def make_beta_domain(self,domain_index,bead_index,x,y,z):

        self.atomic_index += 1
        if (domain_index == 0) and (bead_index == 0):
            atom_type = 2
        elif  (domain_index == (self.no_domains-1) and bead_index == (self.domain_length-1)):
            atom_type = 3
        elif (domain_index%2 == 0) and bead_index  in (4,5,6):
            atom_type = 4
        elif (domain_index%2==1) and bead_index in (4,5,6):
            atom_type = 5
        else:
            atom_type = 1
        self.types.append(atom_type)
        self.molecules.append(self.molecule_index)
        self.atoms.append((self.atomic_index,x+0,y+domain_index*self.domain_spacing,z+bead_index*self.inner_bond_length))

        if domain_index != 0 or bead_index != 0:
            self.bond_index += 1 
            if (domain_index%(bead_index+1)):
                self.bonded_atoms.append((self.bond_index,1,self.atomic_index-1,self.atomic_index))
            else:
                self.bonded_atoms.append((self.bond_index,2,self.atomic_index-1,self.atomic_index))

        if domain_index != 0 or (bead_index != 0 and bead_index != 1):
            self.angle_index += 1
            self.angles.append((self.angle_index,1,self.angle_index,self.angle_index+1,self.angle_index+2))



        
                    
    def create_fibronectin(self,x,y,z):
        self.molecule_index += 1 
        for domain_index in range(self.no_domains):
            if domain_index % 2 == 0:
                for bead_index in range(self.domain_length):
                    self.make_beta_domain(domain_index,bead_index,x,y,z)
                    if bead_index == 2:
                        self.interaction_site_index += 1


            else:
                for bead_index in range(self.domain_length):
                    self.make_beta_domain(domain_index,bead_index,x,y,z)
                    if bead_index == 3:
                        self.interaction_site_index += 1

    def make_system(self):
        
        for i in range(5):
            for j in range(5):
                for k in range(5):
                    self.create_fibronectin(i*self.fibronectin_separation,j*self.fibronectin_separation,k*self.fibronectin_separation)




    def write_topology_file(self):
        
        header="LAMMPS Description\n\n"+"\t"+str(self.atomic_index)+" atoms\n"+"\t"+\
        str(len(self.bonded_atoms))+" bonds\n \t"+" "+str(self.angle_index)+" angles\n \t0 dihedrals\n \t0 impropers\n\n \t5 atom types\n \t2 bond types"+\
       "\n \t1 angle types\n \t0 dihedral types\n \t0 improper types\n\n \t"+str(-self.side_lengthx)+\
       " "+str(self.side_lengthx)+" xlo xhi\n \t"+str(-self.side_lengthy)+\
        " "+str(self.side_lengthy)+" ylo yhi\n \t"+str(-self.side_lengthz)+\
        " "+str(self.side_lengthz)+" zlo zhi\n \t"+"\nMasses\n\n \t1 10.000\n\t2 10.000\n\t3 10.000\n\t4 10.000\n\t5 10.000\n\n"+\
        "Atoms\n\n"

        self.atoms= sorted(self.atoms)
        for i in range(len(self.atoms)):
            header += "\t"+" "+str(self.atoms[i][0])+" "+str(self.molecules[i])+" "+str(self.types[i])+" 0 "+\
            str(self.atoms[i][1])+" "+str(self.atoms[i][2])+" "+str(self.atoms[i][3])+" 0 0 0\n"     
            
        header += "Bonds\n\n"
        
        for i in range(len(self.bonded_atoms)):
            header += "\t"+" "+str(self.bonded_atoms[i][0])+" "+str(self.bonded_atoms[i][1])+" "+ str(self.bonded_atoms[i][2])+" "+\
            str(self.bonded_atoms[i][3])+"\n"   

        header += "Angles\n\n"

        for i in range(len(self.angles)):
            header += "\t"+" "+str(self.angles[i][0])+" "+str(self.angles[i][1])+" "+ str(self.angles[i][2])+" "+\
            str(self.angles[i][3])+" "+str(self.angles[i][4])+"\n"   
            
        with open("topology_fibro.dat",'w') as top_file:
            top_file.write("{}".format(header))
            


    def write_lammps_file(self):


        with open("f.txt",'r') as f:
            force = float(f.readline())

        header = "units\tlj\n"+"boundary\tp p p\n"+ "atom_style\tfull\n"+ "read_data\t topology_fibro.dat\n"+ "neighbor\t0.3 bin\n"+\
        "neigh_modify\tevery 1 delay 1\n"+ "pair_style\tlj/cut 3.0\n"+ "pair_coeff\t* * 1 1 1.0 \n"+"pair_coeff       4 5 3 1 2.50\n"+\
        "pair_modify\tshift yes\n"+\
        "bond_style hybrid fene harmonic\n" + "bond_coeff 1 fene 30.0 1.5 1.0 1.0\n"+ "bond_coeff 2 harmonic 80.0 1.2\n"+"angle_style\tharmonic\n"+"angle_coeff\t 1 5.0 180\n"\
         "group\tfibronectins type 1\n"+"group\tright_ends type 2\n"+"group\tleft_ends type 3\n"+"neigh_modify exclude molecule fibronectins\n""velocity\tfibronectins create 1.0 1\n"+\
        "dump            id all xyz 1000 fibro.xyz\n"+"comm_modify     cutoff 50\n"+\
        "dump_modify     id sort id\n"+"thermo          1000\n"+"thermo_modify   flush yes\n"+\
        "compute cc1 all chunk/atom molecule\n"+"compute myChunk all com/chunk cc1\n"+"fix 10 all ave/time 1000 1 1000 c_myChunk[*] file coms.out mode vector\n"+\
        "minimize 1.0e-4 1.0e-6 100 1000\n"+"fix             1 all nvt temp 1 1 1\n"+"fix             2 all langevin 1 1 10000 12345\n"+\
        "timestep     0.008\n"+"run 200000\n"+\
        "fix  pull    right_ends smd cfor  "+str(force)+" tether 100.0 0.0 0.0 0.0\n"+"fix  pull2    left_ends smd cfor  "+str(force)+" tether -100.0 0.0 0.0 0.0\n"+\
        "run             5000000\n"
        with open("in.fibronectins",'w') as md_file:
            md_file.write("{}".format(header))


if __name__ == "__main__":

    fibronectins = Fibronectins()
