#!/usr/bin/env python3


def write_lammps_file(self):

    header = "units\tlj\n"+"boundary\tp p p\n"+ "atom_style\tfull\n\n"+\
    "read_data\t"+self.topology_file+"\n"+"#restart\t10000 restart.dat\n\n"+"neighbor\t0.3 bin\n"+\
    "neigh_modify\tevery 1 delay 1\n"+"neigh_modify\texclude molecule all\n\n"+\
    "pair_style\tlj/cut 3.0\n"+ "pair_coeff\t* * 1 1 1.0 \npair_coeff\t2 3 10 1 1.5\n"+\
    "pair_modify\tshift yes\n\n"+"bond_style\tharmonic\n"+"bond_coeff\t1 100 1.0"+"\n\n"+\
    "group\tfibronectins type 1 2 3\n"+"group loops type 4 5\n\n"+"velocity\tfibronectins create 1.0 1\n"+\
    "velocity\tloopss create 1.0 1\n"+\
    "dump\tid all xyz 100 "+self.trajectory_file +"\n"+"comm_modify\tcutoff 50\n"+\
    "dump_modify\tid sort id\n"+"\nthermo\t100\n"+"thermo_modify\tflush yes\n"+\
    "minimize 1.0e-4 1.0e-6 100 1000\n\n"+\
    "fix\t1 fibronectins rigid/nve molecule\n"+"fix\t2 fibronectins langevin 1 1 10000 12345\n"+\
    "fix\t3 loops nve\n\n"+\
    "timestep        0.008\n"+"run             5000000\n"
    with open("in.fibronectins",'w') as f:
        f.write("{}".format(header))
