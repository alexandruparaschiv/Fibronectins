#!/usr/bin/env python3


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
