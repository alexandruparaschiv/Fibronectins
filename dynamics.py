#!/usr/bin/env python3

class Dynamics:

    def __init__(self):
        pass



    @staticmethod
    def write_dynamics_file(fibro):

        header = "units\tlj\n" + "boundary\tp p p\n" + "atom_style\tfull\n\n" + \
                 "read_data\t" + fibro.topology_file + "\n" + "#restart\t10000 restart.dat\n\n" + "neighbor\t0.3 bin\n" + \
                 "neigh_modify\tevery 1 delay 1\n" + "neigh_modify\texclude molecule all\n\n" + \
                 "pair_style\tlj/cut 3.0\n" + "pair_coeff\t* * 1 " + str(fibro.loop_epsilon) + " 1.0 \n" + \
                 "pair_coeff\t 2 2 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t2 3 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t 3 3 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_modify\tshift yes\n\n" + "bond_style\tfene\n" + "bond_coeff\t1 30.0 1.5 1.0 1.0" + "\n\n" + \
                 "group\tfibronectins type 1 2 3\n" + "group loops type 4 5 6 7\n\n" + "group left_ends type 6\n" + "group right_ends type 7\n" + "velocity\tfibronectins create 1.0 1\n" + \
                 "velocity\tloops create 1.0 1\n" + \
                 "dump\tid all xyz 100 " + fibro.trajectory_file + "\n" + "comm_modify\tcutoff 50\n" + \
                 "dump_modify\tid sort id\n" + "\nthermo\t100\n" + "thermo_modify\tflush yes\n\n" + \
                 "fix\t1 fibronectins rigid/nve molecule\n" + "fix\t2 fibronectins langevin 1 1 10000 12345\n" + \
                 "fix\t3 loops nvt temp 1.0 1.0 1\n" + "fix       Colvars all colvars config.colvars\n" + "minimize 1.0e-4 1.0e-6 100 1000\n\n" + \
                 "timestep\t" + str(fibro.timestep) + "\n" + "run\t" + str(fibro.duration) + "\n"

        with open("in.fibronectins", 'w') as f:
            f.write("{}".format(header))



