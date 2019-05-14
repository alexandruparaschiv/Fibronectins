#!/usr/bin/env python3

class Dynamics:

    def __init__(self):
        pass

    @staticmethod
    def write_dynamics_file(fibro):

        # defines units, data files, atom styles and neighbor lists
        header = "units\tlj\n" + "boundary\tp p p\n" + "atom_style\tfull\n\n" + \
                 "read_data\t" + fibro.topology_file + "\n" + "#restart\t10000 restart.dat\n\n" + "neighbor\t0.3 bin\n" + \
                 "neigh_modify\tevery 1 delay 1\n" + "neigh_modify\texclude molecule all\n\n"

        # defines the pairwise interactions   and the bonds
        header += "pair_style\tlj/cut 3.0\n" + "pair_coeff\t* * 1 " + str(fibro.loop_epsilon) + " 1.0 \n" + \
                 "pair_coeff\t 2 2 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t2 3 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t 3 3 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" 
        header += "pair_coeff\t 8 8 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t8 9 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t 8 9 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" 
        header += "pair_coeff\t 2 8 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t3 8 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t 2 9 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_coeff\t 3 9 " + str(fibro.hydrogen_bond_epsilon) + " " + str(fibro.patch_diameter) + " " + str(
            fibro.hydrogen_bond_cutoff) + "\n" + \
                 "pair_modify\tshift yes\n\n" + "bond_style\tfene\n" + "bond_coeff\t1 30.0 1.5 1.0 1.0\n\n"
        # defines the groups and output files
        header += "group\tfibronectins type 1 2 3 8 9\n" + "group loops type 4 5 6 7\n\n" + "group left_ends type 6\n" + "group right_ends type 7\n" + \
                 "group molecules union fibronectins loops\n"+"velocity\tfibronectins create 1.0 1\n" + "velocity\tloops create 1.0 1\n" + \
                 "dump\tid all xyz 100 " + fibro.trajectory_file + "\n" + "comm_modify\tcutoff 50\n" + \
                 "dump_modify\tid sort id\n" + "\nthermo\t100\n" + "thermo_modify\tflush yes\n\n"
        
       # defines the fixes and  simulation time steps and duration
        header += "fix\t1 fibronectins rigid/nve molecule\n" + "fix\t2 fibronectins langevin 1 1 10000 12345\n" + \
                 "fix\t3 loops nvt temp 1.0 1.0 1\n" +"minimize 1.0e-4 1.0e-6 100 1000\n\n"  +"run\t"+str(fibro.equilibration)+" \n"+ \
                 "fix       Colvars all colvars config.colvars\n" + "minimize 1.0e-4 1.0e-6 100 1000\n\n" + \
                 "timestep\t" + str(fibro.timestep) + "\n" + str(fibro.duration) + "\n"

        with open("in.fibronectins", 'w') as f:
            f.write("{}".format(header))
