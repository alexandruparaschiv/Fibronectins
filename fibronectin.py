#!/usr/bin/env python3

import math
from system import make_system
from topology import write_topology_file
from settings import write_lammps_file

class Fibronectins:

    def __init__(self):

        # sets the simulation parameters
        self.side_lengthx = self.side_lengthy = self.side_lengthz = 200
        self.timestep = 0.008
        self.duration = 5000000
        self.atom_mass = 1

        # sets the parameters of a fibronectin monomer
        self.no_rods = 3
        self.beads_per_rod = 10
        self.beads_per_loop = 7
        self.bead_spacing = 1
        self.patch_offset = 1.15
        self.loop_epsilon = 1
        self.patch_diameter = 0.5
        self.hydrogen_bond_epsilon = 50
        self.hydrogen_bond_cutoff = 1.12

        # contains information about the monomer topology
        self.atoms = []
        self.molecules = []
        self.bonds = []
        self.atomic_index = 0
        self.molecular_index = 0

        # controls the output
        self.trajectory_file = "fibro.xyz"
        self.topology_file = "topology.dat"
        self.dynamics_file = "in.fibronectins"

        make_system(self)
        write_topology_file(self)
        write_lammps_file(self)
        print("The system was initialised successfully!")


    def make_monomer(self,x,y,z):


        for i in range(self.no_rods):

            bead_type = 2 if i % 2 == 1 else 3
            self.molecular_index += 1

            for j in range(self.beads_per_rod):

                self.atomic_index += 1

                core_x = x
                core_y = y+(i+2)*self.no_rods
                core_z = z+j*self.bead_spacing
                self.atoms.append((self.atomic_index,self.molecular_index,1,core_x,core_y,core_z))

                patch_x = x
                patch_y = core_y-self.patch_offset
                patch_z = core_z
                self.atoms.append((self.atomic_index+self.beads_per_rod,self.molecular_index,bead_type,
                                    patch_x,patch_y,patch_z))

                if j != self.beads_per_rod - 1 :
                    self.bonds.append((self.atomic_index,self.atomic_index+1))
                self.bonds.append((self.atomic_index,self.atomic_index+self.beads_per_rod))


            self.atomic_index += self.beads_per_rod

            for k in range(self.beads_per_loop):

                self.atomic_index += 1
                orientation = 1 if not i%2 else -1
                theta = math.pi/(self.beads_per_loop-1)

                if orientation == 1:

                    loop_x = x
                    loop_y = core_y - 2 * self.bead_spacing*(1-math.cos(k*theta))
                    loop_z = -1 + z - 2 * self.bead_spacing*math.sin(k*theta)

                    self.atoms.append((self.atomic_index, self.molecular_index, 4, loop_x,loop_y,loop_z))
                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-2*self.beads_per_rod))

                    elif k < self.beads_per_loop - 1:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))
                        if i != 0:
                            self.bonds.append((self.atomic_index,self.atomic_index-4*self.beads_per_rod-self.beads_per_loop-k))

                    #code to add one extra loop to the end of the monomer in case of odd number of rods
                    """if i == self.no_rods-1:

                        loop_x = x
                        loop_y = core_y + 2*self.bead_spacing*(1-math.cos(k*theta))
                        loop_z = -1 + z + 2*self.bead_spacing*math.sin(k*theta)+self.beads_per_rod*self.bead_spacing

                        self.atoms.append((self.atomic_index, self.molecular_index, 5, loop_x, loop_y,loop_z))
                        if k == 0:
                            self.bonds.append((self.atomic_index,self.atomic_index-self.beads_per_rod-1))

                        elif k < self.beads_per_loop - 1:
                            self.bonds.append((self.atomic_index,self.atomic_index-1))

                        else:
                            self.bonds.append((self.atomic_index,self.atomic_index-1))
                            self.bonds.append((self.atomic_index,self.atomic_index-3*self.beads_per_rod-2*self.beads_per_loop))"""

                    #new code block ends here

                else:

                    loop_x = x
                    loop_y = core_y - 2*self.bead_spacing*(1-math.cos(k*theta))
                    loop_z = z + 2*self.bead_spacing*math.sin(k*theta)+self.beads_per_rod*self.bead_spacing

                    self.atoms.append((self.atomic_index, self.molecular_index, 5, loop_x, loop_y,loop_z))

                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-self.beads_per_rod-1))

                    elif k < self.beads_per_loop - 1:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))
                        self.bonds.append((self.atomic_index,self.atomic_index-3*self.beads_per_rod-2*self.beads_per_loop))


            #add protein tail here
            if (i == self.no_rods-1) and self.no_rods%2 == 1:
                for k in range(self.beads_per_loop):

                    self.atomic_index += 1

                    loop_x = core_x
                    loop_y = core_y + 2*self.bead_spacing*(1-math.cos(k*theta))
                    loop_z = core_z - self.beads_per_rod*self.bead_spacing + 2*self.bead_spacing*math.sin(k*theta)+self.beads_per_rod*self.bead_spacing

                    self.atoms.append((self.atomic_index, self.molecular_index, 4, loop_x,loop_y,loop_z))
                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-self.beads_per_rod-self.beads_per_loop-1))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))





if __name__ == "__main__":
    fibronectins = Fibronectins()
