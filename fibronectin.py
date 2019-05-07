#!/usr/bin/env python3

import math
from system import System
from topology import Topology
from dynamics import Dynamics

class Fibronectins:

    def __init__(self):

        # sets the simulation parameters
        self.side_lengthx = self.side_lengthy = self.side_lengthz = 300
        self.side_lengthx = 300
        self.side_lengthx = 100
        self.timestep = 0.008
        self.duration = 5000000
        self.atom_mass = 1

        # sets the parameters of a fibronectin monomer
        self.no_rods = 14
        self.beads_per_rod = 10
        self.beads_per_loop = 8
        self.bead_spacing = 1
        self.patch_offset = 1.25
        self.loop_epsilon = 1
        self.rod_spacing = 3.65
        self.patch_diameter = 1.0
        self.hydrogen_bond_epsilon = 5
        self.hydrogen_bond_cutoff = 1.30

        # contains information about the monomer topology
        self.atoms = []
        self.molecules = []
        self.bonds = []
        self.monomer_chunks = []
        self.atomic_index = 0
        self.molecular_index = 0
        self.x_spacing = 1.5
        self.y_spacing = 50
        self.z_spacing = 30

        # controls the output
        self.trajectory_file = "fibro.xyz"
        self.topology_file = "topology.dat"
        self.dynamics_file = "in.fibronectins"

        System.make_lattice(self)
        Topology.write_topology_file(self)
        Dynamics.write_dynamics_file(self)
        print("The system was initialised successfully!")

    def make_monomer(self,x,y,z):

        self.monomer_chunks.append(self.atomic_index)
        for i in range(self.no_rods):

            if i < 4 or i > 11:
                bead_type = 8 if i % 2 == 1 else 9
            else:
                bead_type = 2 if i % 2 == 1 else 3
            self.molecular_index += 1

            for j in range(self.beads_per_rod):

                self.atomic_index += 1

                core_x = x
                core_y = y+i*self.rod_spacing
                core_z = z+j*self.bead_spacing
                self.atoms.append((self.atomic_index, self.molecular_index, 1, core_x, core_y, core_z))

                # place the hydrogen bonding patches on alternate sides of the chain
                patch_x = x
                patch_y = core_y-2*(2.5-bead_type)*self.patch_offset
                patch_z = core_z

                if bead_type in [8,9]:
                    patch_y = core_y-2*(8.5-bead_type)*self.patch_offset
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

                    loop_bead_type = 6 if (k == self.beads_per_loop-1 and i==0) or (k==0 and i == self.no_rods)  else 4
                    self.atoms.append((self.atomic_index, self.molecular_index, loop_bead_type, loop_x,loop_y,loop_z))
                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-2*self.beads_per_rod))

                    elif k < self.beads_per_loop - 1:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))
                        if i != 0:
                            self.bonds.append((self.atomic_index,self.atomic_index-4*self.beads_per_rod-self.beads_per_loop-k))


                else:

                    loop_x = x
                    loop_y = core_y - 2*self.bead_spacing*(1-math.cos(k*theta))
                    loop_z = z + 2*self.bead_spacing*math.sin(k*theta)+self.beads_per_rod*self.bead_spacing

                    loop_bead_type = 5

                    self.atoms.append((self.atomic_index, self.molecular_index, loop_bead_type, loop_x, loop_y,loop_z))

                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-self.beads_per_rod-1))

                    elif k < self.beads_per_loop - 1:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))
                        self.bonds.append((self.atomic_index,self.atomic_index-3*self.beads_per_rod-2*self.beads_per_loop))

            # protein tail is added here
            if (i == self.no_rods-1) and self.no_rods%2 == 1:
                for k in range(self.beads_per_loop):

                    self.atomic_index += 1

                    loop_x = core_x
                    loop_y = core_y + 2*self.bead_spacing*(1-math.cos(k*theta))
                    loop_z = core_z - self.beads_per_rod*self.bead_spacing + 2*self.bead_spacing*math.sin(k*theta)+self.beads_per_rod*self.bead_spacing

                    loop_bead_type = 4
                    self.atoms.append((self.atomic_index, self.molecular_index, loop_bead_type, loop_x,loop_y,loop_z))
                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-self.beads_per_rod-self.beads_per_loop-1))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))

            if (i == self.no_rods-1) and self.no_rods%2 == 0:
                for k in range(self.beads_per_loop):

                    self.atomic_index += 1

                    loop_x = core_x
                    loop_y = core_y + 2*self.bead_spacing*(1-math.cos(k*theta))
                    loop_z = core_z - 2*self.beads_per_rod*self.bead_spacing + 2*self.bead_spacing*math.sin(-k*theta)+self.beads_per_rod*self.bead_spacing

                    loop_bead_type = 4 if k != self.beads_per_loop - 1 else 7
                    self.atoms.append((self.atomic_index, self.molecular_index, loop_bead_type, loop_x,loop_y,loop_z))
                    if k == 0:
                        self.bonds.append((self.atomic_index,self.atomic_index-2*self.beads_per_rod-self.beads_per_loop))

                    else:
                        self.bonds.append((self.atomic_index,self.atomic_index-1))

        self.monomer_chunks.append(self.atomic_index) 


if __name__ == "__main__":
    fibronectins = Fibronectins()
