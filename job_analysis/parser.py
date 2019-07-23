import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import time
import re 
import distutils.util
from multiprocessing import Process,Pool

parser = argparse.ArgumentParser(description='process the resulting file.')
parser.add_argument('--movie_name', nargs="?",type=str, const="fibro.xyz",
    help='the trajectory file name')
parser.add_argument('--last_frame', '-lf', type=distutils.util.strtobool, default='true',
    help='only processes the last frame')
parser.add_argument('--verbose', '-v', type=distutils.util.strtobool, default='true')
parser.add_argument('--nematic_cutoff','-n', type=float, default=0.5,
    help='stores the nematic parameter cutoff')
parser.add_argument('--radial_cutoff', '-r', type=float, default=1.45,
    help='stores the Lennard-Jones interaction cutoff')
args = parser.parse_args()

try:
    from analyser import FibrilAnalyser
    from colvars import Colvars
except ImportError:
    print("Could not import the required classes!")


class Parser:

    def __init__(self,filename):

        self.filename = filename
        self.initial_time = time.time()
        self.progress = 0
        self.network = []
        self.process_pool = Pool(processes=4)

        self._count_atoms()
        self._track_progress()

    def _count_atoms(self):

        with open(self.filename,'r') as infile:
            try:
                self._atoms_number = int(infile.readline())
            except:
                print("Could not read the number of atoms,please make sure the .xyz file has the correct format.")
                exit()

    def _make_last_movie_frame(self):

        os.system("rm last_frame.xyz")
        os.system("rm last.xyz")
        cmd = "tail -"+str(self._atoms_number)+" fibro.xyz >> last_frame.xyz "
        with open("last_frame.xyz",'w') as outfile:
            outfile.write("{}\nAtoms. Timestep: 0\n".format(self._atoms_number))
        os.system(cmd)
        cmd = "cat last_frame.xyz last_frame.xyz >> last.xyz"
        os.system(cmd)

    def _track_progress(self):

        if os.path.exists("log.lammps"):
            os.system("tail -1 log.lammps >> progress.txt")
            with open("progress.txt") as infile:
                try:
                    self.progress = int(re.findall(r'\d+', infile.readline())[0])
                except:
                    print("Could not find the last frame, run the simulation for longer!")
            os.system("rm progress.txt")
        else:
            print("Warning: Could not find the LAMMPS log file.")


    def go_through_frames(self,filename):

        frame_index = 0
        types = []
        xs = []; ys = []; zs = []
        jobs = []
        fibril_analyser = FibrilAnalyser()


        with open(filename, 'r') as infile:
            for item in infile:
                if len(item.split()) == 1 and types:

                    frame_index += 1
                    types = np.asarray(types);xs=np.asarray(xs);ys=np.asarray(ys);zs=np.asarray(zs)
                    if args.verbose and not args.last_frame:
                        print("Processing frame......{}".format(frame_index))
                    self.process_pool.apply_async(fibril_analyser.analyse_fibril_distribution,
                    args=(args.nematic_cutoff,args.radial_cutoff,types,xs,ys,zs))
                    types=[];xs=[];ys=[];zs=[]

                if len(item.split()) > 3:

                    atom_type,x,y,z=item.split()
                    atom_type = str(atom_type)
                    x, y, z = float(x), float(y), float(z)
                    types.append(int(atom_type))
                    xs.append(x);ys.append(y);zs.append(z)

        self.process_pool.close()
        self.process_pool.join()


if __name__ == "__main__":

    parser = Parser(args.movie_name)
    parser._count_atoms()
    parser._track_progress()
    if args.last_frame:
        print("analysing last frame")
        parser._make_last_movie_frame()
        parser.go_through_frames("last.xyz")
    else:
        parser.go_through_frames(args.movie_name)
    print("Successfully analysed all frames!")
    print("{} seconds elapsed.".format(time.time()-parser.initial_time))
