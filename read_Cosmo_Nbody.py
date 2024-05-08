'''Read in particle data from Cosmo Nbody
#
This module, written by Alex Garcia for AST5470 @ UVA, Spring 2024,
#
Usage:
    - Import the necessary packages before using this module.
    - Use the functions provided in this module for specific tasks.
#
Functions Provided:
'''
__author__ = "Alex Garcia"
__copyright__ = "Copyright 2024, The Authors"
__credits__ = ["Alex Garcia"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Alex Garcia"
__email__ = "alexgarcia@virginia.edu"

import re
import os
import csv
import numpy as np
from IPython import display
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

class Particle:
    def __init__(self,out_directory):
        self.directory = out_directory
        self.positions = None

        self.files = sorted([f for f in os.listdir(out_directory) if f.startswith('particle_positions')])
        if len(self.files) < 1:
            raise Exception("You did not run your simulation yet (or did not save output)")
    
        self.filename = out_directory.replace('./output/','').replace('/','.txt')

        self._get_params(self.filename)
        self.get_positions()

        print(f'Found: ')
        sep = max(len(str(self.Nbody)), len(str(self.Nsnaps))) + 1
        print(f"\t{self.Nbody:<{sep}}particles")
        print(f"\t{self.Ndim:<{sep}}dimensions")
        print(f"\t{self.Nsnaps:<{sep}}snapshots")
    
    def _read_particle_positions(self,filename):
        positions = []
        with open(filename, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)  # Skip header row
            for row in csvreader:
                positions.append([float(row[0]), float(row[1]), float(row[2])])
        return positions

    def _get_params(self, filename):        
        self.params = {}
        with open(filename, 'r') as file:
            for line in file:
                # Use regular expressions to match lines with parameter values
                match = re.match(r'^(\w+)\s*:\s*(.+)$', line)
                if match:
                    key = match.group(1)
                    value = match.group(2)
                    self.params[key] = value
    
    def get_positions(self):
        self.Nbody = int(self.params['N_particles'])
        self.Ndim = 3 # X, Y, Z, Mass
        self.Nstep = int(self.params['N_steps']) 
        self.cadence = int(self.params["Snapshot_cadence"])

        self.Nsnaps = len(self.files)
        
        self.positions = np.empty((self.Nbody,self.Ndim,self.Nsnaps))
        
        for f_index, file in enumerate(self.files):
            line = self._read_particle_positions(os.path.join(self.directory, file))
            for index in range(len(line)):
                self.positions[index,:,f_index] = line[index]
        return self.positions

    def _animate(self, i):
        self.ax.clear()
        L_BOX = float(self.params['L_box'])
        self.ax.set_ylim(0, L_BOX)
        self.ax.set_xlim(0, L_BOX)
        
        # Assuming self.positions is a 3D array with shape (Nbody, 3, num_frames)
        positions_at_frame_i = self.positions[:, :, i]
        self.ax.scatter(positions_at_frame_i[:, self.projection1], 
                        positions_at_frame_i[:, self.projection2], marker='o')

    def _animate3d(self, i):
        self.ax.clear()
        L_BOX = float(self.params['L_box'])
        self.ax.set_xlim(0, L_BOX)
        self.ax.set_ylim(0, L_BOX)
        self.ax.set_zlim(0, L_BOX)

        angle = i * 2 * np.pi / 36  # Adjust the rotation speed if needed
        self.ax.view_init(elev=30, azim=angle)  # Set the elevation and azimuth angles
        
        # Assuming self.positions is a 3D array with shape (Nbody, 3, num_frames)
        positions_at_frame_i = self.positions[:, :, i]
        self.ax.scatter(positions_at_frame_i[:, 0], 
                        positions_at_frame_i[:, 1], 
                        positions_at_frame_i[:, 2], marker='o')


    
    def make_movie(self,interval=40,repeat=False,projection='xy'):
        projection = projection.upper()
        different_projections = {
            "XY":[0,1],
            "XZ":[0,2],
            "YZ":[1,2],
            "3D":[0,1]
        }
        self.projection1, self.projection2 = different_projections[projection]
        if projection.upper() == "3D":
            fig = plt.figure()
            self.ax = fig.add_subplot(111, projection='3d')
            ani = FuncAnimation(fig, self._animate3d, frames=self.Nsnaps, interval=interval, repeat=repeat)
        else:
            fig, self.ax = plt.subplots()
            ani = FuncAnimation(fig, self._animate, frames=self.Nsnaps, interval=interval, repeat=repeat)
        video = ani.to_html5_video()
        html = display.HTML(video)
        display.display(html)
        plt.close()