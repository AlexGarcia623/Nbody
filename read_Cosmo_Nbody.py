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
from matplotlib.animation import FuncAnimation

class Particle:
    def __init__(self,out_directory):
        self.directory = out_directory
        self.positions = None

        self.files = sorted([f for f in os.listdir(out_directory) if f.startswith('particle_positions')])
    
        self.filename = out_directory.replace('./output/','').replace('/','.txt')

        self._get_params(self.filename)
        self.get_positions()
    
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
        Ndim = 3 # X, Y, Z, Mass
        self.Nstep = int(self.params['N_steps'])
        self.positions = np.empty((self.Nbody,Ndim,self.Nstep))
        
        for f_index, file in enumerate(self.files):
            line = self._read_particle_positions(os.path.join(self.directory, file))
            for index in range(len(line)):
                self.positions[index,:,f_index] = line[index]
        return self.positions

    def _animate(self,i):
        self.ax.clear()
        self.ax.set_ylim(-15,15)
        self.ax.set_xlim(-15,15)
        for index in range(self.Nbody):
            pos = self.positions[index,:,i]
            self.ax.scatter(pos[0], pos[1], marker='o')
    
    def make_movie(self,interval=20,repeat=False):
        fig, self.ax = plt.subplots()
        ani = FuncAnimation(fig, self._animate, frames=self.Nstep, interval=interval, repeat=repeat)
        video = ani.to_html5_video()
        html = display.HTML(video)
        display.display(html)
        plt.close()