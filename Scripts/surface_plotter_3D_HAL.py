import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import os
import sys
import matplotlib
import scipy
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import copy

def randrange(n, vmin, vmax):
    return (vmax-vmin)*np.random.rand(n) + vmin

def split_process(lines):
#processes from csv format to a float format...
    def try_float(entry):
        try:
            outentry = float(entry)
            return outentry
        except:
            return None


    splitlines = [' '.join(line.split()) for line in lines[1:]]
    splitlines = [line.split(',') for line in splitlines]
    splitlines = [[try_float(entry) for entry in line[1:]] for line in splitlines]   #takes out index
    splitlines = [[entry for entry in line if not entry == None] for line in splitlines]
    splitlines = [line for line in splitlines if not line == []]
    return splitlines


def abs_process(splitlines, index_array):
#absolute value processing
    splitlines = [entry for entry in splitlines if len(entry) >= (2+len(index_array))]
    max_entry_length = max([len(entry) for entry in splitlines])

    splitlines = [entry[:2] + [entry[index] for index in index_array] for entry in splitlines ]
    return splitlines


def stat_print(split_inlines, index):
#prints out the proportion of points where the indexed variable is most important
    abslines = [[abs(subentry) for subentry in entry[2:]] for entry in split_inlines]
    maxcounter = 0
    for entry in abslines:
        if max(entry) == entry[index-2]:
            maxcounter += 1
    print(float(maxcounter)/float(len(abslines)))


extension = ''
interact_string = sys.argv[1]
infiles = [f for f in os.listdir(os.path.join('.', extension)) if '_HAL.csv' in f and not 'Prevalences' in f and len(f.split('_')) == 4] #and interact_string in f]
interact_names = infiles
interact_names.sort()
interact_names = interact_names
#interact_names = [entry for entry in interact_names if len(entry.split('_')) == 2]
print(len(infiles))

print(interact_names)

print_counter = 0
for interact_string in interact_names:
    print(interact_string)
    print_counter += 1
    print(str(print_counter) + ' out of ' + str(len(interact_names)))
    grids = []

    sub_infiles = [interact_string]
    #print(sub_infiles)
    #q()
    
    for infile_name in sub_infiles:
        infile = open(infile_name, 'r')
        inlines = infile.readlines()
        infile.close()

        inlines = split_process(inlines)
        yield_grid = inlines
        grids = yield_grid

    grids = np.array(grids)
    print(grids.shape)
#    deviations = np.mean(grids, axis=0)
    res_up_factor = 1.0
    deviations = grids # for deviations

    yield_max = np.max(deviations[:,-1])
    yield_min = np.min(deviations[:,-1])

    minx = np.min(deviations[:,0])
    maxx = np.max(deviations[:,0])
    miny = np.min(deviations[:,1])
    maxy = np.max(deviations[:,1])

    xs = deviations[:,0]
    ys = deviations[:,1]
    zs = deviations[:,2]
    the_fourth_dimension = deviations[:,-1]

    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111,projection='3d')

    colors = cm.coolwarm((the_fourth_dimension-min(the_fourth_dimension))/(max(the_fourth_dimension)-min(the_fourth_dimension)))

    colmap = cm.ScalarMappable(cmap=cm.coolwarm)
    colmap.set_array(the_fourth_dimension)

    yg = ax.scatter(xs, ys, zs, c=colors, marker='.', alpha=.4)
    cb = fig.colorbar(colmap)

    iss = interact_string.split('_')

    ax.set_xlabel(iss[0])
    ax.set_ylabel(iss[1])
    ax.set_zlabel(iss[2])


    plt.savefig('{0}_3D.png'.format(''.join(interact_string.split('.csv'))) , dpi=600  )   #subtract 2 for the two coordinate columns
    plt.clf()
    plt.close(fig)
