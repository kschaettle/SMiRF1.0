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
evaluated_filename = sys.argv[1]
infile=open(evaluated_filename, 'r')
inlines = infile.readlines()
infile.close()

splitlines = [entry[:-1].split(',') for entry in inlines] #take out newline character
header = splitlines[0]
splitlines = splitlines[1:]
splitlines = [[','.join(subentry.split(' ')) for subentry in entry] for entry in splitlines]
new_splitlines = []
for entry in splitlines:
    new_entry = []
    for subentry in entry:
        if '[' in subentry:
            new_subentry = ''.join(''.join(subentry.split('[')).split(']')).split(',')
            new_entry.append(float(np.average(np.array([float(subsub) for subsub in new_subentry if subsub]  ))))
        else:
            new_entry.append(float(subentry))
    new_splitlines.append(new_entry)

splitlines = new_splitlines
interact_names = [entry for entry in header if len(entry.split('_')) == 3]
interact_names.sort()
print(len(interact_names))

print_counter = 0
for interact_string in interact_names:
    print(interact_string)
    iss = interact_string.split('_')
    print_counter += 1
    print(str(print_counter) + ' out of ' + str(len(interact_names)))

    keep_indices = [header.index(entry) for entry in iss]  #keep these indices for the grid
    keep_indices.append(header.index('Yield'))

    grids = []
    grids = [[line[entry] for entry in keep_indices] for line in splitlines]

    sub_infiles = [interact_string]

    grids = np.array(grids)
    print(grids.shape)
#    deviations = np.mean(grids, axis=0)
    res_up_factor = 1.0
    deviations = grids # for deviations

    iss = interact_string.split('_')  #interaction names

    for dim_index in range(3):
        #one 2d plot per dimension
        #first we bin z data and choose one point from each grid
        num_grids = 30    #number of grids in x and y dimension
        num_grids = num_grids - 1 #corrects for endpoints

        xs = deviations[:, dim_index%3]
        ys = deviations[:, (dim_index+1)%3]
        tfd = deviations[:, (dim_index+2)%3]   #temporary; the fourth dimension

        minx = np.min(xs)
        miny = np.min(ys)
        maxx = np.max(xs)
        maxy = np.max(ys)

        griddict = {}   #takes x,y grid pointers and aims at tfd values
        for iterant in range(len(xs.tolist())):
            grid = (int(num_grids*(xs[iterant]-minx)/(maxx-minx)),int(num_grids*(ys[iterant]-miny)/(maxy-miny)))
            if not grid in griddict:
                griddict[grid] = []
            griddict[grid].append((tfd[iterant], deviations[iterant,-1]))  #also append the yield

        sub_deviations = []
        for grid in griddict:
            newx = grid[0]/num_grids * (maxx-minx) + minx
            newy = grid[1]/num_grids * (maxy-miny) + miny   #reconvert to x and y
            newtfd = np.average(np.array([entry[0] for entry in griddict[grid]]))  #average over binned entries
            newyield = np.average(np.array([entry[1] for entry in griddict[grid]]))
            sub_deviations.append([newx, newy, newtfd, newyield])
        sub_deviations = np.array(sub_deviations)

        xs = sub_deviations[:, 0]
        ys = sub_deviations[:, 1]
        zs = sub_deviations[:,-1]
        the_fourth_dimension = sub_deviations[:, 2]

        fig = plt.figure(figsize=(16,8))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')

        colors = cm.coolwarm((zs-min(zs))/(max(zs)-min(zs)))
        colmap1 = cm.ScalarMappable(cmap=cm.coolwarm)
        colmap1.set_array(zs)
        colors = cm.coolwarm((the_fourth_dimension-min(the_fourth_dimension))/(max(the_fourth_dimension)-min(the_fourth_dimension)))
        colmap2 = cm.ScalarMappable(cmap=cm.coolwarm)
        colmap2.set_array(the_fourth_dimension)

        yg1 = ax1.plot_trisurf(xs, ys, zs, cmap=cm.coolwarm)#colors)
        cb  = fig.colorbar(colmap1)
        yg2 = ax2.plot_trisurf(xs, ys, the_fourth_dimension, cmap=cm.coolwarm)
        cb  = fig.colorbar(colmap2)

        ax1.set_xlabel(iss[dim_index%3])
        ax1.set_ylabel(iss[(dim_index+1)%3])
        ax1.set_zlabel('Yield')

        ax2.set_xlabel(iss[dim_index%3])
        ax2.set_ylabel(iss[(dim_index+1)%3])
        ax2.set_zlabel(iss[(dim_index+2)%3])

        plt.tight_layout()

        plt.savefig('{0}_threeaxis_{1}.png'.format(interact_string, iss[(dim_index+2)%3]) , dpi=600  )   #subtract 2 for the two coordinate columns
        plt.clf()
        plt.close(fig)
