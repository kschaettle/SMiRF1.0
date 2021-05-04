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



interact_names = ['_'.join(entry.split('_')[:3]) for entry in os.listdir('.') if len(entry.split('_')) == 4 and '_HAL.csv' in entry]
interact_names.sort()
print(len(interact_names))

print_counter = 0
for interact_string in interact_names:
    inname = interact_string + '_HAL.csv'
    infile = open(inname, 'r')
    inlines = infile.readlines()
    infile.close()
    header = ['Index'] + interact_string.split('_') + ['Yield']
    splitlines = [[float(''.join(subentry.split('\"'))) for subentry in entry[:-1].split(',')] for entry in inlines[1:]] 

    print(interact_string)
    iss = interact_string.split('_')
    print_counter += 1
    print(str(print_counter) + ' out of ' + str(len(interact_names)))

    keep_indices = [header.index(entry) for entry in iss]  #keep these indices for the grid
    keep_indices.append(header.index('Yield'))   #also keep the yield

    grids = []
    grids = [[line[entry] for entry in keep_indices] for line in splitlines]

    sub_infiles = [interact_string]

    grids = np.array(grids)
    print(grids.shape)
#    deviations = np.mean(grids, axis=0)
    deviations = grids # for deviations

    for dim_index in range(3):
        print(iss[(dim_index+2)%3])
        #one 2d plot per dimension
        #first we bin z data and choose one point from each grid
        num_grids = int(deviations.shape[0]**(1/3))

        xs = deviations[:, dim_index%3]
        ys = deviations[:, (dim_index+1)%3]
        tfd = deviations[:, (dim_index+2)%3]   #temporary; the fourth dimension

        sub_deviations = copy.deepcopy(deviations)

        sub_deviations[:,0], sub_deviations[:,1], sub_deviations[:,2] = xs, ys, tfd    #cycle these three

        sub_deviations = sub_deviations[sub_deviations[:,2].argsort()] #sort on tfd 
        zs = sub_deviations[:,-1]   #for initial color map
        minz = np.min(zs)
        maxz = np.max(zs)

        cutoffs = 4
        for cutoff_iter  in range(cutoffs):

            fig=plt.figure()
            ax=fig.gca(projection='3d')

            total_points = sub_deviations.shape[0]   #number of points
            xs = sub_deviations[int((cutoff_iter/cutoffs) * total_points):int(((cutoff_iter+1)/cutoffs)*total_points) ,  0]
            ys = sub_deviations[int((cutoff_iter/cutoffs) * total_points):int(((cutoff_iter+1)/cutoffs)*total_points) ,  1]
            zs = sub_deviations[int((cutoff_iter/cutoffs) * total_points):int(((cutoff_iter+1)/cutoffs)*total_points) , -1]

            ##########This section just for twocolor plots
            colmin = np.array([0,0,1])
            colmax = np.array([1,0,0])
            white = np.array([1,1,1])
            mylambda = cutoff_iter/(cutoffs-1)
            if mylambda <= 0.5:
                mycolor = colmin + (white - colmin) * (1 - 2 * max(0, 0.5-mylambda))
            else:
                mycolor = colmax + (white - colmax) * (1-  2 * max(0, mylambda-0.5))
            mycolor = tuple(mycolor.tolist())
            ##############This section just for twocolor plots

            #now we average over the points
            grid_dict = {}
            for x in range(len(xs.tolist())):
                point = (xs[x], ys[x])
                if not point in grid_dict:
                    grid_dict[point] = []
                grid_dict[point].append(zs[x])
            for key in grid_dict:
                grid_dict[key] = np.average(np.array(grid_dict[key]))

            new_xs = []
            new_ys = []
            new_zs = []
            for key in grid_dict:
                new_xs.append(key[0])
                new_ys.append(key[1])
                new_zs.append(grid_dict[key])
            new_xs, new_ys, new_zs = np.array(new_xs), np.array(new_ys), np.array(new_zs)

            ax.plot_trisurf(new_xs, new_ys, new_zs,color=mycolor, shade='True' )#colors)
            #plt.colorbar(colmap1)

            ax.set_xlabel(iss[dim_index%3])
            ax.set_ylabel(iss[(dim_index+1)%3])
            ax.set_zlabel('Yield (bushel/acre)')

            plt.tight_layout()

            plt.savefig('{0}_2D_twocolor_{1}_{2}.png'.format(interact_string, iss[(dim_index+2)%3], cutoff_iter+1) , dpi=600  )   #subtract 2 for the two coordinate columns
            plt.clf()
            plt.close(fig)
