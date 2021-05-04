### This file will take a csv formatted file with [lon][lat]...[plotted value] and produce a Matlab plot
import os
import sys
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import copy

#import plotly.plotly as py
#import plotly.graph_objs as go
#import pandas as pd



#column_plot_array = [x for x in range(2,20)]


matplotlib.use('agg')
plt.switch_backend('agg')
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

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

def get_sub_infiles(my_set, cutoff, string):
#returns a subset of the infiles with a cutoff, string
    outnames = []
    for sub_string in my_set:
        if string in sub_string and len(string.split('_')) == (len(sub_string.split('_')) -1):
            subsub = sub_string[:-4]    #cuts out extension
            minus_string = subsub.split('-')[-1]
            plus_string = subsub.split('+')[-1]
            if len(plus_string) < len(minus_string):
                num_string = plus_string
            else:
                num_string = minus_string    #gets the number...
            number = int(num_string)
            if number >= cutoff:
                outnames.append(sub_string)
    return outnames


#extension = 'zzzasurface_backup'
extension = ''
interact_string = sys.argv[1]
infiles = [f for f in os.listdir(os.path.join('.', extension)) if '_HAL.csv' in f and not 'Prevalences' in f and len(f.split('_')) == 3] #and interact_string in f]
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

    xs = np.arange(minx, maxx, (maxx - minx)/(int(len(deviations[:,0])**0.5)-1))
    ys = np.arange(miny, maxy, (maxy - miny)/(int(len(deviations[:,0])**0.5)-1))

    value_dict = {}
    for x in range(len(deviations[:,0])):
        value_dict[(round(deviations[x,0],3), round(deviations[x,1],3))] = deviations[x,-1]  #map to yield

    xs, ys = np.meshgrid(xs, ys)
    new_deviations = np.zeros(xs.shape)
    for x in range(len(xs[:,0])):
        for y in range(len(xs[1,:])):
            new_deviations[x,y] = value_dict[round(xs[x,y],3), round(ys[x,y],3)]  #put in new values

    fig=plt.figure()
    ax=fig.gca(projection='3d')

    contour_shift = [0, 0, 5]
    contour_shift = [0,0,0]
    ax.plot_surface(xs, ys, new_deviations, alpha=1.0, cmap='coolwarm', rstride=1, cstride=1)   #vmax = yield_max, vmin = yield_min
#    cset=ax.contourf(xs, ys, deviations, zdir='z', offset=yield_min-contour_shift[2]   , cmap=cm.coolwarm)
#            cset=ax.contourf(xi, yi, R, zdir='x', offset=min_x-contour_shift[0], cmap=cm.coolwarm)
 #           cset=ax.contourf(xi, yi, R, zdir='y', offset=max_y+contour_shift[1], cmap=cm.coolwarm)

    both_names = interact_string.split('_')   #array with both names  
       
    plt.xlabel('{0}'.format(both_names[0]))
    plt.ylabel('{0}'.format(both_names[1]))

    ax.set_xlim(np.min(xs), np.max(xs))
    ax.set_ylim(np.min(ys), np.max(ys))
#    ax.set_zlabel('Predicted yield (bushel/acre)')   #Local stdev. in 
    ax.set_zlabel('Yield (bushel/acre)')
    ax.set_zlim(yield_min-contour_shift[2], yield_max)

    ax.get_yaxis().get_major_formatter().set_scientific(False) 	#gets rid of scientific notation
    plt.ticklabel_format(useOffset=False)

    for view_iter in range(4):
        ax.view_init(30, 60+view_iter*90)
        plt.savefig('{0}_{1}.png'.format(''.join(interact_string.split('.csv')), view_iter) , dpi=600)	#subtract 2 for the two coordinate columns
    plt.clf()
    plt.close(fig)

