### This file will take a csv formatted file with [lon][lat]...[plotted value] and produce a Matlab plot
### Used for paper images of field residuals
### Also used for making the local importance maps; simply remove the holdout limits and change file extension
import os
import sys
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import griddata
from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize









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
  
    
    splitlines = [' '.join(line.split()) for line in lines]
    splitlines = [line.split(',') for line in splitlines]
    splitlines = [[try_float(entry) for entry in line] for line in splitlines]
    splitlines = [[entry for entry in line if not entry == None] for line in splitlines]
    splitlines = [line for line in splitlines if not line == []]
    return splitlines


def abs_process(splitlines):
#absolute value processing
    splitlines = [entry for entry in splitlines if len(entry) >= 3]
    max_entry_length = max([len(entry) for entry in splitlines])
    if max_entry_length == 3:
        return splitlines

    total_influence_score = [sum([abs(subentry) for subentry in entry[2:]] ) for entry in splitlines]
#    splitlines = [splitlines[x][:2] + [splitlines[x][column_plot_index]*1/total_influence_score[x] ] for x in range(len(splitlines))] #normalized value
    splitlines = [splitlines[x][:2] + [splitlines[x][column_plot_index]] for x in range(len(splitlines))]
    return splitlines


def stat_print(split_inlines, index):
#prints out the proportion of points where the indexed variable is most important
    abslines = [[abs(subentry) for subentry in entry[2:]] for entry in split_inlines]
    maxcounter = 0
    for entry in abslines:    
        if max(entry) == entry[index-2]:
            maxcounter += 1
    print(float(maxcounter)/float(len(abslines)))





inname = sys.argv[1]
infile = open(inname, 'r')
inlines =  infile.readlines()
infile.close()
split_inlines = split_process(inlines)
column_plot_array = [x for x in range(2, len(split_inlines[0]))]
#stat_print(split_inlines)


sublines = [entry[2:] for entry in split_inlines]
maxrange = np.max(abs(np.array(sublines)))    #takes the max value
###just set the maximum to be five for consistency across models...
maxrange = 5
holdout_limits = [621880, 621920] 

for column_plot_index in column_plot_array:
    print(column_plot_index-2)
    #stat_print(split_inlines, column_plot_index)
    inlines = abs_process(split_inlines)

    print(inlines[:10])

    inlines = [entry for entry in inlines if len(entry) == 3]
    inlines = [entry for entry in inlines if entry[0] > holdout_limits[0] and entry[0] < holdout_limits[1]]

    interpolate = True
    xs = np.array([entry[0] for entry in inlines])
    ys = np.array([entry[1] for entry in inlines])
    vals = np.array([entry[2] for entry in inlines])

    gridsize = 0.25
    max_x, min_x = max(xs), min(xs)
    max_y, min_y = max(ys), min(ys)
    

    grid_x = np.linspace(min_x, max_x, int((max_x - min_x)/gridsize)+1 )
    grid_y = np.flip(np.linspace(min_y, max_y, int((max_y - min_y)/gridsize ) +1 ) , 0)


    meshgrid = np.meshgrid(grid_x, grid_y)
    xi, yi = meshgrid
    meshgrid = np.array([[xi[x][y], yi[x][y]]  for x in range(len(xi)) for y in range(len(xi[0])) ])



    xys = np.array([entry[:2] for entry in inlines])

    R = griddata(xys, vals, meshgrid, 'linear')  #cubic
    R = np.array(R)
    R = np.array(R).reshape(len(grid_y), len(grid_x))


    plt.imshow(R, extent=[min_x, max_x, min_y, max_y], vmin=-1*maxrange, vmax=maxrange , cmap= 'coolwarm')   #'Reds' )
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_scientific(False) 	#gets rid of scientific notation
    plt.ticklabel_format(useOffset=False)
   # plt.xticks([621800,622000])
  #  plt.xticks([621800, 622000, 622200])
   # plt.yticks([3808800, 3809000])
  #  plt.yticks([3808000, 3808200])

    #outdex = sys.argv[2]
    #plt.savefig('UAV_extrapolated_north_{0}.tiff'.format(outdex))
    plt.savefig('{1}_residualplot_{0}.png'.format(column_plot_index-2, inname), dpi=1200)	#subtract 2 for the two coordinate columns
    plt.clf()
