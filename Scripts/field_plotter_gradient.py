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
import copy

matplotlib.use('agg')
plt.switch_backend('agg')
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

def get_splitlines(name):
#returns the header and splitlines OR just the splitlines
    infile = open(name, 'r')
    inlines = infile.readlines()
    infile.close()
    header = False
    if inlines[0][0] == '#' or inlines[0][0] == '\ufeff':
        header = True
    
    splitlines = [''.join(entry.split('#')) for entry in inlines]
    splitlines = [entry[:-1].split(',') for entry in splitlines]
    if header:
        header = splitlines[0]
        splitlines = splitlines[1:]
    float_splitlines = [[float(subentry) for subentry in entry] for entry in splitlines]
    if header:
        return header, float_splitlines
    return float_splitlines


remove_names = ['Lime', 'Potash', 'TSP', 'NIR','RedEdge','Red','Green','Blue','NDVI', 'Yield2018', 'Yield'] # 'Conductivity', 'Slope', 'pH', 'Zinc', 'Sulfur', 'Boron', 'Magnesium', 'Manganese', 'Copper', 'CEC', 'OrganicMatter'] 
#remove_cols =  [3,4,5,6,7,8, 25]   # to be removed from original datafile
plot_col_names = ['Potassium', 'Phosphorus']

orig_name = sys.argv[1]
orig_scaled_name = sys.argv[2]
gradient_name = sys.argv[3]

header, orig_lines = get_splitlines(orig_name)
print(header)
os_lines = get_splitlines(orig_scaled_name)
grad_lines = get_splitlines(gradient_name)
remove_cols = [header.index(entry) for entry in remove_names]     #indices of columns to be removed
header = [header[x] for x in range(len(header)) if not x in remove_cols]
plot_cols = [header.index(entry) for entry in plot_col_names]
print(header)
print(plot_cols)
orig_lines = [[entry[x] for x in range(len(entry)) if not x in remove_cols] for entry in orig_lines]  #takes out these columns... now directly comparable

os_array   = np.array(os_lines)
grad_array = np.array(grad_lines)
orig_array = np.array(orig_lines)    #now we can just use mins and maxs to rescale them

print(os_array[0])
print(grad_array[0])
print(orig_array[0])

delta_array= grad_array - os_array   #this is the change scaled from 0 to 1
orig_maxs  = np.max(orig_array, axis=0)
orig_mins  = np.min(orig_array, axis=0)
orig_deltas = orig_maxs-orig_mins         #full scale deltas deltas

print(orig_deltas)

min_array = np.expand_dims(orig_mins, axis=1)   #just add a dimensions
print(min_array.shape)
min_array = np.repeat(min_array, delta_array.shape[0] , axis=1)
min_array = min_array.transpose()

true_os_array    = min_array + orig_deltas * os_array
true_grad_array  = min_array + orig_deltas * grad_array    #these are both now true-valued arrays like orig_array
true_delta_array = delta_array * orig_deltas
prop_delta_array = true_delta_array / true_os_array        #proportion changes; this will include locations but that is not important. will be nan and inf values here.
prop_delta_array = 100*prop_delta_array   #scaled to percentage

prop_delta_array = prop_delta_array.tolist()
true_os_array    = true_os_array.tolist()
maxrange = 100
for p_c in plot_cols:
#each plot ccolumn index
    inline_indices = [x for x in range(len(true_os_array)) if true_os_array[x][1] > 3808650]   #180 650

    scale_factor = 1.0
    min_x = min(np.array([true_os_array[x][0] for x in inline_indices]))
    interpolate = True

    xs =   np.array([ min_x +  (true_os_array[x][0]-min_x)*scale_factor  for x in inline_indices])
    ys =   np.array([true_os_array[x][1] for x in inline_indices])
    vals = np.array([prop_delta_array[x][p_c] for x in inline_indices])

    gridsize = 0.25
    max_x, min_x = max(xs), min(xs)
    max_y, min_y = max(ys), min(ys)

    grid_x = np.linspace(min_x, max_x, int((max_x - min_x)/gridsize)+1 )
    grid_y = np.flip(np.linspace(min_y, max_y, int((max_y - min_y)/gridsize ) +1 ) , 0)

    meshgrid = np.meshgrid(grid_x, grid_y)
    xi, yi = meshgrid
    meshgrid = np.array([[xi[x][y], yi[x][y]]  for x in range(len(xi)) for y in range(len(xi[0])) ])
    print(xi.shape)
    print(yi.shape)

    xys = np.array([ true_os_array[x][:2] for x in inline_indices])

    R = griddata(xys, vals, meshgrid, 'linear')  #cubic
    R = np.array(R)
    R = np.array(R).reshape(len(grid_y), len(grid_x))

    print(R.shape)
    #plt.imshow(R, extent=[min_x, max_x, min_y, max_y], vmin=-1*maxrange, vmax=maxrange , cmap= 'seismic')   #'Reds' )

    plt.contourf(xi, yi, R, vmin=0, vmax=maxrange, cmap='coolwarm')
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_scientific(False)  #gets rid of scientific notation
    plt.ticklabel_format(useOffset=False)
    #plt.xticks([621800,622000])
    plt.xticks([621800, 622000, 622200])
    plt.yticks([3808800, 3809000])
    #plt.yticks([3808000, 3808200])

    #outdex = sys.argv[2]
    #plt.savefig('UAV_extrapolated_north_{0}.tiff'.format(outdex))
    plt.tight_layout()
    plt.savefig('{0}_gradient_north.png'.format(header[p_c]), dpi=600)      #subtract 2 for the two coordinate columns
    print('success')
    plt.clf()


