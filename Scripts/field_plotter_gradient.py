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
import matplotlib.mlab as ml
import matplotlib.tri as tri
from scipy.ndimage import gaussian_filter

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
plot_col_names = ['Phosphorus', 'Potassium'] #['Potassium', 'Phosphorus']

orig_name = sys.argv[1]            #don't need this if not using scaling
gradient_name = sys.argv[2]
orig_scaled_name = orig_name


header, orig_lines = get_splitlines(orig_name)
print(header)
grad_lines = get_splitlines(gradient_name)
remove_cols = [header.index(entry) for entry in remove_names]     #indices of columns to be removed
header = [header[x] for x in range(len(header)) if not x in remove_cols]
plot_cols = [header.index(entry) for entry in plot_col_names]
print(header)
print(plot_cols)
orig_lines = [[entry[x] for x in range(len(entry)) if not x in remove_cols] for entry in orig_lines]  #takes out these columns... now directly comparable
os_lines = copy.deepcopy(orig_lines)


os_array   = np.array(os_lines)
grad_array = np.array(grad_lines)
orig_array = np.array(orig_lines)    #now we can just use mins and maxs to rescale them

arrays=[os_array, grad_array, orig_array]
new_arrays = []
for array in arrays:
    newcol = np.max(array[:,:2]) * array[:,1] + array[:,0]
    newcol = np.reshape(newcol, (newcol.size, 1))
    new_array = np.concatenate((array, newcol), axis=1)
    new_array = new_array[new_array[:,-1].argsort()]
    new_arrays.append(new_array)

os_array, grad_array, orig_array = new_arrays[0], new_arrays[1], new_arrays[2]
header.append('LatLonVal')

print(os_array[0])
print(grad_array[0])
print(orig_array[0])

delta_array = grad_array - os_array
true_os_array    = os_array
true_grad_array  = grad_array    #these are both now true-valued arrays like orig_array
true_delta_array = grad_array - os_array
#prop_delta_array = true_delta_array / true_os_array        #proportion changes; this will include locations but that is not important. will be nan and inf values here.
#prop_delta_array = 100*prop_delta_array   #scaled to percentage
print(true_delta_array[0,:])

true_delta_array = np.concatenate((true_os_array[:,:2], true_delta_array[:,2:]), axis=1)   #preserve the lat lon info
true_delta_array_copy = copy.deepcopy(true_delta_array)

maxrange = 100
cutoff_fraction = 0.005     #cut off outlier values
for p_c in plot_cols:
#each plot ccolumn index
    true_delta_array = copy.deepcopy(true_delta_array_copy)
    true_delta_array = true_delta_array[true_delta_array[:, p_c].argsort()]
    print(true_delta_array[0,:])
    if cutoff_fraction:
        true_delta_array = true_delta_array[int(cutoff_fraction*true_delta_array.shape[0]):int((1-cutoff_fraction)*true_delta_array.shape[0]) , :]

    print(true_delta_array.shape)
    print(np.average(true_delta_array[:,1]))
    inline_indices = [x for x in range(true_delta_array.shape[0]) if float(true_delta_array[x,1]) < 3808500]   #180 650
    print(len(inline_indices))

    scale_factor = 1.0
    min_x = min(np.array([true_delta_array[x,0] for x in inline_indices]))
    interpolate = True

    xs =   np.array([ min_x +  (true_delta_array[x,0]-min_x)*scale_factor  for x in inline_indices])   #scales the longitude
    ys =   np.array([true_delta_array[x,1] for x in inline_indices])
    vals = np.array([true_delta_array[x, p_c] for x in inline_indices])

    gridsize = 0.25
    x_lin = np.linspace(np.min(xs), np.max(xs), int((np.max(xs)-np.min(xs))/gridsize))
    y_lin = np.linspace(np.min(ys), np.max(ys), int((np.max(ys)-np.min(ys))/gridsize))

    triang=tri.Triangulation(xs, ys)
    interpolator = tri.LinearTriInterpolator(triang, vals)
    x_mesh, y_mesh = np.meshgrid(x_lin, y_lin)
    z_mesh = interpolator(x_mesh, y_mesh)

    #smooth out extremal points
    z_mesh = gaussian_filter(z_mesh, 1)

    plt.contour(x_lin, y_lin, z_mesh, cmap='coolwarm', vmin=-30, vmax=30)
    plt.contourf(x_lin, y_lin,z_mesh, cmap='coolwarm', vmin=-30, vmax=30)



    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_scientific(False)  #gets rid of scientific notation
    plt.ticklabel_format(useOffset=False)
    plt.xticks([621800,622000])
    #plt.xticks([621800, 622000, 622200])
    #plt.yticks([3808800, 3809000])
    plt.yticks([3808000, 3808200])
    plt.tight_layout()
    plt.savefig('{0}_gradient_south.png'.format(header[p_c]), dpi=600)      #subtract 2 for the two coordinate columns
    print('success')
    plt.clf()

    continue


