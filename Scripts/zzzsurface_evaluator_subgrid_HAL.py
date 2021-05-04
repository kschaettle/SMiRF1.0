import numpy as np
import os
import sys
import math as m
import copy
import random

sigma = 1.5  #in terms of gridsize ; for gaussian kernel
sample_fraction = 0.1

def get_dataheader(datalines):
    dataheader_line = datalines[0]
    dataheader = dataheader_line[:-1].split(',')   #split on csv delimiter
    return dataheader

def cubify(max_corner, min_corner, subgridsize):
    #generates all points in a cube of arbitrary dimension given (min, min...min) and (max, max, ...,max) corners
    dimension = max_corner.size
    corner_delta = max_corner - min_corner
    points_in_each_dimension = np.ceil(np.abs((max_corner - min_corner))/subgridsize) + np.ones(max_corner.size)
    total_points = int(np.prod(points_in_each_dimension))

    cumulative_products = [1]
    for x in range(dimension-1):
        cumulative_products = [cumulative_products[0]*points_in_each_dimension[-1*x]] + cumulative_products
#cumulative products is used to take modulus of total_points to give points on a grid

    out_array = []
    for point_index in range(total_points):   #number of points in the box
        sum = point_index
        point_array = []
        for x in range(dimension):
            quotient = int(sum/cumulative_products[x])
            remainder = sum - cumulative_products[x] * quotient
            point_array.append(quotient)
            sum = remainder
        point_array = np.array(point_array)
        out_array.append(min_corner + point_array*subgridsize)
    return out_array



dataname = sys.argv[1]
start_index = int(sys.argv[2])
end_index   = int(sys.argv[3])
file_extension = sys.argv[4]
datafile = open(dataname, 'r')
datalines = datafile.readlines()
datafile.close()

outheader = [] #header needs to include rawdata info and evaluated surface values
dataheader = get_dataheader(datalines)
data_splitlines = [[float(subentry) for subentry in entry[:-1].split(',')] for entry in datalines[1:]]
outlines = copy.deepcopy(data_splitlines)

for entry in dataheader:
    outheader.append(entry)    #so each entry is in the header


surface_files = [f for f in os.listdir('.') if f[-4:] == '.csv' and not 'Prevalences' in f and 'HAL' in f]
surface_files.sort()
surface_files = surface_files[start_index:end_index]
#surface_files = [entry for entry in surface_files if len(entry.split('_')) == 4]  #for speedy testing
#print(surface_files)

surface_file_lengths = [len(entry.split('_')) for entry in surface_files]
frequencies = {}
for entry in surface_file_lengths:
    if not entry in frequencies:
        frequencies[entry] = 0
    frequencies[entry] += 1
print(frequencies)

counter = 0
for surface_filename in surface_files:
    counter += 1
    print(counter)
    print(surface_filename)
    surface_file = open(surface_filename, 'r')
    surface_lines = surface_file.readlines()
    surface_file.close()
    new_header_entry = ''.join(surface_filename[:-4].split('_HAL'))   #cuts out .csv
    outheader.append(new_header_entry)

    surface_cols = new_header_entry.split('_')   #new columns to evaluate
    data_sub_splitlines = [[line[dataheader.index(col)] for col in surface_cols] for line in data_splitlines]
       #data but only retaining the relevant columns

    dimension = len(surface_cols)
    num_grids = (len(surface_lines) - 1) ** (1/dimension)
    if num_grids == 1:
        universal_yield = float(surface_lines[-1][:-1].split(',')[-1])  #last entry is the yield
        outlines = [entry.append(universal_yield) for entry in outlines]   #all have same value

    else:
        #here we need to actually evaluate
        surface_splitlines = [[float(subentry) for subentry in line.split(',')[1:]] for line in surface_lines[1:]]
        #need to separaate the yield entry from the other entries
        surface_yield = [entry[-1] for entry in surface_splitlines]
        surface_numpylines = [np.array(entry[:-1]) for entry in surface_splitlines]
        data_sub_numpylines= [np.array(entry) for entry in data_sub_splitlines]     #these will be compared
        gridsize = (surface_numpylines[-1] - surface_numpylines[0])/(num_grids-1)  #corners of cube
    
        if dimension >= 4:
            randoms = [random.random() for x in range(len(surface_numpylines)-2)]
            surface_numpylines = surface_numpylines[:1] + [surface_numpylines[x+1] for x in range(len(randoms)) if randoms[x] < sample_fraction] + surface_numpylines[-1:]
            surface_yield      = surface_yield[:1] + [surface_yield[x+1] for x in range(len(randoms)) if randoms[x] < sample_fraction] + surface_numpylines[-1:]

        two_s_two = 2*(sigma*gridsize)**2  #different for each dimension
        subgrid_cutoff = 2*sigma*gridsize   #subgrid will only compare to points this close
        num_subgrids = (m.ceil(  (len(surface_numpylines)**0.5)**(1/dimension) ) + 1)   #number in each direction
        subgridsize  = gridsize * (num_grids - 1) / (num_subgrids - 1)                  #lenght of subgrids
        subgrids_kept = np.prod(  np.maximum(np.ceil(subgridsize/subgrid_cutoff)*2 - np.ones(dimension) , 3*np.ones(dimension)  ))                     #amount to keep for searching
        subgrid_numpylines = cubify(surface_numpylines[-1], surface_numpylines[0], subgridsize)

        map_dict = {}    #maps subgrid to regular gridpoints
        for x in range(len(surface_numpylines)):
            surface_numpyline = surface_numpylines[x]
            deltas = [(np.sum(np.abs(subgrid_numpylines[y] - surface_numpyline)**2),y) for y in range(len(subgrid_numpylines))]
            deltas.sort()
            if not deltas[0][1] in map_dict:
                map_dict[deltas[0][1]] = []
            map_dict[deltas[0][1]].append(x)    #associated with closest gridpoint
            


        for x in range(len(outlines)):
            if x % 1000 == 0:
                print(x)
            data_array = data_sub_numpylines[x]
            deltas =    [(data_array - subgrid_numpylines[y], y) for y in range(len(subgrid_numpylines))]
            distances = [(np.sum(entry[0]**2), entry[1]) for entry in deltas]   #we will sort these
            distances.sort()
            distances = distances[:int(subgrids_kept)]
            kept_subgrids = [entry[-1] for entry in distances]   #index of kept subgrids
            kept_surface_indices = []
            for key in kept_subgrids:
                if key in map_dict:
                    kept_surface_indices += map_dict[key]  #kept surface line indices
            kept_surface_numpylines = [surface_numpylines[entry] for entry in kept_surface_indices]
            kept_surface_yield = [surface_yield[entry] for entry in kept_surface_indices]

            surface_deltas = [data_array - entry for entry in kept_surface_numpylines]
            deltas2 =[entry[0]**2 for entry in surface_deltas]
            deltas2 =[-np.sum(entry/two_s_two) for entry in deltas2]
            kernel_weights = [np.exp(entry) for entry in deltas2]
            if np.sum(np.array(kernel_weights)) == 0:
                yield_val = np.average(np.array(kernel_weights))
            else:
                yield_val = np.sum(np.array(kernel_weights) * np.array(kept_surface_yield))/np.sum(np.array(kernel_weights))
            outlines[x].append(yield_val)

    outfile = open('evaluated_surfaces_Prevalences{0}.csv'.format(file_extension), 'w')
    outfile.write(','.join(outheader)+ '\n')
    outlines = [[str(subentry) for subentry in entry] for entry in outlines]
    for entry in outlines:
        outfile.write(','.join(entry) + '\n')
    outfile.close()
