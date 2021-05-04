import os
import sys
import copy
import numpy as np

inname = sys.argv[1]
infile = open(inname, 'r')
inlines = infile.readlines()
infile.close()

header = inlines[0]
inlines = inlines[1:]
header = header[:-1].split(',')

splitlines = [entry[:-1].split(',') for entry in inlines]
splitlines = [[float(subentry) for subentry in entry] for entry in splitlines]
splitlines = np.array(splitlines)

bins = 30
colnames = ['Slope','Conductivity','Zinc','Sulfur','Boron','Magnesium', 'Manganese','Copper','CEC','OrganicMatter', 'pH', 'Potassium', 'Phosphorus']
for colname in colnames:
    print(colname)
    ###here we bin the values
    bin_dict = {}
    minval = np.min(splitlines[:, header.index(colname)])
    maxval = np.max(splitlines[:, header.index(colname)])
    delta = (maxval-minval)/bins   #will need to put last val properly
    for row_iter in range(splitlines.shape[0]):
        newval = splitlines[row_iter, header.index(colname)]
        mybin = min(bins-1, int((newval-minval)/delta))
        if not mybin in bin_dict:
            bin_dict[mybin] = []
        bin_dict[mybin].append(splitlines[row_iter, header.index('Yield')])
    for key in bin_dict:
        bin_dict[key] = np.average(np.array(bin_dict[key]))   #point to average

    mysum = 0
    num_keys = 0
    for key in bin_dict:
        mysum += bin_dict[key]
        num_keys += 1
    my_avg = mysum/num_keys

    outfile = open(colname + '_1D_observed_Prevalences.csv', 'w')
    for x in range(bins):
        newstring = str(x) + ','
        if x in bin_dict:
            newstring += str(bin_dict[x]) + '\n'
            outfile.write(newstring)
        else:
            newstring += str(my_avg) + '\n'
            outfile.write(newstring)
    outfile.close()
