import os
import sys
import copy

innames = [f for f in os.listdir('.') if '_pointshal.csv' in f and len(f.split('_')) < 5 ]

file_counter = 0
for inname in innames:
    print(inname)
    file_counter += 1
    infile = open(inname, 'r')
    inlines = infile.readlines()
    infile.close()

    splitlines = [''.join(line.split('\"')) for line in inlines]
    splitlines = [line[:-1].split(',') for line in splitlines]
    splitlines = [entry[1:] for entry in splitlines]
    header, splitlines = splitlines[0], splitlines[1:]

    if file_counter == 1:
        outheader = header[:-1]
        outsplit  = [entry[:-1] for entry in splitlines]
    outheader.append(inname.split('_pointshal.csv')[0])
    outsplit = [outsplit[x] + [splitlines[x][-1]] for x in range(len(splitlines))]

outfile = open('zzzevaluated_lines_Prevalences_HAL.csv', 'w')
outfile.write(','.join(outheader) + '\n')
for entry in outsplit:
    outfile.write(','.join(entry) + '\n')
