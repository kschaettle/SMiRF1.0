import os
import sys

filenames = []
counter = 0
while True:
    counter += 1
    try:
        inname = sys.argv[counter]
        filenames.append(inname)
    except:
        break
print(filenames)

running_lines = []
for filename in filenames:
    file = open(filename, 'r')
    file_lines = file.readlines()
    file.close()

    if not running_lines:
        running_lines = file_lines[:]
        running_lines = [line[:-1].split(',') for line in running_lines]
    else:
        running_lines = [running_lines[x] + file_lines[x][:-1].split(',') for x in range(len(running_lines))]
print(len(running_lines[0]))

#get first instance of each column_entry
unique_cols = list(set(running_lines[0]))  #then get first one
min_dict = {}
for entry in unique_cols:
    min_dict[entry] = len(running_lines[0])
for x in range(len(running_lines[0])):
    entry = running_lines[0][x]
    min_dict[entry] = min(min_dict[entry], x)   #find minimum

retain_cols = []
for key in min_dict:
    retain_cols.append(min_dict[key])
retain_cols.sort()   #retain the columns we want

print(len(retain_cols))
running_lines = [[line[entry] for entry in retain_cols] for line in running_lines]
outfile = open('zzzevaluated_lines_Prevalences_combined_HAL.csv', 'w')
for entry in running_lines:
    outfile.write(','.join(entry) + '\n')
outfile.close()
