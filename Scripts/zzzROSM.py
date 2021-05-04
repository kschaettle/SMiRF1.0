import os
import numpy as np
import scipy
import sys
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn import linear_model
from sklearn.linear_model import Lasso
import statsmodels.api as sm
import pandas
import math as m
import pandas as pd
import random

from scipy.interpolate import griddata
from numpy import ma
import copy


want_minmax =  False   #converts percentile to minmax format for gradient climbing; should be set to FALSE for MARS (already in minmax format)
want_yieldpenalty = False  #uses arctan penalty against the desired yield min; balances ROI - yield tradeoff
want_uncertainty = False

want_covariance = False     #use covariance mat
downsampling_factor = 101    #reduce ALL mats by this factor to work with cov. array dimensions

if want_yieldpenalty:
    yield_setpoint = 71.2  #in bushels/acre
    yield_charscale = 0.0001   #characteristic scale of arctan change; should be < 1 bushel/acre

def get_covariances(percentile_mat, column_titles, cov_array, interact_list, beta_hat, cov_downsample=1):
#covariance downsample is essentially fraction of points in list to look at to reduce size
#need to evaluate covariance based on the surfaces
    beta_hat_noconst = beta_hat.tolist()[1:]
    beta_mat = np.array([[beta_hat_noconst[x] * beta_hat_noconst[y] for x in range(len(beta_hat_noconst))] for y in range(len(beta_hat_noconst))])
    #this is the scaling between the surface covariance and model covariance
    num_downsampled_points = int(percentile_mat.shape[0]/cov_downsample)
    print(num_downsampled_points)


    length_limit = cov_array.shape[-1]   #this is the length of the matrix
    
    outmat = [[0 for x in range(num_downsampled_points)] for y in range(num_downsampled_points)]
    for p1_index in range(num_downsampled_points):
        print(p1_index)
        true_p1_index = p1_index * cov_downsample
        for p2_index in range(p1_index, num_downsampled_points):
            true_p2_index = p2_index * cov_downsample
            covariance_grid = np.zeros(np.array(beta_mat).shape)   #store covariance values between surfaces
            for interact_index_1 in range(len(interact_list)):
                interact_1_name = interact_list[interact_index_1]
                interact_1_split = ''.join(''.join(interact_1_name.split('-')).split('+'))  #cuts out signs
                interact_1_split = interact_1_split.split('_')   #
                interact_1_cols = [column_titles.index(entry) for entry in interact_1_split]  #indices of cols
                interact_1_indices = [ min(int(percentile_mat[true_p1_index][entry]*length_limit), length_limit-1)   for entry in interact_1_cols]  #indices to use on surface-surface covmati
                for interact_index_2 in range(interact_index_1, len(interact_list)):
                    interact_2_name = interact_list[interact_index_2]   #names of interactions
                    interact_2_split = ''.join(''.join(interact_2_name.split('-')).split('+'))
                    interact_2_split = interact_2_split.split('_')   #
                    interact_2_cols = [column_titles.index(entry) for entry in interact_2_split]
                    interact_2_indices = [ min(int(percentile_mat[true_p2_index][entry]*length_limit), length_limit-1)   for entry in interact_2_cols]
                    covariance_grid[interact_index_1][interact_index_2] = cov_array[interact_index_1][interact_index_2][interact_1_indices[0]][interact_1_indices[1]][interact_2_indices[0]][interact_2_indices[1]]
                    covariance_grid[interact_index_2][interact_index_1] = covariance_grid[interact_index_1][interact_index_2]   #repeat by symmetry...
            total_covariance = np.sum(np.array(covariance_grid) * np.array(beta_mat))  #by linearity
            outmat[p1_index][p2_index] = total_covariance
            outmat[p2_index][p1_index] = total_covariance
    return outmat

def shrink(data, shrink_array):
#shrinks each dimension according to array 
    data = np.array(data)
    shrink_nparray = np.array(shrink_array)
    datashape = np.array(list(data.shape))    #need to subtract
    overages = np.mod(datashape, shrink_nparray)   #amount of extra entries in each dim.
    limits = (datashape - overages).tolist()
    for dimindex in range(len(data.shape)):
        data = np.take(data, [x for x in range(limits[dimindex])] , dimindex)   #truncates on each dim.

    rows, cols = shrink_array[0], shrink_array[1]

    return data.reshape(rows, int(data.shape[0]/rows), cols, int(data.shape[1]/cols)).sum(axis=1).sum(axis=2)

def reshape(array, interact_list, sample_mat):
#reshapes array to be len(interact_list)^2 x length^4 for covariance tensor
#unwrap from 2D back to 6D; this will be MUCH easier to work with in numpy later on
    num_interactions = len(interact_list)
    side_length = np.array(sample_mat).shape[0]    #assumes same in each dim.
    newmat = np.zeros((num_interactions, num_interactions, side_length, side_length, side_length, side_length))

    s_l_2 = side_length ** 2   #this will get used a few times

    for x in range(num_interactions):
        print('On interaction {0} out of {1}'.format(x+1, num_interactions))
        first_index_x = x * s_l_2    #first index of matrix
        for y in range(num_interactions):
            surface_surface_array = np.zeros((side_length, side_length, side_length, side_length))
            first_index_y = y * s_l_2   #use this + s_l**2 to get cutoffs
            submat = array[first_index_x:(first_index_x+s_l_2), first_index_y:(first_index_y+s_l_2)]          #this is the (still unwrapped) surface-surface interaction
            for x2 in range(side_length):
                x2sl = x2 * side_length        #need this for coordinate
                for y2 in range(side_length):
                #Take subpart of matrix, then upwrap AGAIN
                    surface_2_vec = submat[x2sl + y2][:]    #now we unwrap down to 2d
                    surface_2_vec = np.reshape(surface_2_vec, surface_2_vec.size)   #make sure this is row
                    inner_mat = np.zeros((side_length, side_length))
                    for inner_index in range(surface_2_vec.size):   #get innermost mat
                        inner_mat[int(inner_index/side_length)][inner_index % side_length] = surface_2_vec[inner_index]   #innermost index
                    newmat[x][y][x2][y2][:][:] = copy.deepcopy(inner_mat)

    print(newmat.shape)
    return newmat


def dumpout(array, filename='zzzDumped_array.csv'):
#creates textfile from a 2D array
    listarray = array.tolist()
    listarray = [[str(subentry) for subentry in entry] for entry in listarray]
    listarray = [','.join(entry) + '\n' for entry in listarray]
    outfile = open(filename, 'w')
    for entry in listarray:
        outfile.write(entry)
    outfile.close()

def array_print(nparray, filename):
#similar to function above; creates textfile from np array
    outfile = open(filename, 'w')
    outarray = nparray.tolist()
    for x in range(len(outarray)):
        entry = outarray[x]   #should ALSO be a list...
        entry = [str(subentry) for subentry in entry]
        outfile.write(','.join(entry)  + '\n')
    outfile.close()

def sort_fb_beta(val_array, new_name_array, old_name_array):
#sort the values between the old and new name arrays
    out_val_array = []
    for x in range(len(old_name_array)):
        old_name = old_name_array[x]
        try:
            new_index = new_name_array.index(old_name)
            value = val_array[new_index]
            out_val_array.append(value)
        except:
            out_val_array.append(0)
    return np.array(out_val_array)

def stepwise_selection(X, y, 
                       initial_list=[], 
                       threshold_in=0.001, 
                       threshold_out = 0.05, 
                       verbose=True):
    """ Perform a forward-backward feature selection 
    based on p-value from statsmodels.api.OLS
    Arguments:
        X - pandas.DataFrame with candidate features
        y - list-like with the target
        initial_list - list of features to start with (column names of X)
        threshold_in - include a feature if its p-value < threshold_in
        threshold_out - exclude a feature if its p-value > threshold_out
        verbose - whether to print the sequence of inclusions and exclusions
    Returns: list of selected features 
    Always set threshold_in < threshold_out to avoid infinite looping.
    See https://en.wikipedia.org/wiki/Stepwise_regression for the details
    """
    included = list(initial_list)
    dropped = []
    while True:
        changed=False
        # forward step
        excluded = list(set(X.columns)-set(included))
        new_pval = pd.Series(index=excluded)
        for new_column in excluded:
            model = sm.OLS(y, sm.add_constant(pd.DataFrame(X[included+[new_column]]))).fit()
            new_pval[new_column] = model.pvalues[new_column]
        best_pval = new_pval.min()
        if best_pval < threshold_in:
            best_feature = new_pval.argmin()
            if not excluded[best_feature] in dropped:
                included.append(excluded[best_feature])   #included.append(best_feature)
                changed=True
                if verbose:
                    print('Add  {:30} with p-value {:.6}'.format(excluded[best_feature], best_pval))

        # backward step
        model = sm.OLS(y, sm.add_constant(pd.DataFrame(X[included]))).fit()
        # use all coefs except intercept
        pvalues = model.pvalues.iloc[1:]
        worst_pval = pvalues.max() # null if pvalues is empty
        worst_feature = pvalues.argmax()
        if worst_pval > threshold_out and excluded[worst_feature] in included:
            changed=True
            included.remove(excluded[worst_feature])
            dropped.append(excluded[best_feature])
            if verbose:
                print('Drop {:30} with p-value {:.6}'.format(excluded[worst_feature], worst_pval))
        if not changed:
            break
    return included, model


def datprint(X, Y, colnames):
#creates a textfile from a matrix X and a vector Y
    outfile = open('zzzXARRAY.dat', 'w')
    Xlist = X.tolist()
    Ylist = Y.tolist()
    outfile.write('\t'.join(colnames) + '\t' + 'yield' + '\n')
    for x in range(len(X)):
        outstring = ''
        for y in range(len(X[0])):
            outstring = outstring + str(Xlist[x][y]) + '\t'
        outstring = outstring + str(Ylist[x]) + '\n'
        outfile.write(outstring)
    outfile.close()


def sort_print(valarray, namearray, outname):
    #sorts value array & name array, assuming they correspond, and then puts into textfile
    tuple_array = [tuple([abs(valarray[x]), namearray[x], valarray[x]]) for x in range(len(valarray))]
    tuple_array.sort()
    tuple_array = tuple_array[::-1]
    outfile = open(outname, 'w')
    for entry in tuple_array:
        outstring = str(entry[1]) + '\t' + str(entry[2]) + '\n'
        outfile.write(outstring)
    outfile.close()

def pad(data):
#interpolates over convex hull of data in a matrix despite the presence of nan... useful for spotty data
    bad_indexes = np.isnan(data)
    good_indexes = np.logical_not(bad_indexes)
    good_data = data[good_indexes]
    interpolated = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data)
    data[bad_indexes] = interpolated
    return data


def get_lines(string, sample_fraction=1.0):
#returns column titles and floated, split inlines from a csv
    infile = open(string, 'r')
    inlines = infile.readlines()
    infile.close()


    splitlines = [line[:-1].split('#')[-1] for line in inlines]   #for first line
    splitlines = [''.join(''.join(line.split('[')).split(']')) for line in splitlines]
    splitlines = [line.split(',') for line in splitlines]
    splitlines = [[subentry.split(' ')[0] for subentry in entry] for entry in splitlines]

    remove_cols = []
    for col_index in range(len(splitlines[0])):
        for row_index in range(len(splitlines)):
            if splitlines[row_index][col_index] == 'nan':
                remove_cols.append(col_index)
                break
    print(remove_cols)
    splitlines = [[entry[x] for x in range(len(splitlines[0])) if not x in remove_cols] for entry in splitlines]

    column_titles, splitlines = splitlines[0], splitlines[1:]   #separates
    splitlines = [[float(subentry) for subentry in entry] for entry in splitlines]
    random.shuffle(splitlines)
    splitlines = splitlines[:int(sample_fraction*len(splitlines))]
    return column_titles, splitlines, remove_cols   #returns these portions separately

def get_yield_vec(lines, column):
#just returns a vector of the yield values given a column ID
    return np.array([entry[column] for entry in lines])

def remove_columns(titles, lines, cols=False):
#columns to remove from the data after storing elsewhere
#should NOT be negative or otherwise outside the limits
    if not cols:
        return titles, lines   #just returns original input
    cols = [entry % len(titles) for entry in cols]   #modulo the length
    titles = [titles[x] for x in range(len(titles)) if not x in cols]
    lines = [[entry[x] for x in range(len(entry)) if not x in cols] for entry in lines]
    return titles, lines

def percentile_format(Xmat, Y):
#uses percentile format from R output for X; sorts Y on same basis
#can use MINMAX version instead of percentile version... linear interpolation between points

    xy_mat = np.zeros((Xmat.shape[0], Xmat.shape[1]+1))
    percentile_vec = np.array([x/(Xmat.shape[0]-1) for x in range(Xmat.shape[0])])
    for col in range(xy_mat.shape[1]-1 ):   #for each col except last
        xy_mat[:, col] = Xmat[:, col]
    xy_mat[:,-1] = Y[:]    #last col is Y
    for col in range(xy_mat.shape[1]-1 ):   #now we sort on these cols and then reassess vals
        xy_mat = xy_mat[xy_mat[:, col].argsort()]   #sorts on this column
        xy_mat[:, col] = percentile_vec[:]          #dumps these vals into the column...
    out_Y = xy_mat[:,-1]
    out_X = xy_mat[:, [dummyvar for dummyvar in range(xy_mat.shape[-1]-1)]]
    if not want_minmax:
        return out_Y, out_X    

    xy_mat = np.zeros((Xmat.shape[0], Xmat.shape[1]+1))   #don't sort values
    for col in range(xy_mat.shape[1]-1):
        xy_mat[:, col] = Xmat[:, col]   #dummy values 
    xy_mat[:, -1] = Y[:]    #last col is Y
    for col in range(xy_mat.shape[1]-1):
        minval = np.min(xy_mat[:, col])
        maxval = np.max(xy_mat[:, col])
        coldelta = maxval - minval
        xy_mat[:, col] = (xy_mat[:, col] - minval)/coldelta    #puts between zero and 1; can't be constant
    out_Y = xy_mat[:, -1]
    out_X = xy_mat[:, [dummyvar for dummyvar in range(xy_mat.shape[-1]-1)]]
    return out_Y, out_X     #these are unshuffled
     

def get_interpolated_mat(string, substring='average'):  
#opens the file from the string... right now version is for two features
#This function returns an interpolated matrix
    #downsmaple by global variable
    def cast_to_n_dims(array, newdims, charlen):
    #turns array of one dimension into MORE dimensions
        zeros_array = np.zeros(tuple(charlen for x in range(newdims))) #zero array of right size
        for index in range(len(array)):   #need to unwrap
            index_backup = index
            pos_tuple = []
            pos_tuple = tuple(int(index/charlen**(newdims-dim_index-1)) % charlen for dim_index in range(newdims))
            zeros_array[pos_tuple] = array[index]   #add this entry here...   ###karl
        return zeros_array.tolist()
            
            
    ###Get floated, split lines
    print(string)
    inname = [f for f in os.listdir('.') if string in f and not (('_' + string) in f) and not ((string+'_') in f )  and '.csv' in f][0]
    infile = open(inname, 'r')
    inlines = infile.readlines()
    infile.close()
    splitlines = [line[:-1].split(',') for line in inlines[1:]]
    splitlines = [[float(subentry) for subentry in entry[1:]] for entry in splitlines]  #cut out the first entry

    numpy_splitlines = [np.array(entry) for entry in splitlines] #includes the yield
    numpy_coords     = [np.array(entry[:-1]) for entry in splitlines]
    total_dims = len(string.split('_'))   #need special protocol for unwrapping
 #   char_len = int(len(splitlines)**(1/total_dims))      #lenght of each dimension...
    char_len = 20
    corner_delta = numpy_coords[-1] - numpy_coords[0]    #extremal corners
    grid_delta   = corner_delta/(char_len-1)             #distance between corners

    return_array = np.zeros(tuple([char_len]*total_dims))
    for x in range(len(numpy_splitlines)):
        return_array_coord = np.rint((numpy_coords[x] - numpy_coords[0])/grid_delta).astype(np.int)  #position coordinates
        return_array[return_array_coord] = numpy_splitlines[x][-1]

    ###DOWNsample the return array by downsampling_factor
#    return_array = shrink(return_array, [downsampling_factor, downsampling_factor])
    return_array = np.array(return_array)

    print(return_array.shape)
    print(np.average(return_array))
    return return_array


def get_interpolated_mat2(interaction_name, surface_splitlines, surface_column_titles, feature_column_titles, splitlines, substring='average'):
    #need to transform surfaces 
    print(interaction_name)
    return [line[surface_column_titles.index(interaction_name)] for line in surface_splitlines]
    return None

def get_responses(submat, interact_mat):
#interpolates between responses of interaction given an array of values (submat)
    lower_ints = [[int(subentry) for subentry in entry] for entry in submat]
    upper_ints = [[int(np.ceil(subentry))  for subentry in entry] for entry in submat]
    li_np = np.array(lower_ints)
    ui_np = np.array(upper_ints)
    lambdas = (ui_np - np.array(submat))    #this is fraction of lower to take... 
    om_lambdas = 1-lambdas                 #one minus lambda... for interpolation

    super_int_mat = [[ [lower_ints[x][y], upper_ints[x][y]] for y in range(len(lower_ints[0]))] for x in range(len(lower_ints))]     #for interpolation
    super_lambda_mat = [[ [lambdas[x][y], om_lambdas[x][y]] for y in range(len(lambdas[0]))] for x in range(len(lambdas))]

###Find corners given the dimension of the array
    out_array = np.zeros(lambdas.shape[0])   #to be added to for sum
    for corner_index in range(2**lambdas.shape[-1]):   #num of cols...
        corner_tuple = tuple(int(corner_index/2**(lambdas.shape[-1]-1-dim_index) )  % 2 for dim_index in range(lambdas.shape[-1]))   #just binarizes
      #  res_corner = np.array([interact_mat
        corner_indices = [[super_int_mat[x][y][corner_tuple[y]] for y in range(len(super_int_mat[x]))] for x in range(len(submat))] 
        corner_indices = [tuple(entry) for entry in corner_indices]
        res_corner = np.array([interact_mat[entry] for entry in corner_indices])
        lambda_list = [[super_lambda_mat[x][y][corner_tuple[y]] for y in range(len(super_lambda_mat[x]))] for x in range(len(submat))]    #list of lambdas
        lambda_list = np.array([np.prod(np.array(entry)) for entry in lambda_list])   #product for weighting
        weighted_response = res_corner * lambda_list
        out_array += weighted_response    #add weighed response using array...

   # print(out_array)

    return out_array   #linear interpolation
   
def get_row_subset(array, rows):
#index of rows to choose
    list_array = array.tolist()
    list_array = np.array([list_array[x] for x in range(len(list_array)) if x in rows])
    return list_array 

def get_colvec(array):
#turns 1D numpy array in to true colvec
#This could be done in one line but cuts down on syntax later on
    listvec = array.tolist()
    return np.array([[entry] for entry in listvec])

def get_overall_derivative(beta, pd_name, varnames, interaction_mat_dict, X_percentile, column_names_signed, uncertainty_switch=False):
    #beta is the set of global coefficients, pd_name is the name of the variable
    #vector is the vector of values where we are evaluating the derivative
    #this is one of the main functions for calculating the gradient
    #vector should be scaled from 0 to 1...
    pd_index = varnames.index(pd_name)    #index where we evaluate the derivative
    outsum = 0

    uncertainty_power = 1
    if uncertainty_switch:
        uncertainty_power = 2    #for form of variance calculation

###Collect the partial derivative matrices and multiply by the global coefficients
    for key in interaction_mat_dict.keys():
        if pd_name in key:      #then it will also be in the lower directories
            #print(key)
            deriv_mats = interaction_mat_dict[key][-1]
            if pd_name in deriv_mats.keys():
                deriv_mat = deriv_mats[pd_name] 
                #now get the value of the derivative..
                interaction_name = key
                interaction_name_split = [''.join(''.join(entry.split('-')).split('+')) for entry in interaction_name.split('_')]
                colnames = interaction_name_split      #names of columns that belong in interaction
                col_indices = [varnames.index(entry) for entry in colnames]    #col indices for each variable
                submatrix = X_percentile[:, col_indices]   #list form
                submatrix = submatrix * (interaction_mat.shape[0]-1)
                responses = get_responses(submatrix, deriv_mat)     #derivative mat...
                outsum += responses * beta[column_names_signed.index(interaction_name)]**uncertainty_power  #multiplied by linear factor for interaction
                #go over the interactions
    return outsum


def evaluate_response(interaction_lines, X_array, Y_vec, column_names_signed, interaction_mat_dict, column_titles):
#record the response of the ROSM over a percentile matrix
    X_signed = np.ones((X_array.shape[0],1+len(interaction_lines) ))
    column_names_signed = ['CONSTANT'] + interaction_lines     #this will be the order we use...

    Y_vec, X_percentile =  Y_vec, X_array    #re-assigns names
    Y_average = np.mean(Y_vec)

    for x in range(1, X_signed.shape[1]):
        interaction_name = column_names_signed[x]
        interaction_mat = interaction_mat_dict[interaction_name][0]

        interaction_name_split = [''.join(''.join(entry.split('-')).split('+')) for entry in interaction_name.split('_')]  #take out sign
        colnames = interaction_name_split      #names of columns that belong in interaction
        col_indices = [column_titles.index(entry) for entry in colnames]    #col indices
        submatrix = X_percentile[:, col_indices]   #list form
        submatrix = submatrix * (interaction_mat.shape[0]-1)
        
        responses = get_responses(submatrix, interaction_mat)
        X_signed[:, x] = responses[:]   #adds the response...
    return X_signed         


def get_minmax_version(interaction_mat, interaction_name, X_array, column_titles):
#returns a minmax version of the array instead of percentile-based version 

    def get_percentile_from_minmax(minmax_val, X_array, col_index):
    #get percentile... of each entry based on the minmax value
        colvals = X_array
        true_val = colvals[0] + (colvals[-1] - colvals[0]) * minmax_val    #convert from 0 to 1 back to ordinary units
        
        index_below = 0
        while colvals[index_below] < true_val:
            index_below = index_below + 1
        if colvals[index_below] == true_val:
            return index_below/(colvals.shape[0] - 1)    #goes between 0 and 1 (if exact)

        index_below = index_below - 1
        overage = true_val - colvals[index_below]
        delta = colvals[index_below+1] - colvals[index_below]
        if delta == 0:
            return index_below

        return  (index_below + overage/delta)/(colvals.shape[0] - 1)    #between 0 and 1...
        
    #remove signs from names
    interaction_name_unsigned = ''.join(interaction_name.split('+'))
    interaction_name_unsigned = ''.join(interaction_name_unsigned.split('-'))
    interaction_name_split = interaction_name_unsigned.split('_')     #split with no signs
    col_indices = [column_titles.index(entry) for entry in interaction_name_split]  #indices of cols

    #need to convert minmax limits to a percentile format first
    col_mins  = [np.min(X_array[entry]) for entry in col_indices]
    col_maxes = [np.max(X_array[entry]) for entry in col_indices]
    col_deltas = (np.array(col_maxes) - np.array(col_mins)).tolist()    #

    numdims = len(col_mins)   #number of dimensions
    dim_size = interaction_mat.shape[0]   #number of entries (observations)
    out_mat = np.zeros(interaction_mat.shape)    #copy but with zeros

    #sort the X_array columns in the copy X_array_sorted2
    X_array_sorted = [X_array[:, col_indices[x]].tolist() for x in range(numdims)]
    X_array_sorted2 = []
    for entry in X_array_sorted:
        entry2 = copy.deepcopy(entry)
        entry2.sort()
        X_array_sorted2.append(np.array(entry2))
    X_array_sorted = X_array_sorted2


    for num_index in range(out_mat.size):    #total entries in the numpy mat...  
        minmax_tuple = []
        num_index_copy = num_index    #backup
        for dim_index in range(numdims):    #go over the dims
            new_num = (num_index_copy % (dim_size-1))/(dim_size-1)   #between 0 and 1
            num_index_copy = num_index_copy - (num_index_copy % (dim_size-1))
            num_index_copy = num_index_copy/dim_size   #(dimension size)
            minmax_tuple.append(new_num)
        minmax_tuple = tuple(minmax_tuple)   #tuple converting to minmax format; has proper dimensions for array shape
        percentile_format = tuple([get_percentile_from_minmax(minmax_tuple[x], X_array_sorted[x], col_indices[x]) for x in range(len(minmax_tuple)) ])

###evaluate the percentile formatted values
        percentile_format = list(percentile_format)
        submat = [[entry*(dim_size-1) for entry in percentile_format]]
        matrix_response = get_responses(submat, interaction_mat)  #evaluate using the old format 
        minmax_list = [int(entry*(dim_size)) for entry in minmax_tuple]
        out_mat[minmax_list] = matrix_response
    print(out_mat)
    return out_mat
         
def get_unique_sign_names(inlist):
#returns a list w/ unique elements ONLY, ignoring signs in names
    unsigned_list = [ ( ''.join( (''.join(entry.split('-'))).split('+')) , entry) for entry in inlist]
    unsigned_list.sort()
#    print(unsigned_list)
    outlist = [inlist[0]]
    for x in range(1, len(inlist)):  
        if not unsigned_list[x][0] == unsigned_list[x-1][0]:
            outlist.append(unsigned_list[x][1])    #original entry
    return outlist



###################################################FUNCTIONS##############################
holdout_limits = [621880, 621920]   #lower and upper longitude values for strip trial

remove_col_names = ['NIR','RedEdge','Red','Green','Blue','NDVI', 'Yield2018', 'Lime', 'Potash', 'TSP']  #names of cols to be removed
#remove_cols = [3,4,5,6,7,8, 25]   #removes lat, lon, and UAV columns from textfile
target_col_name = 'Yield'
#target_col = 13  #BEFORE removing remove_cols; this is the column id of the yield
sample_fraction = 0.1   #fraction of rows to keep (global)
sample_fraction2 = 0.1 #fraction of observations to keep for training

inname = sys.argv[1]
column_titles, splitlines, remove_cols = get_lines(inname, sample_fraction)
remove_cols = [column_titles.index(entry) for entry in remove_col_names] + remove_cols  #get index of these columns
target_col = column_titles.index(target_col_name)

train_rows = [x for x in range(len(splitlines)) if splitlines[x][0] < holdout_limits[0] or splitlines[x][0] > holdout_limits[1]]
train_rows = random.sample(train_rows, int(sample_fraction2 * len(train_rows)))    #subsample for speed
test_rows  = [x for x in range(len(splitlines)) if splitlines[x][0] >= holdout_limits[0] and splitlines[x][0] <= holdout_limits[1]]
print(len(train_rows))
print(len(test_rows))


###get yield vector and testing lines for model fitting
Y_vec = get_yield_vec(splitlines, target_col)   #get yield vector
yield_mean = np.mean(Y_vec)
test_lines = get_row_subset(np.array(splitlines), test_rows)
test_positions = np.array([entry[:2] for entry in test_lines.tolist()])
print(test_positions[0])

column_titles, splitlines = remove_columns(column_titles, splitlines, remove_cols + [target_col])
##now we need to split into rawdata matrix and evaluated surface matrix
feature_column_titles, surface_column_titles = [entry for entry in column_titles if not '_' in entry], column_titles[:2] + [entry for entry in column_titles if '_' in entry]    #retain the latitude and longitude
surface_splitlines = [[line[x] for x in range(len(line)) if column_titles[x] in surface_column_titles] for line in splitlines]
splitlines         = [[line[x] for x in range(len(line)) if column_titles[x] in feature_column_titles] for line in splitlines]
column_titles = feature_column_titles[:]    #only the features; 

X_array = np.array(splitlines)
surface_array = np.array(surface_splitlines)

###covariance must be pre-computed between surfaces
if want_covariance:
    cov_name = sys.argv[2]
    cov_file = open(cov_name, 'r')
    cov_lines = cov_file.readlines()
    cov_file.close()
    cov_array = [line[:-1].split(',') for line in cov_lines]
    cov_array = [[float(subentry) for subentry in entry] for entry in cov_array]
    cov_array = np.array(cov_array)    #covariance point-point array  #this needs to be reshaped



###This is all pre-processing; now we add MORE dimensions to the array for the linear nth order models
####Model 1: ALL 2nd order interactions are fitted and evaluated...
num_vars = X_array.shape[-1]    #number of variables...
ones_array = np.ones((X_array.shape[0], 1))  #constant value for constant in "linear" model
second_order_array = np.ones((X_array.shape[0], int(num_vars*(num_vars+1)/2 ) )) #second order var
second_order_titles = ['CONSTANT'] + column_titles
new_titles = []
col_counter = 0
for x in range(len(column_titles)):
    for y in range(x , len(column_titles)):
        new_title = column_titles[x] + '_' + column_titles[y]
        new_titles.append(new_title)
        second_order_array[:,col_counter] = X_array[:,x] * X_array[:,y]   #duplicates data in col
        col_counter += 1



X1_array = np.concatenate((ones_array,  np.array([line[2:] for line in X_array.tolist()])    ), axis=1)
X1_train = get_row_subset(X1_array, train_rows)
Y_train = get_row_subset(Y_vec, train_rows)
X1_test = get_row_subset(X1_array, test_rows)
Y_test = get_row_subset(Y_vec, test_rows)
beta_hat = np.linalg.lstsq(X1_train, Y_train)[0]
predicted = np.dot(X1_test, beta_hat)
residuals = Y_test - predicted
ave_residual = np.mean(np.abs(residuals))
print(X1_array.shape)
print(ave_residual)
print(beta_hat)
residuals = get_colvec(residuals)
sort_print(beta_hat, ['CONSTANT'] + column_titles[2:], 'linear_model.dat')
array_print(np.concatenate((test_positions, residuals), axis=1), 'linear_model_residuals.txt')
print('Above: linear model')



second_order_titles = second_order_titles + new_titles
X2_array = np.concatenate((ones_array, X_array, second_order_array), axis=1)

   #training and testing rows
X2_train = get_row_subset(X2_array, train_rows)
Y_train = get_row_subset(Y_vec, train_rows)
X2_test = get_row_subset(X2_array, test_rows)
Y_test = get_row_subset(Y_vec, test_rows)

beta_hat = np.linalg.lstsq(X2_train, Y_train)[0]   #train on train rows...
predicted = np.dot(X2_test, beta_hat)
residuals = Y_test - predicted
ave_residual = np.mean(np.abs(residuals))
print(X2_array.shape)
print(ave_residual)
#print(beta_hat)
residuals = get_colvec(residuals)
	#This is the second-order model.



####Model 2: similar to above, but ONLY use interactions that are IMPORTANT from RF
interaction_lines = surface_column_titles[2:]   #length of matrix will exclude the lat-lon info
X_signed = np.array([[1.0] + line[2:] for line in surface_splitlines])  #cut out the lat-lon and add a constant term
column_names_signed = ['CONSTANT'] + interaction_lines
#Y_vec, X_percentile = percentile_format(X_array, Y_vec)  #NEWEDIT

X_signed_train = get_row_subset(X_signed, train_rows)
Y_train = get_row_subset(Y_vec, train_rows)
X_signed_test = get_row_subset(X_signed, test_rows)
Y_test = get_row_subset(Y_vec, test_rows)

###Compute the Linear ROSM
beta_hat = np.linalg.lstsq(X_signed_train, Y_train)[0]
signed_predicted = np.dot(X_signed_test, beta_hat)
signed_residuals = Y_test - signed_predicted
print(np.mean(np.abs(signed_residuals) ))
print(X_signed.shape)
print(len(column_names_signed))
print(beta_hat)
sort_print(beta_hat, column_names_signed, 'linear_ROSM.dat')
signed_residuals = get_colvec(signed_residuals)
array_print(np.concatenate((test_positions ,signed_residuals), axis=1), 'Linear_important.txt')
print('Above: linear regression over interactions')


###Ridge and Lasso regression for building ROSMs
yield_mean = np.mean(Y_vec)    #mean of the yield; put in place of constant for penalized regr.

alpha_value = 0.5

clf = Ridge(alpha=alpha_value)

clf.fit(X_signed_train, Y_train)
Y_pred = clf.predict(X_signed_test)
print(np.mean(np.abs(Y_pred-Y_test)))
print(clf.coef_)
sort_print(clf.coef_, column_names_signed, 'Ridge_ROSM_{0}.dat'.format(alpha_value))
array_print(np.concatenate((test_positions, get_colvec(Y_test - Y_pred)), axis=1), 'Ridge_residuals.txt')
print('Above: Ridge regression')


clf2 = Lasso(alpha=alpha_value)
clf2.fit(X_signed_train, Y_train)
Y_pred2 = clf2.predict(X_signed_test)
print(np.mean(np.abs(Y_pred2-Y_test)))
print(clf2.coef_)
sort_print(clf2.coef_, column_names_signed, 'Lasso_ROSM_{0}.dat'.format(alpha_value))
print('Above: Lasso regression')
array_print(np.concatenate((test_positions, get_colvec(Y_test - Y_pred2)), axis=1) , 'Lasso_residuals.txt')


#X_signed_train = X_signed_train.tolist()
#X_signed_train = np.array([line[1:] for line in X_signed_train])   #for removing constant
#print(X_signed_train.shape)



X_signed_pd = pandas.DataFrame(data=X_signed_train, columns=column_names_signed) #try without constant
y_pd = pandas.DataFrame(data=Y_train, columns = ['Yield'])
datprint(X_signed, Y_vec, column_names_signed)


result, stepwise_model = stepwise_selection(X_signed_pd, y_pd)
print(stepwise_model.params.tolist())
fb_beta_hat = sort_fb_beta(stepwise_model.params.tolist(), result, column_names_signed)
predicted = np.dot(X_signed_test, fb_beta_hat)
residual = Y_test - predicted
print('forward-backward residual:')
print(np.mean(np.abs(residual)))

print('resulting features:')
result.sort()
print(result)
print(len(result))

#sort_print(stepwise_model.params.tolist(), result, 'ForwardBackward.txt')
#array_print(np.concatenate((test_positions, get_colvec(residual)), axis=1), 'stepwise_residuals.txt')


#####################################################################
Y_vec, X_percentile = percentile_format(X_array, Y_vec)
interaction_mat_dict = {}   #{name:[interaction_mat, derivs], ...}; derivs == {name:object, ...}
uncertainty_mat_dict = {}   #these dictionaries store the derivatives of the interactions
X_signed_uncertainty = copy.deepcopy(X_signed)
for x in range(1, X_signed.shape[1]):
    interaction_name = column_names_signed[x]
    interaction_mat = get_interpolated_mat(interaction_name)   #gets matrix version of the feature...

    print(np.array(interaction_mat).shape)   # test
    
    if want_minmax:
        interaction_mat = get_minmax_version(interaction_mat, interaction_name, X_array, column_titles)
        if want_uncertainty:
            uncertainty_mat = get_minmax_version(uncertainty_mat, interaction_name, X_array, column_titles)
    interaction_name_split = [''.join(''.join(entry.split('-')).split('+')) for entry in interaction_name.split('_')]  #take out sign
    colnames = interaction_name_split      #names of columns that belong in interaction
    colnames = [''.join( ''.join(subentry.split('-')).split('+')) for subentry in colnames]
    print(colnames)
    print(column_titles)
    col_indices = [column_titles.index(entry) for entry in colnames]    #col indices
    submatrix = X_percentile[:, col_indices]   #list form
    submatrix = submatrix * (interaction_mat.shape[0]-1)

    #Now we will compute and organize the derivative matrices...
    interaction_mat_dict[interaction_name] = [interaction_mat, {}]

    gradient_mat_array = np.gradient(interaction_mat)
    gradient_mat_array = [np.array(entry) for entry in gradient_mat_array]  #list of 1D gradients
    print([entry.shape for entry in gradient_mat_array])

    for y in range(len(colnames)):
        colname = colnames[y]      #name and index
        interaction_mat_dict[interaction_name][-1][colname] = gradient_mat_array[y]   #set up the structure first... 

    responses = get_responses(submatrix, interaction_mat)
    X_signed[:, x] = responses[:]   #adds the response...
    
    #we must compute the gradient of the uncertainty for each response surface as well
    if want_uncertainty:
        uncertainty_mat = get_interpolated_mat(interaction_name, 'deviations')
        uncertainty_mat = uncertainty_mat**2    #squares the read-in plot
        if want_minmax:
            uncertainty_mat = get_minmax_version(uncertainty_mat, interaction_name, X_array, column_titles) 
        print(uncertainty_mat)
        uncertainty_mat_dict[interaction_name] = [uncertainty_mat, {}]
        gradient_var_array = np.gradient(uncertainty_mat)
        gradient_var_array = [np.array(entry) for entry in gradient_var_array]
        ###store the gradients in the dictionary uncertainty_mat_dict
        for y in range(len(colnames)):
            uncertainty_mat_dict[interaction_name][-1][colname] = gradient_var_array[y]
        var_responses = get_responses(submatrix, uncertainty_mat)
        X_signed_uncertainty[:, x] = var_responses[:]   #variance responses
        X_signed_uncertainty[:, 0] = 0      #constant column uncertainty set to 0 


if want_covariance:
#now we reshape the array in this case...
    cov_array = reshape(cov_array, column_names_signed[1:], interaction_mat)  #cuts out constant


#use variance mats


beta_hat = fb_beta_hat
beta_hat_squared = (np.array(beta_hat)**2).tolist()
tolerance = 0.00000000001     #hill climbing tolerance...
integrating_factor = 100    #another hyperparameter for the inner gradient loop; optimizing ROI analogue
current_delta = tolerance * 10000
old_ROI = 0.0
portfolio_lambda = 0    #lambda for effective utility function
new_X_percentile = X_percentile.tolist()
new_X_percentile = [[min(subentry*1.0000001, 1.0) for subentry in entry] for entry in new_X_percentile]
new_X_percentile = np.array(new_X_percentile)
dumpout(new_X_percentile, 'orig_values.csv')


unadjustables = [0, 1]    #columns to not change; alternatively only optimization columns will be changed
old_model_eval = np.dot(X_signed, beta_hat)   #clf.predict(X_signed)    #this is f(v0)
print(np.mean(old_model_eval))
prices = [1/X_signed.shape[0] for entry in column_titles]  #this will have to be adjusted... price per unit per area will be important. just use 1 as price "per unit" for now
#for entry in unadjustables:
#    prices[entry] = 1


###prices must be in standardized units, for example $/lb
###Use hard-coded column ids for now
#prices[3] = 30/X_signed.shape[0]  * (np.max(X_array[3]) - np.min(X_array[3])) * 0.00001
prices[5] = 318/X_signed.shape[0] * (np.max(X_array[4]) - np.min(X_array[4])) * 0.00001
prices[4] = 425/X_signed.shape[0] * (np.max(X_array[5]) - np.min(X_array[5])) * 0.00001

###This section is for setting near global optimum phosphorus and potassium
new_X_percentile = new_X_percentile.tolist()
new_X_percentile = [entry[:4] +[0.319, 0.352] + entry[6:] for entry in new_X_percentile]
new_X_percentile = np.array(new_X_percentile)
#Lime, Potash, TSP

ones_mat =  np.ones(new_X_percentile.shape)
zeros_mat = np.zeros(new_X_percentile.shape)
optimization_columns = [4,5]     #just these columns; adjust by 1 for the constant
print(column_titles)

##establish lower and upper bounds
lower_bounds = np.array([0 for x in range(ones_mat.shape[-1])])
upper_bounds = np.array([100000000 for x in range(ones_mat.shape[-1])])   #these must be converted to 0 - 1 range format

mins = np.min(X_array, axis=0)
maxs = np.max(X_array, axis=0)
deltas = maxs - mins
lower_bounds = (lower_bounds - mins)/deltas
upper_bounds = (upper_bounds - mins)/deltas    #set to min and max
new_lower_bounds = np.zeros(ones_mat.shape)
new_upper_bounds = np.zeros(ones_mat.shape)
for column_index in range(ones_mat.shape[-1]):  #set these new matrix col vals to old ones...
    new_lower_bounds[:, column_index] = lower_bounds[column_index]
    new_upper_bounds[:, column_index] = upper_bounds[column_index] 


loop_counter = 0
while current_delta > tolerance and loop_counter < 30000:
    if want_covariance:   
#we will use this to calculate field-scale uncertainty, and other features...                                                                                       
        all_covariances = get_covariances(new_X_percentile, column_titles, cov_array, column_names_signed[1:], beta_hat, cov_downsample=50)         
        all_covariances = np.array(all_covariances)
        all_covariances = all_covariances.tolist()
        all_covariances = [[str(subentry) for subentry in entry] for entry in all_covariances]
        all_covariances = [','.join(entry)+'\n' for entry in all_covariances]



    loop_counter += 1
    partial_derivs = [get_overall_derivative(beta_hat, entry, column_titles, interaction_mat_dict, new_X_percentile, column_names_signed)
                     for entry in column_titles]   #derivative for each adjustable variable... might need to change which columns to use

    #print(partial_derivs)    

    response_vector = evaluate_response(interaction_lines, new_X_percentile, Y_vec, column_names_signed, interaction_mat_dict, column_titles) 
    new_model_eval = np.dot(response_vector, beta_hat)  #This is f(v+v0)

    if want_uncertainty:
        for x in range(1, len(column_names_signed)):
         ###get uncertainty at each location
            interaction_name = column_names_signed[x]
            uncertainty_mat = uncertainty_mat_dict[interaction_name][0]
            interaction_name_split = [entry[:-1] for entry in interaction_name.split('_')]  #
            colnames = interaction_name_split      #names of columns that belong in interaction
            col_indices = [column_titles.index(entry) for entry in colnames]    #col indices
            submatrix = new_X_percentile[:, col_indices]   #list form
            submatrix = submatrix * (uncertainty_mat.shape[0]-1) 
            X_signed_uncertainty[:, x] =  get_responses(submatrix, uncertainty_mat)
        new_model_variance = np.dot(X_signed_uncertainty, beta_hat_squared) 

    print('Current yield: '+str(np.mean(new_model_eval))+' Orig.: '+str(np.mean(old_model_eval)))
    model_delta = new_model_eval - old_model_eval  #this is mean f(v0+v) - f(v0)
    price_change = np.sum(np.dot((new_X_percentile - X_percentile), np.array(prices)))   #total price of new treatment
    print(price_change)
    pc_100000 = price_change+100000
    gradient = [0 for x in range(len(partial_derivs))]
    for x in optimization_columns:   #these are the ONLY ones we change
        if not want_yieldpenalty or np.mean(new_model_eval) > yield_setpoint:
         #prevent model from getting stuck at very low changes in yield; yield_setpoint should be slightly greater than original yield for proper performance
            gradient[x] = (partial_derivs[x]*(pc_100000) - model_delta*prices[x])/((pc_100000)**2+0.0000001) 
        else:
            gradient[x] = partial_derivs[x]*1/(abs(price_change)+0.00000001)   #just use the yield itself as the target funciton in this case... make sure derivatives are scaled accordingly

    for x in range(len(gradient)):
        if not x in optimization_columns:
            gradient[x] = np.zeros(gradient[optimization_columns[0]].shape)   #not all columns will be adjusted.
    gradient = np.array([entry.tolist() for entry in gradient])
    gradient = np.transpose(gradient)    #transposition 


    if want_uncertainty:   #portfolio optimization; first compute uncertainty, then apply lambda
        partial_derivs = [get_overall_derivative(beta_hat, entry, column_titles, uncertainty_mat_dict, new_X_percentile, column_names_signed, want_uncertainty) for entry in column_titles]
        print(partial_derivs)
        gradient2 = [0 for x in range(len(partial_derivs))]
        for x in optimization_columns:
            gradient2[x]  = partial_derivs[x]
        for x in range(len(gradient2)):
            if not x in optimization_columns:
                gradient2[x] = np.zeros(gradient2[optimization_columns[0]].shape)
        for x in range(len(gradient2)):
            if np.sum(np.abs(np.array(gradient2[x]))) == 0:
                gradient2[x] = np.zeros(gradient2[optimization_columns[0]].shape)
        gradient2 = np.array([gradient2])     #see above
        gradient2 = np.transpose(gradient2)
        #print(gradient2)
        print('Gradients1 and 2: ' + str(np.sum(gradient)) + ' ' + str(np.sum(gradient2)))
        gradient2 = np.squeeze(gradient2)
        print('Portfolio lambda: ' + str(portfolio_lambda))
        gradient = gradient - portfolio_lambda * gradient2    #climbs combined optimization function


    current_delta = np.mean(np.abs(gradient))

    dumpout(new_X_percentile)

    #update the observations using the newly computed gradient
    new_X_percentile = new_X_percentile + integrating_factor*gradient
    new_X_percentile = np.minimum(new_X_percentile, ones_mat)
    new_X_percentile = np.maximum(new_X_percentile, zeros_mat)     #keep between 0 and 1
    #    new_X_percentile = np.maximum(new_X_percentile, X_percentile)  #positive prescription
    new_X_percentile = np.minimum(new_X_percentile, new_upper_bounds)
    new_X_percentile = np.maximum(new_X_percentile, new_lower_bounds)  #upper and lower bounds
    ROI = np.mean(model_delta)/(price_change+0.00000000001)
    current_delta = abs(ROI - old_ROI)
    old_ROI = ROI

#    new_X_percentile = apply_physical_equations(new_X_percentile, X_percentile)

    print('ROI: ' + str(ROI))
    print('Avg. of gradient: ' + str(current_delta))
    print('Total price: ' + str(price_change/model_delta.shape[0]))
    if want_uncertainty:
        print('Average yield variance: {0} (bushel/acre)^2'.format(np.mean(new_model_variance)))
    print()
