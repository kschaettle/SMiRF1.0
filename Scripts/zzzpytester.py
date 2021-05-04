import numpy as np
import time
import random
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator as RGI

dimension = 4
cutoff = 100000

points = np.array([[random.random() for x in range(dimension+1)] for y in range(cutoff)])
points2 =np.array([[random.random() for x in range(dimension  )] for y in range(cutoff)])


time1 = time.time()
outvals = griddata(points[:, :dimension], points[:, -1], points2, method='nearest')
time2 = time.time()
print(time2-time1)


grids = [30 for dim_iter in range(dimension)]
####Now we test grid linear interpolant
linspaces = [np.linspace(0,10,grids[x]) for x in range(dimension)]
values = np.array([random.random() for x in range(int(np.prod(np.array(grids)))) ])
values = np.reshape(values, tuple(grids))
print(values.shape)

time1 = time.time()
my_interpolant = RGI(tuple(linspaces), values, method='linear')
time2 = time.time()
print(time2-time1)

time3 = time.time()
outvals = my_interpolant(points2)
time4 = time.time()
print(time4-time3)
print(outvals.shape)
