#! /usr/bin/python
# -*: coding:utf-8 -*-

import scipy.interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
"""
f = open('./out_hp_2.csv')
buffers = f.readlines()
mat = []
for tmp in buffers:
	mat.append(tmp.split(','))

numat = []## must change to array
for col in mat:
	numat.append(map(lambda x: float(x), col))

## interpolation

## plotting

fig = plt.figure()
ax = Axes3D(fig)
y = [i for i in range(len(numat))]
x = [i for i in range(len(numat[0]))]
X, Y = np.meshgrid(x,y)
# 行列ではなく点の集合としてプロットするのがscatter
"""
#ax.scatter(X, Y, numat)
#plt.show()
"""
#R = np.sqrt(X**2 + Y**2)
Z = np.array(numat)
print X
print Y
print Z
"""
#data = np.loadtxt("./plot.py", delimiter=',')
"""
x = np.arange(0, len(data), 1)
y = np.arange(0, len(data[0]), 1)
X, Y = np.meshgrid(x,y)
"""
plt.plotfile('./plot.csv', (0,1), delimiter=',', names=('ppv', 'sen'))
plt.show()

