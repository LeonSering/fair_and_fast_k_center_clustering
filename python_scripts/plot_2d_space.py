#!/usr/bin/python3

import sys
import matplotlib.pyplot as plt
from matplotlib import colors

colormap = ['k','r','g','c','m','y','b',"orange","lime","pink"]


x = []
y = []
c = []

with open(sys.argv[1],'r') as reader:
    lines = reader.read().splitlines()
    for line in lines:
        line_content = line.split(',')
        x.append(float(line_content[0]))
        y.append(float(line_content[1]))
        c.append(int(line_content[2]))
colors = [colormap[i] for i in c]
print("points: ", list(zip(x, y)))
print("colors: ", c)

fig = plt.figure()

plt.scatter(x, y, c = colors)
plt.axis('equal')

if len(sys.argv) > 2: # if there is a third parameters, plot centers
    center_indices = []
    with open(sys.argv[2],'r') as reader:
        line = reader.readline()
        content = line.split(',')
        for c in content:
            center_indices.append(int(c))
        print(center_indices)
    x_center = [x[j] for j in center_indices]
    y_center = [y[j] for j in center_indices]
    colors_center = [colors[j] for j in center_indices]
    plt.scatter(x_center,y_center, c = colors_center, marker = 'x', s = 300)   
plt.show()
