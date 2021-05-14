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
    path = sys.argv[2]

    if(path.endswith(".centers")):
        centers = []
        with open(path,'r') as reader:
            line = reader.readline()
            content = line.split(',')
            for c in content:
                centers.append(int(c))
            print(centers)

    if(path.endswith(".clustering")):
        centers = []
        cluster = {} 
        with open(path,'r') as reader:
            lines = reader.read().splitlines()
            for line in lines:
                content = line.split(':')
                center = int(content[0])
                centers.append(center)
                points = content[1].split(',')
                print(points)
                if points[0] != '':
                    cluster[center] = [int(x) for x in points]
                else:
                    cluster[center] = []
        print(centers)
        print(cluster)
        for i in centers:
            for j in cluster[i]:
                print(i,j)
                plt.plot([x[i],x[j]], [y[i],y[j]], c = 'lightgrey', linewidth = 0.5, zorder = -1.0)

    
    x_center = [x[j] for j in centers]
    y_center = [y[j] for j in centers]
    colors_center = [colors[j] for j in centers]
    plt.scatter(x_center,y_center, c = colors_center, marker = 'x', s = 300, zorder = 1.0)
plt.show()
