#/usr/bin/python3

import sys
import matplotlib.pyplot as plt
from matplotlib import colors

import tkinter as tk
from tkinter import filedialog


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])

class point:
    def __init__(self, x, y, c, mpl_point):
        self.x = x
        self.y = y
        self.c = c
        self.mpl_point = mpl_point

points = []

colormap = ['k','r','g','c','m','y','b',"orange","lime","pink"]
def onclick(event):
    x_value = round(event.xdata,2)
    y_value = round(event.ydata,2)

    if current_color == 'delete mode':
        delete_threshold = 1
        current_dist_squared = delete_threshold + 1
        for p in points:
            dist_squared = (x_value - p.x)**2 + (y_value - p.y)**2
            if dist_squared < current_dist_squared:
                current_dist_squared = dist_squared
                current_point = p
        if current_dist_squared < delete_threshold:
            current_point.mpl_point.remove()
            points.remove(current_point)
            print("removed point", (current_point.x, current_point.y, current_point.c), "dist_squared:", current_dist_squared)
        else:
            print("no point removed.")
    else:
        mpl_point = plt.scatter(x_value, y_value, c = colormap[current_color])

        points.append(point(x_value, y_value, current_color,mpl_point))
        print("plot point", (x_value, y_value), "with color:", current_color)

    fig.canvas.draw()

def on_press(event):
    global current_color
    key = event.key
    if key == 'delete':
        current_color = 'delete mode'
        print("current color changed to:", current_color)
    if key.isdigit():
        current_color = int(key)
        print("current color changed to:", current_color)

def on_close(event):
    print("Close figure!")
    files = [('2D-Spaces', '*.2dspace'),
            ('All Files', '*')]
    file = filedialog.asksaveasfile(filetypes = files, defaultextension = ".2dspace")
    file_path = file.name
    print(file_path)
    with open(file_path, 'a') as writer:
        for i in range(len(points)):
            writer.write(f"{points[i].x},{points[i].y},{points[i].c}\n")


global current_color 
current_color = 0

print("""Use digit-keys (0 to 9) to set colors. (default: 0). Click on the canvas to create a new point. 
Use 'del' to go into delete mode: Click nearby a point to delete it. Exit delete mode by choosing a color.
2dspace can be written to a file after closing the windows.""")

fig.canvas.mpl_connect('key_press_event', on_press)
fig.canvas.mpl_connect('button_press_event', onclick)
fig.canvas.mpl_connect('close_event', on_close)
plt.show()
