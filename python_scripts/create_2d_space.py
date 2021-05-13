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


x = []
y = []
c = []

colormap = ['k','r','g','c','m','y','b',"orange","lime","pink"]
def onclick(event):
    x_value = round(event.xdata,2)
    y_value = round(event.ydata,2)
    x.append(x_value)
    y.append(y_value)
    c.append(current_color)
    print("plot point", (x_value, y_value), "with color:", current_color)
    plt.scatter(x_value, y_value, c = colormap[current_color])
    fig.canvas.draw()

def on_press(event):
    global current_color
    key = event.key
    if key.isdigit():
        current_color = int(key)
        print("current color changed to:", current_color)

def on_close(event):
    print("Close figure!")
    files = [('2D-Spaces', '*.2dspace'),
            ('All Files', '*')]
    file = filedialog.asksaveasfile(filetypes = files, defaultextension = ".2dspace")
    file_path = file.name;
    print(file_path)
    with open(file_path, 'a') as writer:
        for i in range(len(x)):
            writer.write(f"{x[i]},{y[i]},{c[i]}\n")


global current_color 
current_color = 0

fig.canvas.mpl_connect('key_press_event', on_press)
fig.canvas.mpl_connect('button_press_event', onclick)
fig.canvas.mpl_connect('close_event', on_close)
plt.show()
