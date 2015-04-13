# scatterplot.py

import numpy as np
import pylab as pl

# Maxe an array of x values
x = [1, 2, 3, 4, 5]
# Make an array of y values
y = [1, 4, 9, 16, 25]

# Use pylab to plot x and y as red circles
pl.plot(x, y, 'ro')

# show the plot on the screen
pl.show()
