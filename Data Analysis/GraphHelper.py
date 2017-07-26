# imports
# scipy/anaconda imports
import pandas
from scipy import stats
import numpy
from matplotlib import pyplot
# python standard library imports
import tkinter
from tkinter import filedialog
import math
import collections
import time

color_sequence = ['black', 'red', 'yellow', 'blue', '#666666', '#664400', '#ee6600', 'green', '#6600ee', '#aaaaaa']

def multi_plot(name_list, *data, linestyle='-', x_label='x', y_label='y', legend_location='upper right'):
	"""
For example:

multi_plot(["series 1", "series 2", "series 3"], x1, y1, x2, y2, x3, y3)
	"""
	#WARNING: pyplot code is UGLY
	figure, axese = pyplot.subplots()
	axese.set_ylabel(y_label)
	axese.set_xlabel(x_label)
	# add the data
	for n in range(0,len(name_list)):
		color = color_sequence[n % len(color_sequence)]
		axese.plot(data[2*n], data[2*n+1], linestyle=linestyle, color=color, label=name_list[n])
	#
	handles, labels = axese.get_legend_handles_labels()
	figure.legend(handles=handles, labels=labels, loc=legend_location)
	
	pyplot.ioff() # ioff() to make figure.show() block the main thread, ion() to prevent blocking
	pyplot.show()

def plot_data(y_data):
	plot_xy_data( list( range(1,len(y_data)+1) ), y_data )

def plot_xy_data(x_data, y_data, linestyle='-', color='black'):
	"""
Displays a simple plot of the provided data

x_data: list of x values

y_data: list of y values

	"""
	#WARNING: pyplot code is UGLY
	figure, axese = pyplot.subplots()
	axese.plot(x_data, y_data, linestyle=linestyle, color=color)
	pyplot.ioff() # ioff() to make figure.show() block the main thread, ion() to prevent blocking
	pyplot.show()
	
	

def ask_user_for_slopes(x_data, y_data, fit_margin): #WARNING: pyplot code is UGLY
	"""
Shows the data to the user and then the user can pick where the slopes should 
be calculated by double-clicking on the graph. Right-clicking removes the last 
line that was drawn. When done, hitting enter closes the window and returns 
the fitted lines.

x_data: list of x values

y_data: list of y values

fit_margin: number of points in either direction from the click point to use 
            to perform the linear fits

Returns a list of tuples, where each tuple represents a linear fit in the 
format of (slope, intercept, r-squared, P-value, std-err)
	"""
	fit_lines = []
	line_plots = []
	global user_done
	user_done	= False
	def click(event):
		if event.button == 1 and event.dblclick == True:
			#print("adding line to graph: ",end="")
			x_coord = event.xdata
			x_index = 0
			while(x_index < len(x_data) and x_data[x_index] < x_coord):
				x_index += 1
			left = x_index - fit_margin
			if(left < 0):
				left = 0
			right = x_index + fit_margin
			if(right > len(x_data)):
				right = len(x_data)
			slope, intercept, r, prob, sterrest = stats.linregress(x_data[left:right],y_data[left:right])
			#print(str.format("Y = {m}X + {b}", m=slope, b=intercept))
			fit_y = []
			for nx in x_data[left:right]:
				fit_y.append(slope * nx + intercept)
			line, = axese.plot(x_data[left:right], fit_y, linestyle='-', color='red')
			figure.canvas.draw()
			line_plots.append(line)
			fit_lines.append( (slope, intercept, r*r, prob, sterrest) )
		if event.button > 1:
			#print("removing last line from graph")
			if(len(fit_lines) > 0):
				line_plots[len(line_plots)-1].remove()
				del line_plots[-1:]
				del fit_lines[-1:]
				figure.canvas.draw()
	def keyboard(event):
		if(str(event.key) == 'enter'):
			global user_done
			#print("done")
			pyplot.close(figure)
			user_done = True
	def closed(event):
		global user_done
		pyplot.close()
		user_done = True
	figure, axese = pyplot.subplots()
	figure.canvas.mpl_connect('button_press_event', click)
	figure.canvas.mpl_connect('key_press_event', keyboard)
	figure.canvas.mpl_connect('close_event', closed)
	axese.plot(x_data, y_data, linestyle='-', color='black')
	pyplot.ion() # ioff() to make figure.show() block the main thread, ion() to prevent blocking
	pyplot.show()

	
	# wait until user is done (pyplot is really stupid about its multi-threading)
	try:
		while not user_done:
			time.sleep(0.01)
			if not user_done:
				figure.canvas.flush_events()
			if not user_done:
				figure.canvas.draw()
			#figure.canvas.get_tk_widget().update() # will raise a _tkinter.TclError as soon as the window closes
	except Exception:
		pass # do nothing
	
	return fit_lines