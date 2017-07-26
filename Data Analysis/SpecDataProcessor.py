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
import sys
import os
import time
# local file imports
import MathHelper # MathHelper.py
import GraphHelper # GraphHelper.py

spec_flr_offset = 0.12

flr_fs_start = 75
flr_fs_end = 100
flr_fm_start = 150
flr_fm_end = 200
flr_f0p_start = 550
flr_f0p_end = 600

ecs_baseline_start = 200
ecs_baseline_end = 250
ecs_decay_start = 250
ecs_decay_end = 275
ecs_decay_end_dcmu = 500

p700_fit_start = 2500
p700_fit_end = 2550
p700_pt = 6

vO2_probe_interval = 1 # 1 measurement per second
vO2_fit_margin = 30
vO2_phi2_timepoints = [1,3]

def process_experiment(experiment_dictionary):
	"""
	experiment_dictionary needs to be in the following format:
{
	"folder":folder path,
	"name":name of experiment,
	"experiment":base name of experiment (excludes the rep. #),
	"repetition":rep number of experiment,
	"script_notes":spec script notes (includes spec script itself)
	"O2_data":pandas DataFrame of the O2 probe CSV file,
	"timepoints":[list of timepoints that have 1 or more spec measurements],
	"flr_data":python dictionary{
		0:t0 flr dat file as pandas DataFrame,
		1:t1 flr dat file as pandas DataFrame,
		2:t2 flr dat file as pandas DataFrame,
		... etc.
	},
	"ecs_data":python dictionary{
		1:t1 ECS dat file as pandas DataFrame,
		3:t3 ECS dat file as pandas DataFrame,
		... etc.
	},
	"p700_data":python dictionary{
		6:t6 p700 dat file as pandas DataFrame,
		... etc.
	}
}
	The following will be added to the experiment_dictionary:
{
	...
	"vO2_per_Phi2":[list of O2 per phi2 values, NaN means no data],
	"vO2":[list of vO2 values],
	"vO2_light":[list of vO2 values],
	"vO2_dark":[list of vO2 values],
	"Phi2":[list of Phi2 values, NaN means no data],
	"Fs":[list of Fs values, NaN means no data],
	"Fm":[list of Fm values, NaN means no data],
	"F0p":[list of F0' values, NaN means no data],
	"vH+":[list of vH+ values, NaN means no data],
	"gH+":[list of gH+ values, NaN means no data],
	"ECSt":[list of ECS tau values, NaN means no data],
	"ECSamp":[list of ECS amplitude values, NaN means no data],
	"p700t":p700 re-reduction first-order decay fit, tau parameter,
	"p700k":p700 re-reduction first-order decay fit, k parameter (1/tau),
	"p700amp":p700 re-reduction amplitude of change
}
	"""
	# vO2 (needs user interaction)
	y = extract_column_from_dataframe(1,experiment_dictionary["O2_data"])
	x = MathHelper.scalar_multiply( list(range(0,len(y))), vO2_probe_interval )
	vO2_lines = GraphHelper.ask_user_for_slopes(x, y, vO2_fit_margin)
	vO2 = []
	vO2_light = []
	vO2_dark = []
	vO2_counter = 0
	for n in range(0, len(vO2_lines), 2):
		if(n+1 == len(vO2_lines)):
			# odd number of lines (assume last one is an accident)
			continue
		L = vO2_lines[n][0]
		D = vO2_lines[n+1][0]
		vO2_light.append(L)
		vO2_dark.append(D)
		vO2.append(L-D)
	experiment_dictionary["vO2"] = vO2
	experiment_dictionary["vO2_light"] = vO2_light
	experiment_dictionary["vO2_dark"] = vO2_dark
	# fluorescence
	num_points = len(experiment_dictionary["timepoints"])
	nan = float("nan")
	Phi2 = [nan] * num_points
	vO2perPhi2 = [nan] * num_points
	Fs = [nan] * num_points
	Fm = [nan] * num_points
	F0p = [nan] * num_points
	vH = [nan] * num_points
	gH = [nan] * num_points
	ECSt = [nan] * num_points
	ECSamp = [nan] * num_points
	p700t = [nan] * num_points
	p700k = [nan] * num_points
	p700amp = [nan] * num_points
	# time points
	for n in experiment_dictionary["timepoints"]:
		# check if this timepoint has a fluorescence measurement
		if n in experiment_dictionary["flr_data"]:
			# it does
			x = extract_column_from_dataframe(0, experiment_dictionary["flr_data"][n])
			y = extract_column_from_dataframe(1, experiment_dictionary["flr_data"][n])
			Fs[n] = numpy.mean(y[flr_fs_start:flr_fs_end])
			Fm[n] = numpy.mean(y[flr_fm_start:flr_fm_end])
			F0p[n] = numpy.mean(y[flr_f0p_start:flr_f0p_end])
			Phi2[n] = (Fm[n] - Fs[n]) / Fm[n]
			if n in vO2_phi2_timepoints:
				vO2perPhi2[n] = vO2[vO2_counter] / Phi2[n]
				vO2_counter += 1
		if n in experiment_dictionary["ecs_data"]:
			x = extract_column_from_dataframe(0, experiment_dictionary["ecs_data"][n])
			y = extract_column_from_dataframe(3, experiment_dictionary["ecs_data"][n])
			start = ecs_decay_start
			end = ecs_decay_end
			if(n == p700_pt):
				# ECS is super slow with DCMU
				end = ecs_decay_end_dcmu
			x0 = x[start]
			# define a function for fitting
			def f(x, coefficients):
				# Y=A*e^((-1/tau)*(x-x0))+y0:
				# coefficients: [A, tau, y0]
				return coefficients[0] * exp((-1/coefficients[1]) * (x - x0)) + coefficients[2]
			coefficient_vector = [ (y[start] - y[end]), (x[end] - x[start]) / 10, y[end] ]
			fitter_x = x[start:end]
			coefficient_vector, precision_vector, iteration_count = MathHelper.fit_function_bruteforce(fitter_x, y[start:end], f, coefficient_vector, 250)
			fitter_y = []
			for fit_x in fitter_x:
				fitter_y.append(f(fit_x,coefficient_vector))
			GraphHelper.multi_plot([str(experiment_dictionary["name"])+" ECS data","fit"], x, y, fitter_x, fitter_y)
			experiment_dictionary["ecs_data_fit"] = {}
			experiment_dictionary["ecs_data_fit"][n] = pandas.DataFrame({0:fitter_x, 1:fitter_y})
			#
			ECSt[n] = coefficient_vector[1]
			ECSamp[n] = coefficient_vector[0]
			gH[n] = 1.0 / ECSt[n]
			vH[n] = ECSamp[n] / ECSt[n]
		if n in experiment_dictionary["p700_data"]:
			GraphHelper.multi_plot([str(experiment_dictionary["name"])+" p700 flr. control"], extract_column_from_dataframe(0, experiment_dictionary["flr_data"][n]), extract_column_from_dataframe(1, experiment_dictionary["flr_data"][n]))
			x = extract_column_from_dataframe(0, experiment_dictionary["p700_data"][n])
			y = extract_column_from_dataframe(3, experiment_dictionary["p700_data"][n])
			x0 = x[p700_fit_start]
			# define a function for fitting
			def f(x, coefficients):
				# Y=A*e^((-1/tau)*(x-x0))+y0:
				# coefficients: [A, tau, y0]
				return -1 * coefficients[0] * exp((-1/coefficients[1]) * (x - x0)) + coefficients[2]
			coefficient_vector = [ (y[p700_fit_end] - y[p700_fit_start]), (x[p700_fit_end] - x[p700_fit_start]) / 5, y[p700_fit_end] ]
			fitter_x = x[p700_fit_start:p700_fit_end]
			coefficient_vector, precision_vector, iteration_count = MathHelper.fit_function_bruteforce(fitter_x, y[p700_fit_start:p700_fit_end], f, coefficient_vector, 250)
			fitter_y = []
			for fit_x in fitter_x:
				fitter_y.append(f(fit_x,coefficient_vector))
			GraphHelper.multi_plot([str(experiment_dictionary["name"])+" p700 data","fit"], x, y, fitter_x, fitter_y)
			experiment_dictionary["p700_data_fit"] = {}
			experiment_dictionary["p700_data_fit"][n] = pandas.DataFrame({0:fitter_x, 1:fitter_y})
			p700t[n] = coefficient_vector[1]
			p700k[n] = 1 / p700t[n]
			p700amp[n] = coefficient_vector[0]
	# data processed. Store it....
	experiment_dictionary["vO2_per_Phi2"] = vO2perPhi2
	experiment_dictionary["vO2"] = vO2
	experiment_dictionary["vO2_light"] = vO2_light
	experiment_dictionary["vO2_dark"] = vO2_dark
	experiment_dictionary["Phi2"] = Phi2
	experiment_dictionary["Fs"] = Fs
	experiment_dictionary["Fm"] = Fm
	experiment_dictionary["F0p"] = F0p
	experiment_dictionary["vH+"] = vH
	experiment_dictionary["gH+"] = gH
	experiment_dictionary["ECSt"] = ECSt
	experiment_dictionary["ECSamp"] = ECSamp
	experiment_dictionary["p700t"] = p700t
	experiment_dictionary["p700k"] = p700k
	experiment_dictionary["p700amp"] = p700amp
	# return the dictionary
	return experiment_dictionary

def average_column_data(dataframe, column_index, start, end):
	"""
averages together a range of rows in a column from a pandas DataFrame
	"""
	sum = 0
	column = extract_column_from_dataframe(column_index, df)
	for row in range(start,end):
		sum += column[row]
	return sum / (end - start)

def extract_column_from_dataframe(column, df):
	"""
Returns a single column frmo a pandas DataFrame as a list
	"""
	if isinstance( column, int ):
		column_address = list(df)[column]
	else:
		column_address = column
	num_rows = df.shape[0]
	column_data = []
	for row in range(0,num_rows):
		column_data.append(df[column_address][row])
	return column_data
	
def exp(d):
	# math.exp can throw an exception for large errors
	try:
		return math.exp(d)
	except:
		return float("inf")