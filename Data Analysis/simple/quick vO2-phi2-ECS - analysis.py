"""
This python script is self-contained (does not require other .py files),
but requires the scientific python packages (SciPy, Numpy, Pandas, matplotlib)
"""

if __name__ == "__main__":
	# this little code block is just to show a GUI thing to make it clear 
	# that it is working (the scipy library imports may take a whole minute
	import tkinter
	tk_init_window = tkinter.Tk()
	label = tkinter.Label( tk_init_window, text="Loading scientific python libraries..." )
	label.pack()
	tk_init_window.update()

# imports
print('Loading...', end="")
# scipy/anaconda imports
import pandas
from scipy import stats
import numpy
from matplotlib import pyplot
from matplotlib import animation
# python standard library imports
import math
import copy
import statistics
import tkinter
from tkinter import filedialog
from tkinter import messagebox
import collections
import sys
import os
import time
import traceback
import re # regex library

nan = float('nan')
#
print('...Done!')

### Constants ###
flr_num_points = sum( [100,100,200,100,100] )
flr_column = 1 # sample detector channel
flr_fs_start = 75
flr_fs_end = 100
flr_fm_start = 150
flr_fm_end = 200
flr_f0p_start = 550
flr_f0p_end = 600
flr_blank_offset = 0.12

ecs_num_points = sum( [250,250,250] )
ecs_column = 3 # delta-A channel
ecs_baseline_start = 200
ecs_baseline_end = 250
ecs_decay_start = 250
ecs_decay_end_no_dcmu = 275
ecs_decay_end_dcmu = 500

p700_num_points = sum( [500,1000,1000,1000] )
p700_column = 3 # delta-A channel
p700_fit_start = 2500
p700_fit_end = 2600
p700_baseline_start = 2400
p700_baseline_end = 2499
p700_darkSS_start = 2600
p700_darkSS_end = 2700
p700_graph_start = p700_fit_start - 250
p700_graph_end = p700_darkSS_end + 100

vO2_column = 1
vO2_peak_search_start = 14*60
vO2_peak_search_end = vO2_peak_search_start + 600
vO2_peak_search_size = 30
vO2_peak_fit_size = 60
vO2_light_prod_relative_to_peak = -120
vO2_light_resp_relative_to_peak = 80
vO2_dark_resp_relative_to_peak = -330

#################


def main():
	#root_folder = os.path.abspath('.') # current folder (might be different from script file location)
	#root_folder = os.path.dirname(os.path.abspath(__file__)) # folder with this script file in it
	root_folder = filedialog.askdirectory() # ask user for location
	
	print('Scanning folder: ',root_folder)
	#
	for folder_path, directories_list, file_list in os.walk(root_folder):
		for file_name in file_list:
			try:
				file_path = str(os.path.join(folder_path, file_name))
				if(file_path.endswith('0001.dat')):
					# fluorescence Phi2 measurement
					process_phi2_file(file_path)
				if(file_path.endswith('0003.dat')):
					# ECS DIRK measurement
					process_ECSDIRK_file(file_path)
				if(file_path.endswith('0005.dat')):
					# p700 LIRK measurement
					process_p700_file(file_path)
				if(file_path.endswith('.csv')):
					# vO2 measurement
					process_vO2_file(file_path)
			except Exception as ex:
				print("\tError: "+str(ex))
				print(traceback.format_exc())
	#
	exit()

def process_phi2_file(file_path):
	print('Processing fluorescence data file',file_path)
	data = read_data_file(file_path)
	num_measurements = len(data[flr_column]) // int(flr_num_points)
	_Fs = []
	_Fm = []
	_F0p = []
	_Phi2 = []
	_index = []
	strain, rep, hr = get_meta_from_name(file_path)
	for measurement_index in range(0,num_measurements):
		m_start = measurement_index * flr_num_points
		m_end = (measurement_index + 1) * flr_num_points
		data_slice = data[flr_column][ m_start : m_end ]
		Fs = numpy.mean(data_slice[flr_fs_start:flr_fs_end]) - flr_blank_offset
		Fm = numpy.mean(data_slice[flr_fm_start:flr_fm_end]) - flr_blank_offset
		F0p = numpy.mean(data_slice[flr_f0p_start:flr_f0p_end]) - flr_blank_offset
		Phi2 = (Fm - Fs) / Fm
		_Fs.append(Fs)
		_Fm.append(Fm)
		_F0p.append(F0p)
		_Phi2.append(Phi2)
		_index.append(measurement_index)
	data_dict = {}
	data_dict['strain'] = [strain]*num_measurements
	data_dict['rep'] = [rep]*num_measurements
	data_dict['hour'] = [hr]*num_measurements
	data_dict['file'] = [file_path]*num_measurements
	data_dict['index'] = _index
	data_dict['Fs'] = _Fs
	data_dict['Fm'] = _Fm
	data_dict['F0p'] = _F0p
	data_dict['Phi2'] = _Phi2
	write_data_dictionary(data_dict, str(file_path)+'_out.txt')

def process_ECSDIRK_file(file_path):
	print('Processing ECS DIRK data file',file_path)
	data = read_data_file(file_path)
	strain, rep, hr = get_meta_from_name(file_path)
	num_measurements = len(data[ecs_column]) // int(ecs_num_points)
	_ECS_amp = []
	_ECS_tau = []
	_Rsqr = []
	_vH = []
	_gH = []
	_index = []
	for measurement_index in range(0,num_measurements):
		m_start = measurement_index * ecs_num_points
		m_end = (measurement_index + 1) * ecs_num_points
		y = data[ecs_column][ m_start : m_end ]
		t = data[0][ m_start : m_end ]
		x0 = t[ecs_decay_start]
		# define a function for fitting
		ecs_decay_end = ecs_decay_end_no_dcmu
		if('dcmu' in file_path):
			ecs_decay_end = ecs_decay_end_dcmu
		def f(x, coefficients):
			# Y=A*e^((-1/tau)*(x-x0))+y0:
			# coefficients: [A, tau, y0]
			return coefficients[0] * exp((-1/coefficients[1]) * (x - x0)) + coefficients[2]
		coefficient_vector = [ (y[ecs_decay_start] - y[ecs_decay_end]), (t[ecs_decay_end] - t[ecs_decay_start]) / 10, y[ecs_decay_end] ]
		coefficient_vector, precision_vector, r_sq, iteration_count = fit_function_bruteforce(t[ecs_decay_start:ecs_decay_end], y[ecs_decay_start:ecs_decay_end], f, coefficient_vector, 1000)
		ECS_tau = coefficient_vector[1]
		ECS_amp = coefficient_vector[0]
		gH = 1.0 / ECS_tau
		vH = ECS_amp / ECS_tau
		_ECS_tau.append(ECS_tau)
		_ECS_amp.append(ECS_amp)
		_Rsqr.append(r_sq)
		_vH.append(vH)
		_gH.append(gH)
		_index.append(measurement_index)
		fy = []
		fx = t[ecs_decay_start:ecs_decay_end]
		for _x in fx:
			fy.append(f(*(_x, coefficient_vector)))
		save_fit_graph(t, y, fx, fy, str(file_path)+'_fit#'+str(measurement_index+1)+'.pdf')
	data_dict = {}
	data_dict['strain'] = [strain]*num_measurements
	data_dict['rep'] = [rep]*num_measurements
	data_dict['hour'] = [hr]*num_measurements
	data_dict['file'] = [file_path]*num_measurements
	data_dict['index'] = _index
	data_dict['ECSamp'] = _ECS_amp
	data_dict['ECSt'] = _ECS_tau
	data_dict['R-squared'] = _Rsqr
	data_dict['gH+'] = _gH
	data_dict['vH+'] = _vH
	write_data_dictionary(data_dict, str(file_path)+'_out.txt')

def process_p700_file(file_path):
	print('Processing p700 LIRK data file',file_path)
	data = read_data_file(file_path)
	strain, rep, hr = get_meta_from_name(file_path)
	num_measurements = len(data[p700_column]) // int(p700_num_points)
	_p700_tau = []
	_p700_k = []
	_p700_amp = []
	_Rsqr = []
	_p700_deltaA = []
	_index = []
	for measurement_index in range(0,num_measurements):
		m_start = measurement_index * p700_num_points
		m_end = (measurement_index + 1) * p700_num_points
		y = data[p700_column][ m_start : m_end ]
		t = data[0][ m_start : m_end ]
		_index.append(measurement_index)
		# adjust for drift
		slope, y_intercept, r_squared, p_value, std_err = fit_line(t[p700_baseline_start:p700_baseline_end], y[p700_baseline_start:p700_baseline_end])
		t, y = subtract_line(t, y, slope, y_intercept)
		# fit exponential decay (upside-down)
		x0 = t[p700_fit_start]
		# define a function for fitting
		def f(x, coefficients):
			# Y=A*e^((-1/tau)*(x-x0))+y0:
			# coefficients: [A, tau, y0]
			return -1 * coefficients[0] * exp((-1/coefficients[1]) * (x - x0)) + coefficients[2]
		coefficient_vector = [ (y[p700_fit_end] - y[p700_fit_start]), (t[p700_fit_end] - t[p700_fit_start]) / 5, y[p700_fit_end] ]
		coefficient_vector, precision_vector, r_sq, iteration_count = fit_function_bruteforce(t[p700_fit_start:p700_fit_end], y[p700_fit_start:p700_fit_end], f, coefficient_vector, 1000)
		p700_tau = coefficient_vector[1]
		p700_k = 1 / p700_tau
		p700_amp = coefficient_vector[0]
		dA0 = numpy.mean(y[p700_baseline_start:p700_baseline_end])
		dA1 = numpy.mean(y[p700_darkSS_start:p700_darkSS_end])
		p700_deltaA = dA1 - dA0
		_p700_tau.append(p700_tau)
		_p700_k.append(p700_k)
		_p700_amp.append(p700_amp)
		_p700_deltaA.append(p700_deltaA)
		_Rsqr.append(r_sq)
		fy = []
		fx = t[p700_fit_start:p700_fit_end]
		for _x in fx:
			fy.append(f(*(_x, coefficient_vector)))
		save_fit_graph(t[p700_graph_start:p700_graph_end], y[p700_graph_start:p700_graph_end], fx, fy, str(file_path)+'_fit#'+str(measurement_index+1)+'.pdf')
	data_dict = {}
	data_dict['strain'] = [strain]*num_measurements
	data_dict['rep'] = [rep]*num_measurements
	data_dict['hour'] = [hr]*num_measurements
	data_dict['file'] = [file_path]*num_measurements
	data_dict['index'] = _index
	data_dict['p700-t'] = _p700_tau
	data_dict['p700-k'] = _p700_k
	data_dict['p700-amp'] = _p700_amp
	data_dict['R-squared'] = _Rsqr
	data_dict['dA703'] = _p700_deltaA
	write_data_dictionary(data_dict, str(file_path)+'_out.txt')

def process_vO2_file(file_path):
	print('Processing vO2 data file',file_path)
	data = read_data_file(file_path, separator=',')
	if(data[vO2_column][0] == 'Oxygen'):
		# oxygen data
		y = data[vO2_column][1:]
		x = list(range(0,len(y)))
		data_dict = {}
		strain, rep, hr = get_meta_from_name(file_path)
		data_dict['strain'] = [strain]
		data_dict['rep'] = [rep]
		data_dict['hour'] = [hr]
		data_dict['file'] = [file_path]
		# find the maximum right after the light pulse
		delta = int(vO2_peak_search_size / 2)
		search_index = vO2_peak_search_start
		peak_index = search_index
		peak_value = 0
		while(search_index < vO2_peak_search_end):
			ave = numpy.mean(y[search_index-delta:search_index+delta])
			if(ave > peak_value):
				peak_index = search_index
				peak_value = ave
			search_index += 1
		# fit the slopes relative to the peak
		delta = int(vO2_peak_fit_size / 2)
		dark_resp_start = peak_index + vO2_dark_resp_relative_to_peak - delta
		dark_resp_end = peak_index + vO2_dark_resp_relative_to_peak + delta
		light_prod_start = peak_index + vO2_light_prod_relative_to_peak - delta
		light_prod_end = peak_index + vO2_light_prod_relative_to_peak + delta
		light_resp_start = peak_index + vO2_light_resp_relative_to_peak - delta
		light_resp_end = peak_index + vO2_light_resp_relative_to_peak + delta
		# dark(-ish) acclimated respiration
		slope1, y_intercept1, r_squared1, p_value1, std_err1 = fit_line(x[dark_resp_start:dark_resp_end], y[dark_resp_start:dark_resp_end])
		data_dict['dark-resp rate']=[slope1]
		data_dict['dark-resp R-squared']=[r_squared1]
		# light oxygen production
		slope2, y_intercept2, r_squared2, p_value2, std_err2 = fit_line(x[light_prod_start:light_prod_end], y[light_prod_start:light_prod_end])
		data_dict['light-prod rate']=[slope2]
		data_dict['light-prod R-squared']=[r_squared2]
		# post-light respiration
		slope3, y_intercept3, r_squared3, p_value3, std_err3 = fit_line(x[light_resp_start:light_resp_end], y[light_resp_start:light_resp_end])
		data_dict['light-resp rate']=[slope3]
		data_dict['light-resp R-squared']=[r_squared3]
		# save graph
		fit_y = []
		for i in range(0,len(x)):
			t = x[i]
			if(i >= dark_resp_start and i < dark_resp_end):
				fit_y.append(t * slope1 + y_intercept1)
			elif(i >= light_prod_start and i < light_prod_end):
				fit_y.append(t * slope2 + y_intercept2)
			elif(i >= light_resp_start and i < light_resp_end):
				fit_y.append(t * slope3 + y_intercept3)
			else:
				fit_y.append(nan)
		save_fit_graph(x, y, x, fit_y, str(file_path)+'_fit.pdf')
		# save results
		data_dict['strain'] = [strain]
		data_dict['rep'] = [rep]
		data_dict['hour'] = [hr]
		data_dict['file'] = [file_path]
		write_data_dictionary(data_dict, str(file_path)+'_out.txt')
	else:
		# not oxygen data, do nothing
		return


color_sequence = ['black', 'red', 'yellow', 'blue', '#666666', '#664400', '#ee6600', 'green', '#6600ee', '#aaaaaa']
def save_fit_graph(x_data, y_data, fit_x_data, fit_y_data, file_path, x_label='X', y_label='Y', legend_location='upper right'):
	print('Saving data fit figure', file_path)
	figure, axese = pyplot.subplots()
	axese.set_ylabel(y_label)
	axese.set_xlabel(x_label)
	# add the data
	color = color_sequence[0 % len(color_sequence)]
	axese.plot(x_data, y_data,     linestyle='-', color=color, label='data')
	color = color_sequence[1 % len(color_sequence)]
	axese.plot(fit_x_data, fit_y_data, linestyle='-', color=color, label='fit')
	#
	handles, labels = axese.get_legend_handles_labels()
	figure.legend(handles=handles, labels=labels, loc=legend_location)
	
	pyplot.savefig(file_path)
	pyplot.close(figure)
	
def get_meta_from_name(file_path):
	strain = '?'
	rep = '?'
	hr = '?'
	base = str(os.path.basename(file_path))
	if('rep' in base):
		i1 = base.index('rep') + 3
		if('_' in base[i1:]):
			i2 = base.index('_',i1)
			rep = base[i1:i2]
	if('hr' in base):
		i1 = base.index('hr') + 2
		if('_' in base[i1:]):
			i2 = base.index('_',i1)
			hr = base[i1:i2]
		elif('.' in base[i1:]):
			i2 = base.index('.',i1)
			hr = base[i1:i2]
	if(re.search('[cC]+.*1009',base) != None):
		strain = 'CC-1009'
	if(re.search('[cC]+.*2343',base) != None):
		strain = 'CC-2343'
	return strain, rep, hr

def read_data_file(filepath, separator='\t'):
	"""
	reads a tab-delimited file and returns a 2D array in the format [column][row]
	"""
	col_count = 0
	row_count = 0
	rows = []
	with open(filepath, 'r') as f_in:
		for ln in f_in:
			cells = ln.split(separator)
			if(len(cells) > col_count):
				col_count = len(cells)
			rows.append(cells)
			row_count += 1
	# pivot from [row][col] to [col][row]
	data_table = []
	for col in range(0,col_count):
		column = []
		for row in range(0,row_count):
			if(len(rows[row]) > col):
				try:
					column.append(float(rows[row][col].strip()))
				except:
					column.append(rows[row][col].strip())
			else:
				column.append(nan)
		data_table.append(column)
	return data_table

def write_data_dictionary(data_dict, file_path, separator='\t', newline='\n'):
	print('Writing output to file',file_path)
	num_rows = 0
	for key in data_dict:
		if(len(data_dict[key]) > num_rows):
			num_rows = len(data_dict[key])
	with open(file_path, 'w') as f_out:
		first = True
		for key in data_dict:
			if(first == False):
				f_out.write(separator)
			first = False
			f_out.write(key)
		f_out.write(newline)
		row = 0
		while(row < num_rows):
			first = True
			for key in data_dict:
				if(first == False):
					f_out.write(separator)
				first = False
				if(row < len(data_dict[key])):
					f_out.write(str(data_dict[key][row]))
			f_out.write(newline)
			row += 1

def exp(n):
	# overflow error proof exponential
	try:
		return math.exp(n)
	except:
		return nan

def subtract_line(x_data, y_data, line_slope, line_intercept):
	"""
Subtracts a line from a data set, returning the new x and y values
	"""
	new_x = []
	new_y = []
	for i in range(0,len(x_data)):
		offset = line_slope * x_data[i] + line_intercept
		new_x.append(x_data[i])
		new_y.append(y_data[i] - offset)
	return new_x, new_y

def fit_line(x_data, y_data):
	"""
performs a linear fit to the data and return the slope, y-intercept, R-squared 
value, P value, and standard error for the fit. Use as follows:

// x and y are lists of numbers with more than 75 elements
// fitting points 25 through 75 in from the data
start = 25
end = 75
slope, y_intercept, r_squared, p_value, std_err = MathHelper.fit_line(x[start:end], y[start,end])
print( str.format("Fitted formula: Y = {a}X + {b}", a=slope, b=y_intercept))
print( str.format("\tR-squared = {r2}", r2=r_squared))

	"""
	slope, y_intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
	r_squared = r_value * r_value
	return slope, y_intercept, r_squared, p_value, std_err

def chi_squared_of_fit(y_data, fitted_y_data):
	"""
This function calculates and returnsthe Chi-squared value of a generated set 
of Y values (fitted_y_data) against the observed Y values (y_data).
	"""
	# note that in the typical definition of Chi-squared, it is assumed that 
	# nature is wrong and our formula is theoretically perfect, but here we 
	# are testing a model against emperical data, so the "expected" values 
	# are the data we measured and the "observed" values are the values 
	# generated by a fitting formula, and therefore we swap these variables 
	# in the stats.chisquare(...) invocation
	return stats.chisquare(fitted_y_data,y_data)


def r_squared(y_data,test_y_data):
	average = 0
	iii = 0
	size = len(y_data)
	for n in range(0,size):
		average += y_data[n]
		iii += 1
	average = average / iii
	sumResidual = 0
	sumTotal = 0
	for n in range(0,size):
		d = y_data[n]
		sim = test_y_data[n]
		sumResidual += (d - sim)*(d - sim)
		sumTotal += (d - average)*(d - average)
	return 1 - (sumResidual / sumTotal)


def _fit_err(x_data, y_data, formula_function, coefficient_vector):
	"""
quickly calculates a simple fitting error value (this is NOT standard error!)
	"""
	sum = 0;
	m = 1.0 / len(y_data)
	for n in range(0,len(x_data)):
		t = x_data[n]
		obs = y_data[n]
		sim = formula_function(*( t, coefficient_vector ) )
		er = sim - obs
		try:
			sum += er * er * m
		except ValueError:
			return float('inf')
	return sum

def fit_function_bruteforce(x_data, y_data, formula_function, coefficient_vector, max_iterations ):
	"""
This function uses a brute-force guess-and-check

x_data: the x values of the data set

y_data: the y values of the data set

formula_function: a function whose input parameters are (x, coefficient_vector)
                  and returns a single number by applying a formula to x with 
                  coefficients defined in coefficient_vector. For example, 
                  here's the fuction definition for an exponential decay fit 
                  where the coefficients are [A, tau, C] in the formula 
                  Y=A*e^(-x/tau)+C:

	def exponential_decay_function(x, cv_list):
		return cv_list[0] * math.exp( (-1 / cv_list[1]) x ) + cv_list[2]
	coefficients, precisions, iterations = fit_function_bruteforce(x, y, exponential_decay_function, [1,10,0], 1000000)
	print(str.format("fit: Y={A}*e^(-x/{tau})+{C}",A=coefficients[0], tau=coefficients[1], C=coefficients[2]))

coefficient_vector: list of numbers corresponding to coefficients in the 
                    formula which fit_function_bruteforce(...) will manipulate 
                    to try to fit he formula to the data. The starting values 
                    should be a best guess close to the actual value. If these 
                    values are too far off, the formula may get stick in a 
                    local maxima or edge-case

max_iterations: the maximum number of iterations through the fitting formula

returns coefficient_vector, precision, r_sq, iteration_count: 
	returns the coefficient_vector after adjusting it to achieve the best 
	possible fit within the allowed number of iterations, the +/- precision 
	of the fit for each coefficient (also as a list), the R-squared value of 
	the final fit to the data, and the actual number of iterations used to 
	perform the fit
	"""
	iterations = 0
	# initialize deltas to the coefficients
	delta_vector = scalar_multiply(copy.deepcopy(coefficient_vector), 0.25)
	while(iterations < max_iterations):
		# scale-down the deltas by a factor of 2 each time
		delta_vector = scalar_multiply(delta_vector, 0.5)
		new_cv, jiggles = _improveFit(x_data, y_data, formula_function, coefficient_vector, delta_vector, max_iterations - iterations)
		coefficient_vector = new_cv
		iterations += jiggles
	# done
	final_fit_y = []
	for t in x_data:
		final_fit_y.append(formula_function(*( t, coefficient_vector ) ))
	r_sq = r_squared(y_data, final_fit_y)
	return coefficient_vector, delta_vector, r_sq, iterations

def _improveFit(x,y,formula,cvec,delta, maxIterations):
	"""
jiggles the coefficients to improve the formula fit a little bit

x: x data

y: y data

formula: the fitting formula (see description for fit_function_bruteforce(...) )

cvec: coefficient vector (see description for fit_function_bruteforce(...) )

delta: list of jiggle sizes corresponding to cvec

maxIterations: maximum number of jiggles allowed before returning
	"""
	# adjust the variables by the delta amount in decrease the error value
	iterations = 0
	while True: # python does not support do-while loops
		lastErr = _fit_err(x,y,formula,cvec)
		for i in range(len(cvec)):
			oldC = cvec[i]
			upC = cvec[i]+delta[i]
			downC = cvec[i]-delta[i]
			# current fit error
			currentErr = _fit_err(x,y,formula,cvec)
			# increase the coefficient a little and check again
			cvec[i] = upC
			errPlus = _fit_err(x,y,formula,cvec)
			# decrease the coefficient a little and check again
			cvec[i] = downC
			errMinus = _fit_err(x,y,formula,cvec)
			if(errPlus < currentErr and errPlus < errMinus):
				# increase the variable
				cvec[i] = upC;
			elif(errMinus < currentErr):
				# decrease the variable
				cvec[i] = downC
			else:
				# no change
				cvec[i] = oldC
		iterations += 1
		if(lastErr <= _fit_err(x,y,formula,cvec) or iterations >= maxIterations):
			break
	return cvec, iterations

def scalar_multiply(vector_list, scalar):
	"""
Multiplies a vector (represented as a list of numbers) by a scalar value and 
returns the new vector (original vector values are not changed)s
	"""
	v = copy.deepcopy(vector_list)
	for i in range(0, len(v)):
		v[i] = v[i] * scalar
	return v

def p_value(set1, set2):
	"""
returns the T-test P-value for two independent sets of data
	"""
	s, p = stats.ttest_ind(set1, set2)
	return p

def mean(data_set):
	try:
		return statistics.mean(data_set)
	except:
		return nan

def stdev(data_set):
	if(len(data_set) < 2):
		return nan
	else:
		try:
			return statistics.stdev(data_set)
		except:
			return nan

if __name__ == "__main__":
	# clean-up of pop-up code
	tk_init_window.withdraw()
	tk_init_window.quit()
	#
	main()
