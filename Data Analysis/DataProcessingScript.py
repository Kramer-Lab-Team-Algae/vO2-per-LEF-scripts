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
# local file imports
import FileHelper # FileHelper.py
import MathHelper # MathHelper.py
import SpecDataProcessor # SpecDataProcessor.py
import GraphHelper # GraphHelper.py
#
print('...Ready!')

if __name__ == "__main__":
	# clean-up of pop-up code
	tk_init_window.withdraw()
	tk_init_window.quit()

def main():
	# main function is called at end of script so that the main body of the 
	# script can be written above the function definitions
	setup()
	try: # try-finally block ensures graceful exit even if there's an error
		root_folder = filedialog.askdirectory()
		# scan through all of the sub-folders of the chosen folder and figure out 
		# which folders contain data, then add those folders to a list of folders 
		# to process
		experiment_folder_list = []
		for folder_path, directories_list, files_list in os.walk(root_folder):
			if(is_experiment_folder(folder_path)):
				experiment_folder_list.append(folder_path)
		#
		# now we have a list of experiment folders, now to parse the file 
		# names to generate a list of experiments to process
		experiment_list = []
		for experiment_folder in experiment_folder_list:
			experiment_list = experiment_list + get_experiments_from_folder(experiment_folder)
		# each experiment is a python dictionary with the following keys and values:
		"""
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
	}
	"ecs_data":python dictionary{
		1:t1 ECS dat file as pandas DataFrame,
		3:t3 ECS dat file as pandas DataFrame,
		... etc.
	}
	"p700_data":python dictionary{
		6:t6 p700 dat file as pandas DataFrame,
		... etc.
	}
}
		"""
		# data processing
		show_message("Sloper Picker Instructions","Double-click to add slope fit line, right-click to undo last line placement, and hit 'enter' when done.")
		results_list = []
		for experiment in experiment_list:
			try: 
				# use try-except to prevent an error in one experiment 
				# from stopping the others from being processed
				results = process_experiment(experiment)
				results_list.append(results)
			except Exception as ex:
				print("\tError: "+str(ex))
				print(traceback.format_exc())
		# collation and higher level analysis
		# collect the results by experiment name
		results_by_name = {}
		for results in results_list:
			name = results["experiment"]
			if not name in results_by_name:
				results_by_name[name] = []
			results_by_name[name].append(results)
		#####
		# output the analysis
		output_file_base = str(filedialog.asksaveasfilename())
		if output_file_base == "":
			return
		combined_output = open(output_file_base + "_averaged.txt","w")
		individual_output = open(output_file_base + ".txt","w")
		individual_output.write("folder\texperiment\treplicate\tvO2-ambient-light\tvO2-ambient-dark\tvO2-lowO2-light\tvO2-lowO2-dark\tvO2-ambient\tvO2-lowO2\tFv/Fm\tPhi2-ambient\tPhi2-lowO2\tvO2/Phi2-ambient\tvO2/Phi2-lowO2\tvECStau-ambient\tECStau-lowO2\tECSamp-ambient\tECSamp-lowO2\tgH+-ambient\tgH+-lowO2\tvH+-ambient\tvH+-lowO2\tp700t\tp700k\tp700amp\t\n")
		combined_output.write("experiment\tnum_replicate\tvO2-ambient-light\tSD\tvO2-ambient-dark\tSD\tvO2-lowO2-light\tSD\tvO2-lowO2-dark\tSD\tvO2-ambient\tSD\tvO2-lowO2\tSD\tP-value\tFv/Fm\tSD\tPhi2-ambient\tSD\tPhi2-lowO2\tSD\tP-value\tvO2/Phi2-ambient\tSD\tvO2/Phi2-lowO2\tSD\tP-value\tvECStau-ambient\tSD\tECStau-lowO2\tSD\tP-value\tECSamp-ambient\tSD\tECSamp-lowO2\tSD\tP-value\tgH+-ambient\tSD\tgH+-lowO2\tSD\tP-value\tvH+-ambient\tSD\tvH+-lowO2\tSD\tP-value\tp700t\tSD\tp700k\tSD\tp700amp\tSD\t\n")
		vO2_phi2_timepoints = SpecDataProcessor.vO2_phi2_timepoints
		ambO2_pt = vO2_phi2_timepoints[0]
		lowO2_pt = vO2_phi2_timepoints[1]
		p700_pt = 6
		#
		for name in results_by_name:
			try: 
				# lists for averaging and stats
				# (also order of columns)
				vO2_light_amb_list = []
				vO2_dark_amb_list = []
				vO2_light_low_list = []
				vO2_dark_low_list = []
				vO2_amb_list = []
				vO2_low_list = []
				# P-value for difference between ambient and low O2
				Phi2_FvFm_list = []
				Phi2_amb_list = []
				Phi2_low_list = []
				# P-value for difference between ambient and low O2
				vO2_per_Phi2_amb_list = []
				vO2_per_Phi2_low_list = []
				# P-value for difference between ambient and low O2
				ECSt_amb_list = []
				ECSt_low_list = []
				# P-value for difference between ambient and low O2
				ECSamp_amb_list = []
				ECSamp_low_list = []
				# P-value for difference between ambient and low O2
				gH_amb_list = []
				gH_low_list = []
				# P-value for difference between ambient and low O2
				vH_amb_list = []
				vH_low_list = []
				# P-value for difference between ambient and low O2
				p700t_list = []
				p700k_list = []
				p700amp_list = []
				#experiment
				combined_output.write(str(name)+"\t")
				rep_count = len(results_by_name[name])
				for xp in results_by_name[name]:
					individual_output.write(str(xp["folder"])+"\t")
					individual_output.write(str(name)+"\t")
					# rep
					rep_name = xp["repetition"]
					individual_output.write(str(rep_name)+"\t")
					#vO2
					vO2_light_amb = xp["vO2_light"][0]
					vO2_light_amb_list.append(vO2_light_amb)
					individual_output.write(str(vO2_light_amb)+"\t")
					vO2_dark_amb = xp["vO2_dark"][0]
					vO2_dark_amb_list.append(vO2_dark_amb)
					individual_output.write(str(vO2_dark_amb)+"\t")
					vO2_light_low = xp["vO2_light"][1]
					vO2_light_low_list.append(vO2_light_low)
					individual_output.write(str(vO2_light_low)+"\t")
					vO2_dark_low = xp["vO2_dark"][1]
					vO2_dark_low_list.append(vO2_dark_low)
					individual_output.write(str(vO2_dark_low)+"\t")

					vO2_amb = xp["vO2"][0]
					vO2_amb_list.append(vO2_amb)
					individual_output.write(str(vO2_amb)+"\t")
					vO2_low = xp["vO2"][1]
					vO2_low_list.append(vO2_low)
					individual_output.write(str(vO2_low)+"\t")
					#Phi2
					Phi2_FvFm = xp["Phi2"][0]
					Phi2_FvFm_list.append(Phi2_FvFm)
					individual_output.write(str(Phi2_FvFm)+"\t")
					Phi2_amb = xp["Phi2"][ambO2_pt]
					Phi2_amb_list.append(Phi2_amb)
					individual_output.write(str(Phi2_amb)+"\t")
					Phi2_low = xp["Phi2"][lowO2_pt]
					Phi2_low_list.append(Phi2_low)
					individual_output.write(str(Phi2_low)+"\t")
					#vO2/Phi2
					vO2_per_Phi2_amb = xp["vO2_per_Phi2"][ambO2_pt]
					vO2_per_Phi2_amb_list.append(vO2_per_Phi2_amb)
					individual_output.write(str(vO2_per_Phi2_amb)+"\t")
					vO2_per_Phi2_low = xp["vO2_per_Phi2"][lowO2_pt]
					vO2_per_Phi2_low_list.append(vO2_per_Phi2_low)
					individual_output.write(str(vO2_per_Phi2_low)+"\t")
					# ECS tau
					ECSt_amb = xp["ECSt"][ambO2_pt]
					ECSt_amb_list.append(ECSt_amb)
					individual_output.write(str(ECSt_amb)+"\t")
					ECSt_low = xp["ECSt"][lowO2_pt]
					ECSt_low_list.append(ECSt_low)
					individual_output.write(str(ECSt_low)+"\t")
					# ECS amplitude
					ECSamp_amb = xp["ECSamp"][ambO2_pt]
					ECSamp_amb_list.append(ECSamp_amb)
					individual_output.write(str(ECSamp_amb)+"\t")
					ECSamp_low = xp["ECSamp"][lowO2_pt]
					ECSamp_low_list.append(ECSamp_low)
					individual_output.write(str(ECSamp_low)+"\t")
					# gH+
					gH_amb = xp["gH+"][ambO2_pt]
					gH_amb_list.append(gH_amb)
					individual_output.write(str(gH_amb)+"\t")
					gH_low = xp["gH+"][lowO2_pt]
					gH_low_list.append(gH_low)
					individual_output.write(str(gH_low)+"\t")
					# vH+
					vH_amb = xp["vH+"][ambO2_pt]
					vH_amb_list.append(vH_amb)
					individual_output.write(str(vH_amb)+"\t")
					vH_low = xp["vH+"][lowO2_pt]
					vH_low_list.append(vH_low)
					individual_output.write(str(vH_low)+"\t")
					# p700
					p700t = xp["p700t"][p700_pt]
					p700t_list.append(p700t)
					individual_output.write(str(p700t)+"\t")
					p700k = xp["p700k"][p700_pt]
					p700k_list.append(p700k)
					individual_output.write(str(p700k)+"\t")
					p700amp = xp["p700amp"][p700_pt]
					p700t_list.append(p700amp)
					individual_output.write(str(p700amp)+"\t")
					# done
					individual_output.write("\n")
				# combined output
				combined_output.write(str(rep_count)+"\t")
				# vO2
				vO2_light_amb_ave = MathHelper.mean(vO2_light_amb_list)
				vO2_light_amb_dev = MathHelper.stdev(vO2_light_amb_list)
				combined_output.write(str(vO2_light_amb_ave)+"\t"+str(vO2_light_amb_dev)+"\t")
				vO2_dark_amb_ave = MathHelper.mean(vO2_dark_amb_list)
				vO2_dark_amb_dev = MathHelper.stdev(vO2_dark_amb_list)
				combined_output.write(str(vO2_dark_amb_ave)+"\t"+str(vO2_dark_amb_dev)+"\t")
				vO2_light_low_ave = MathHelper.mean(vO2_light_low_list)
				vO2_light_low_dev = MathHelper.stdev(vO2_light_low_list)
				combined_output.write(str(vO2_light_low_ave)+"\t"+str(vO2_light_low_dev)+"\t")
				vO2_dark_low_ave = MathHelper.mean(vO2_dark_low_list)
				vO2_dark_low_dev = MathHelper.stdev(vO2_dark_low_list)
				combined_output.write(str(vO2_dark_low_ave)+"\t"+str(vO2_dark_low_dev)+"\t")
				vO2_amb_ave = MathHelper.mean(vO2_amb_list)
				vO2_amb_dev = MathHelper.stdev(vO2_amb_list)
				combined_output.write(str(vO2_amb_ave)+"\t"+str(vO2_amb_dev)+"\t")
				vO2_low_ave = MathHelper.mean(vO2_low_list)
				vO2_low_dev = MathHelper.stdev(vO2_low_list)
				combined_output.write(str(vO2_low_ave)+"\t"+str(vO2_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(vO2_amb_list, vO2_low_list))+"\t")
				# Phi2
				Phi2_FvFm_ave = MathHelper.mean(Phi2_FvFm_list)
				Phi2_FvFm_dev = MathHelper.stdev(Phi2_FvFm_list)
				combined_output.write(str(Phi2_FvFm_ave)+"\t"+str(Phi2_FvFm_dev)+"\t")
				Phi2_amb_ave = MathHelper.mean(Phi2_amb_list)
				Phi2_amb_dev = MathHelper.stdev(Phi2_amb_list)
				combined_output.write(str(Phi2_amb_ave)+"\t"+str(Phi2_amb_dev)+"\t")
				Phi2_low_ave = MathHelper.mean(Phi2_low_list)
				Phi2_low_dev = MathHelper.stdev(Phi2_low_list)
				combined_output.write(str(Phi2_low_ave)+"\t"+str(Phi2_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(Phi2_amb_list, Phi2_low_list))+"\t")
				# vO2/Phi2
				vO2_per_Phi2_amb_ave = MathHelper.mean(vO2_per_Phi2_amb_list)
				vO2_per_Phi2_amb_dev = MathHelper.stdev(vO2_per_Phi2_amb_list)
				combined_output.write(str(vO2_per_Phi2_amb_ave)+"\t"+str(vO2_per_Phi2_amb_dev)+"\t")
				vO2_per_Phi2_low_ave = MathHelper.mean(vO2_per_Phi2_low_list)
				vO2_per_Phi2_low_dev = MathHelper.stdev(vO2_per_Phi2_low_list)
				combined_output.write(str(vO2_per_Phi2_low_ave)+"\t"+str(vO2_per_Phi2_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(vO2_per_Phi2_amb_list, vO2_per_Phi2_low_list))+"\t")
				# ECS tau
				ECSt_amb_ave = MathHelper.mean(ECSt_amb_list)
				ECSt_amb_dev = MathHelper.stdev(ECSt_amb_list)
				combined_output.write(str(ECSt_amb_ave)+"\t"+str(ECSt_amb_dev)+"\t")
				ECSt_low_ave = MathHelper.mean(ECSt_low_list)
				ECSt_low_dev = MathHelper.stdev(ECSt_low_list)
				combined_output.write(str(ECSt_low_ave)+"\t"+str(ECSt_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(ECSt_amb_list, ECSt_low_list))+"\t")
				# ECS amplitude
				ECSamp_amb_ave = MathHelper.mean(ECSamp_amb_list)
				ECSamp_amb_dev = MathHelper.stdev(ECSamp_amb_list)
				combined_output.write(str(ECSamp_amb_ave)+"\t"+str(ECSamp_amb_dev)+"\t")
				ECSamp_low_ave = MathHelper.mean(ECSamp_low_list)
				ECSamp_low_dev = MathHelper.stdev(ECSamp_low_list)
				combined_output.write(str(ECSamp_low_ave)+"\t"+str(ECSamp_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(ECSamp_amb_list, ECSamp_low_list))+"\t")
				# gH+
				gH_amb_ave = MathHelper.mean(gH_amb_list)
				gH_amb_dev = MathHelper.stdev(gH_amb_list)
				combined_output.write(str(gH_amb_ave)+"\t"+str(gH_amb_dev)+"\t")
				gH_low_ave = MathHelper.mean(gH_low_list)
				gH_low_dev = MathHelper.stdev(gH_low_list)
				combined_output.write(str(gH_low_ave)+"\t"+str(gH_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(gH_amb_list, gH_low_list))+"\t")
				# vH+
				vH_amb_ave = MathHelper.mean(vH_amb_list)
				vH_amb_dev = MathHelper.stdev(vH_amb_list)
				combined_output.write(str(vH_amb_ave)+"\t"+str(vH_amb_dev)+"\t")
				vH_low_ave = MathHelper.mean(vH_low_list)
				vH_low_dev = MathHelper.stdev(vH_low_list)
				combined_output.write(str(vH_low_ave)+"\t"+str(vH_low_dev)+"\t")
				combined_output.write(str(MathHelper.p_value(vH_amb_list, vH_low_list))+"\t")
				# p700
				p700t_ave = MathHelper.mean(p700t_list)
				p700t_dev = MathHelper.stdev(p700t_list)
				combined_output.write(str(p700t_ave)+"\t"+str(p700t_dev)+"\t")
				p700k_ave = MathHelper.mean(p700k_list)
				p700k_dev = MathHelper.stdev(p700k_list)
				combined_output.write(str(p700k_ave)+"\t"+str(p700k_dev)+"\t")
				p700amp_ave = MathHelper.mean(p700amp_list)
				p700amp_dev = MathHelper.stdev(p700amp_list)
				combined_output.write(str(p700amp_ave)+"\t"+str(p700amp_dev)+"\t")
				# done
				combined_output.write("\n")
			except Exception as ex:
				print("\tError: "+str(ex))
				print(traceback.format_exc())
		# close the files
		individual_output.close()
		combined_output.close()
	finally:
		done()

def process_experiment(experiment):
	print("processing %s..." % experiment["name"], end="")
	results = SpecDataProcessor.process_experiment(experiment)
	print("done!")
	return results


def get_experiments_from_folder(folder_path):
	# using the _notes.txt file to identify experiment names
	# assumes that all experiments start with the same prefix
	experiments = []
	files = FileHelper.get_files_in_folder(folder_path)
	for file in files:
		filename = str(os.path.basename(file))
		if(filename.endswith("_notes.txt")):
			# found notes file
			experiment_name = filename.replace("_notes.txt","")
			# do we have O2 data and spec data?
			regex_name = re.escape(experiment_name) # handles special regex characters that might be in the name
			has_O2_data = is_regex_in_list(regex_name+".*"+re.escape(".csv"), files)
			has_flr_data = is_regex_in_list(regex_name+".*_t.*_flr_\\d+\\.dat", files)
			if(has_O2_data and has_flr_data):
				# there's enough data to process this experiment
				exp_files = FileHelper.filter_files_by_prefix(experiment_name, files)
				experiments.append( create_experiment_from_files(experiment_name, folder_path, exp_files) )
	return experiments

			

def create_experiment_from_files(exp_name, folder_path, file_list):
	"""
	creates a python dictionary object representing the experiment with the 
	following format:
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
	}
	"ecs_data":python dictionary{
		1:t1 ECS dat file as pandas DataFrame,
		3:t3 ECS dat file as pandas DataFrame,
		... etc.
	}
	"p700_data":python dictionary{
		6:t6 p700 dat file as pandas DataFrame,
		... etc.
	}
}
	thus you could reference the t2 flr data like this:
	flr2_dataframe = data_dict["flr_data"][2]
	"""
	print("parsing "+str(exp_name)+" from folder '"+str(folder_path)+"'")
	experiment_dict = {} # create data dictionary object
	experiment_dict["folder"] = folder_path
	experiment_dict["name"] = exp_name
	# extract rep number
	rep_str, exp_name = regex_extract("rep\\d+",exp_name)
	if( rep_str == ''):
		rep_num = 0
	else:
		rep_num = int(rep_str[3:])
	experiment_dict["repetition"] = rep_num
	experiment_dict["experiment"] = exp_name
	print("\tfound experiment "+str(exp_name)+" rep. #"+str(rep_num))
	for file in file_list:
		# get name of file
		file_name = str(os.path.basename(file))
		rep_str, file_name = regex_extract("rep\\d+",file_name)
		while (file_name.startswith("_")):
			file_name = file_name[1:]
		if(file_name.endswith("_notes.txt")): # add notes
			experiment_dict["script_notes"] = FileHelper.read_text_file(file)
		if(file_name.endswith(".csv")): # add O2 data
			experiment_dict["O2_data"] = pandas.read_csv(file,sep=',',header=1)
		if(file_name.endswith(".dat")): # parse spec data file
			components = file_name.split("_") # split the name by underscores
			comp_len = len(components)
			time_point_str = components[comp_len - 3] # 3rd from entry is t#
			time_point_num = int(time_point_str[1:]) # remove the 't' and parse as integer
			if( not "timepoints" in experiment_dict): # add timepoints list if not yet present
				experiment_dict["timepoints"] = []
			if( not time_point_num in experiment_dict["timepoints"]): 
				# add timepoint if it is not yet in the list
				experiment_dict["timepoints"].append(time_point_num)
				experiment_dict["timepoints"].sort()
			#
			measure_type = components[comp_len - 2] # 2nd from end is the measurement type
			#
			keyname = measure_type+"_data" # key in experiment dictionary
			if( not keyname in experiment_dict):
				# add dict if not yet added
				experiment_dict[keyname] = {}
			#
			# read data file to data frame and store it in the experiment dictionary
			experiment_dict[keyname][time_point_num] = pandas.read_csv(file,sep='\t',header=None)
	return experiment_dict

def is_experiment_folder(folder):
	# return true if the folder contains 1 or more spec experiment files
	# and contains an oxygen probe data file
	files = FileHelper.get_files_in_folder(folder)
	has_spec_files = False
	has_O2_files = False
	for file in files:
		if is_spec_experiment_file(file):
			has_spec_files = True
		if is_O2_experiment_file(file):
			has_O2_files = True
		
	return has_spec_files and has_O2_files

def is_spec_experiment_file(filepath):
	# see https://docs.python.org/3/library/re.html for details on regex
	# see also the append_base commands in the spec script
	regex_pattern = ".*_t.*_(ecs|flr|p700)_\\d+\\.dat"
	if( re.search( regex_pattern, str(filepath) ) != None ):
		# is spec experiment file
		return True
	else:
		return False

def is_O2_experiment_file(filepath):
	return str(filepath).endswith(".csv")

def is_regex_in_list(regex_pattern, list_of_things):
	for thing in list_of_things:
		if( re.search( regex_pattern, str(thing) ) != None ):
			# found a hit
			return True
	return False


def regex_extract(regex_pattern, text):
	"""
	regex_pattern: a regular expression (aka regex)
	text: a string
	
	returns ex, t where ex is the found regular expression and t is the 
			original text without the matched expression
	"""
	# uses a regex to pull out a sub-string and return both what was pulled-
	# out and the rest of the string
	found = re.findall(regex_pattern, str(text))
	if(len(found) < 1):
		# not found in the text
		return '', text
	return found[0], str(text).replace(found[0],"")

def show_message(title, message):
	messagebox.showinfo(title, message)

def setup():
	# initialize tkinter and hide the ugly default window
	global tk_root
	tk_root = tkinter.Tk()
	tk_root.withdraw()
	
def done():
	# de-initialize
	global tk_root
	tk_root.quit()

# run main function only if this script is the one we are running
# (as opposed to calling this script file from another script)
if __name__ == "__main__":
	main()
# end of script