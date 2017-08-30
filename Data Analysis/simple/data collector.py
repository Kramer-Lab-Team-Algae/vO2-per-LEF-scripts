
import os
import tkinter
from tkinter import filedialog
from tkinter import messagebox

def main():
	#root_folder = os.path.abspath('.') # current folder (might be different from script file location)
	#root_folder = os.path.dirname(os.path.abspath(__file__)) # folder with this script file in it
	root_folder = filedialog.askdirectory() # ask user for location
	
	print('Scanning folder: ',root_folder)
	#
	output_files = []
	for folder_path, directories_list, file_list in os.walk(root_folder):
		for file_name in file_list:
			try:
				file_path = str(os.path.join(folder_path, file_name))
				if(file_path.endswith('_out.txt')):
					output_files.append(file_path)
			except Exception as ex:
				print("\tError: "+str(ex))
				print(traceback.format_exc())
	#
	total_content = ""
	for file_path in output_files:
		total_content += file_path
		total_content += '\n'
		with open(file_path,'r') as f_in:
			total_content += f_in.read()
	print(total_content)
	with filedialog.asksaveasfile(mode='w', defaultextension=".txt") as f_out:
		if(f_out != None):
			f_out.write(total_content)
	#
	exit()

if __name__ == "__main__":
	#
	main()
