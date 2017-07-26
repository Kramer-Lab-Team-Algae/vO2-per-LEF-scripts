import os


def get_files_in_folder(folder):
	"""
	returns a list of just the files in a folder
	"""
	only_files = []
	for filename in os.listdir(folder):
		full_path = os.path.join(folder, filename)
		if os.path.isfile(full_path):
			only_files.append(full_path)
	return only_files

def read_text_file(file_path):
	"""
	reads the entire contents of a text file and returns it as a string
	"""
	with open(file_path, 'r') as file_in:
		content = file_in.read()
	return content

def filter_files_by_prefix(prefix, file_list):
	output = []
	for f in file_list:
		if( str(os.path.basename(f)).startswith(prefix) ):
			output.append(f)
	return output

def filter_files_by_suffix(suffix, file_list):
	output = []
	for f in file_list:
		if( str(os.path.basename(f)).endswith(suffix) ):
			output.append(f)
	return output