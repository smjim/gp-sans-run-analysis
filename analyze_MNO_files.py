import mcstasHelper as mc
import numpy as np
import pandas as pd
import os

import settings

def weighted_mean(y, yerr):
	mean = np.sum(y/ yerr**2) / np.sum(1/ yerr**2)
	mean_err = np.sqrt(1/ np.sum(1/ yerr**2))
	return mean, mean_err

def calculate_rms(data, mean):
	"""
	Calculate the root mean square (RMS) of a 1D array with error bars.
	
	Parameters:
		data (array-like): The 1D array of data points.
		errors (array-like): The 1D array of corresponding error bars.
	
	Returns:
		float: The root mean square (RMS) value.
	"""

	residuals = data - mean
	rms = np.sqrt(np.sum(np.square(residuals))/ (data.size - 1))
	
	return rms

# Calculate Sum within ROI for given ROI and input file 
def calculate_roi_sum(file, square=None, circle=None, oval=None, show=False):

	I, sigI, N, dataHeader, L = mc.extractMcStasData(file)
	
	# Find extent from file header
	extent = np.array(dataHeader['xylimits'].split(),dtype=float)
	
	# Calculate the spacing between values in the array
	dx = (extent[1] - extent[0]) / (N.shape[1] - 1)
	dy = (extent[3] - extent[2]) / (N.shape[0] - 1)
	
	# Create the mask for the ROI
	mask = np.zeros_like(N, dtype=bool)

	roi_area = 0

	if square:
		x0, x1, y0, y1 = square

		# Square ROI defined by 2 corners
		x_indices = np.where((x0 <= extent[0] + np.arange(N.shape[1]) * dx) & (extent[0] + np.arange(N.shape[1]) * dx < x1))
		y_indices = np.where((y0 <= extent[2] + np.arange(N.shape[0]) * dy) & (extent[2] + np.arange(N.shape[0]) * dy < y1))
		#y_indices = np.where((-y1 <= extent[2] + np.arange(N.shape[0]) * dy) & (extent[2] + np.arange(N.shape[0]) * dy < -y0))	# if y is flipped

		# Set the corresponding elements in the mask to True
		mask[y_indices[0][:, np.newaxis], x_indices[0]] = True

		roi_area = (x1-x0)*(y1-y0)

	if circle: 
		x0, y0, radius = circle 

		# Circle ROI defined by center (x0, y0) and radius
		x_indices, y_indices = np.meshgrid(np.arange(N.shape[1]), np.arange(N.shape[0]))
		mask = ((extent[0] + x_indices * dx - x0) ** 2 + (extent[2] + y_indices * dy - y0) ** 2 <= radius ** 2)

		roi_area = np.pi*radius**2

	if oval:
		x0, y0, sigx, sigy = oval

		# Oval ROI defined by center (x0, y0) and (sigx, sigy)
		r_x = 2*sigx
		r_y = 2*sigy
		x_indices, y_indices = np.meshgrid(np.arange(N.shape[1]), np.arange(N.shape[0]))
		mask = (((extent[0] + x_indices * dx - x0) / r_x) ** 2 + ((extent[2] + y_indices * dy - y0) / r_y) ** 2 <= 1)

		roi_area = np.pi*r_x*r_y
	
	# Apply the mask and calculate the sum within the ROI
	roi_sum = np.sum(N[mask])
	roi_sum_err = np.sqrt(np.sum(np.square(N[mask])))

	if show:
		import matplotlib.pyplot as plt
		from matplotlib.patches import Rectangle, Circle, Ellipse
		import re

		unit1 = re.findall(r"\[(.*?)\]", dataHeader['xlabel'])
		unit2 = re.findall(r"\[(.*?)\]", dataHeader['ylabel'])
	
		# Show ROI mask
		plt.imshow(mask, cmap='binary', extent=extent)
		plt.colorbar()
		plt.title('ROI')
		plt.show()
		# Show data with mask applied
		plt.imshow(mask*N, cmap='plasma', extent=extent)
		plt.colorbar()
		plt.title('Counts within ROI')
		plt.show()
	
		# Show data with mask outlined
		fig, ax = plt.subplots()
	
		img = ax.imshow(np.flipud(N), extent=extent, cmap='plasma')
		#img = ax.imshow(np.flipud(N), extent=extent, cmap='plasma', norm='log', vmin=10, vmax=5e6)
		ax.set_title(f"{dataHeader['component']}; ({dataHeader['position']})m")
		ax.set_xlabel(dataHeader['xlabel'])
		ax.set_ylabel(dataHeader['ylabel'])
		cbar = fig.colorbar(img, ax=ax)
		cbar.set_label(dataHeader['zvar']+'/ '+"{:.2e}".format(dx*dy)+' ['+unit1[0]+'*'+unit2[0]+']')
		#cbar.set_label('$n \cdot s^2$'+'/ '+"{:.2e}".format(dx*dy)+' ['+unit1[0]+'*'+unit2[0]+']')
	
		# Add patch for ROI outline on plot
		if square:
			square = Rectangle((x0, y0), (x1 - x0), (y1 - y0), fill=False, color='red', linewidth=2)
			ax.add_patch(square)
	
		if circle:
			circle = Circle((x0, y0), radius, fill=False, color='red', linewidth=2)
			ax.add_patch(circle)
	
		if oval:
			# Create an Ellipse patch for the oval
			# Define the width and height of the ellipse
			width = 2 * sigx
			height = 2 * sigy
			# Create the ellipse centered at (x0, y0) with the given width, height, and angle (0 for horizontal)
			ellipse = Ellipse((x0, y0), width, height, angle=0, fill=False, color='red', linewidth=2)
			# Add the ellipse to the plot
			ax.add_patch(ellipse)
	
		roi_info = f"Sum within ROI: {roi_sum:.2e} ± {roi_sum_err:.2e}\n"

		ax.annotate(roi_info,
				xy=(0.05, 0.95),
				xycoords='axes fraction',
				ha='left',
				va='top',
				fontsize=10,
				color='black',
				bbox=dict(facecolor='white', edgecolor='black', pad=3.5))
				#bbox=dict(facecolor='white', edgecolor='black', boxstyle='round', alpha=0.7, pad=0.5))
	
		plt.tight_layout()
		plt.show()

	
	return roi_sum, roi_sum_err

def extract_from_MNO(MNO_filename):
	"""
	Extract BM intensity data from a TSV .dat file.

	Parameters:
		MNO_filename (str): Path to the TSV .dat file.

	Returns:
		np.ndarray: 2D array containing time and counts.
		np.ndarray: 2D array containing time and |B| readback value [Gauss].
		np.ndarray: 2D array containing time and counts in GPSANS tubes. 
	"""
	# Initialize lists to store time and intensity values
	time_list_1 = []
	intensity_list = []

	time_list_2 = []
	B_list = []

	time_list_3 = []
	GPSANS_list = []

	# Open the file and read line by line
	with open(MNO_filename, 'r') as file:
		# Skip the header lines
		for _ in range(5):
			next(file)

		# Iterate through each line in the file
		for line in file:
			# Split the line by comma to extract values
			values = line.strip().split(',')
			tube = int(values[2].strip())
			time = float(values[1].strip())
			value = float(values[3].strip())

			# Check if the Detector ID is 2048 (BM countrate per 16667us)
			if tube == 2048:
				time_list_1.append(time)
				intensity_list.append(value)
			# Check if the Detector ID is 2050 (Magnetic Field Readback value in Gauss)
			if tube == 2050:
				time_list_2.append(time)
				B_list.append(value)
			# Check if the Detector ID is in [0, 191] (GP-SANS tubes)
			if 0 <= tube and tube <= 191:
				time_list_3.append(time)
				GPSANS_list.append(1)
				#GPSANS_list.append(value) # this will give the y position of the detection

	# Convert lists to numpy arrays
	time_array_1 = np.array(time_list_1)
	intensity_array = np.array(intensity_list)
	intensity_err_array = np.sqrt(intensity_array) # Statistical error from countrate

	time_array_2 = np.array(time_list_2)
	B_array = np.array(B_list)

	time_array_3 = np.array(time_list_3)
	GPSANS_array = np.array(GPSANS_list)

	# Combine arrays column-wise
	intensity_array = np.column_stack((time_array_1, intensity_array, intensity_err_array))
	B_array = np.column_stack((time_array_2, B_array))
	GPSANS_array = np.column_stack((time_array_3, GPSANS_array))

	return intensity_array, B_array, GPSANS_array

# Calculate sum within ROI for set of "selected runs" (Pandas DataFrame)
def find_ROI_counts(selected_runs, ROI_list, data_dir):
	"""
	Function to select runs from the run summary DataFrame based on given |B| values.

	Parameters:
		selected_runs (DataFrame): DataFrame containing selected run summary data.
		ROI_list (array): Contains [Rectangular ROI limits [x0, x1, y0, y1]] [DataFrame column header for ROI sum output] [DataFrame column header for ROI sum err output]
		data_dir (string): Name of directory containing McStas-format data files to be analyzed.

	Returns:
		selected_runs (DataFrame): DataFram now has extra column containing calculated ROI counts.
	"""

	selected_runs_appended = selected_runs.copy()

	selected_files = []

	# Iterate through selected runs to find corresponding MNO files
	for run_number in selected_runs['Run #']:
		# Construct file path for each run
		#file_path = os.path.join(data_dir, f'{run_number}_output.dat')
		file_path = os.path.join(data_dir, f'reduced_xbin_{run_number}_output.dat')
		# Check if file exists
		if os.path.exists(file_path):
			selected_files.append(file_path)
		else:
			print(f"File not found for run {run_number}: {file_path}")
	
	print(ROI_list)
	for (ROI, roi_string, roi_err_string) in ROI_list:
	
		print(f' ======== Analyzing {roi_string} ========')
		# Find counts within rectangular ROI limits for each MNO file
		for file in selected_files:
			ROI_sum, ROI_sum_err = calculate_roi_sum(file, square=ROI, show=False)
			# Extract run number from file name
			run_number = os.path.basename(file).split('_')[2]
			# Update selected_runs_appended with ROI sum for corresponding run
			selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), roi_string] = ROI_sum
			selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), roi_err_string] = ROI_sum_err
			print(f'{run_number} {roi_string}: {ROI_sum}, {ROI_sum_err}')
		print(f' ======== {roi_string} Done. ========')
	
	return selected_runs_appended

# Calculate sum within ROI for set of "selected runs" (Pandas DataFrame) Using MNO datafiles instead of preprocessed files
def find_ROI_counts_MNO(selected_runs, ROI_list, data_dir):
	"""
	Function to find counts in different ROI's for all selected runs

	Parameters:
		selected_runs (DataFrame): DataFrame containing selected run summary data.
		ROI_list (array): Contains [Rectangular ROI limits [x0, x1, y0, y1]] [DataFrame column header for ROI sum output] [DataFrame column header for ROI sum err output]
		data_dir (string): Name of directory containing McStas-format data files to be analyzed.

	Returns:
		selected_runs (DataFrame): DataFram now has extra column containing calculated ROI counts.
	"""

	selected_runs_appended = selected_runs.copy()

	selected_files = []

	# Iterate through selected runs to find corresponding MNO files
	for run_number in selected_runs['Run #']:
		# Construct file path for each run
		file_path = os.path.join(settings.mno_dir, f'MNO_GPSANS_{run_number}.txt')
		# Check if file exists
		if os.path.exists(file_path):
			selected_files.append(file_path)
		else:
			print(f"File not found for run {run_number}: {file_path}")
	
	print(ROI_list)
	for (ROI, roi_string, roi_err_string) in ROI_list:
	
		print(f' ======== Analyzing {roi_string} ========')
		# Find counts within rectangular ROI limits for each MNO file
		for file in selected_files:
			# Extract run number from file name
			run_number = (os.path.basename(file).split('_')[2]).split('.')[0]

			processed_outfile = os.path.join(settings.processed_outdir, f'reduced_xbin_{run_number}.txt')
			if not os.path.exists(processed_outfile):
				# Preprocess MNO datafiles
				run_duration = os.system(f"python3 preprocess/extract_duration.py {settings.run_summary_file} {run_number}")
				detIDmap_file = os.path.join(settings.data_dir, 'DetIDmap.csv')

				os.system(f"./MNO_to_intensity.out {file} {detIDmap_file} {processed_outfile} {run_duration}")
			
			ROI_sum, ROI_sum_err = calculate_roi_sum(processed_outfile, square=ROI, show=False)

			# Update selected_runs_appended with ROI sum for corresponding run
			selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), roi_string] = ROI_sum
			selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), roi_err_string] = ROI_sum_err
			print(f'{run_number} {roi_string}: {ROI_sum}, {ROI_sum_err}')
		print(f' ======== {roi_string} Done. ========')
	
	return selected_runs_appended

# Calculate average BM intensity 
# MNO tube 2048 gives counts/ 16667us
def find_BM_intensity_MNO(selected_runs):
	"""
	Function to find BM average intensity for all selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing selected run summary data.
		data_dir (string): Name of directory containing McStas-format data files to be analyzed.

	Returns:
		selected_runs (DataFrame): DataFram now has extra column containing calculated average BM intensity. 
	"""

	selected_runs_appended = selected_runs.copy()

	selected_files = []

	# Iterate through selected runs to find corresponding MNO files
	for run_number in selected_runs['Run #']:
		# Construct file path for each run
		file_path = os.path.join(settings.mno_dir, f'MNO_GPSANS_{run_number}.txt')
		# Check if file exists
		if os.path.exists(file_path):
			selected_files.append(file_path)
		else:
			print(f"File not found for run {run_number}: {file_path}")
	
	print(f' ======== Analyzing BM Intensity ========')
	for file in selected_files:
		# Extract run number from file name
		run_number = (os.path.basename(file).split('_')[2]).split('.')[0]
		
		# -------------
		# BM Intensity 
		# -------------
		BM_intensity, B_value, GPSANS_counts = extract_from_MNO(file)	# Intensity per 16667us 
		BM_intensity, BM_intensity_err = BM_intensity[:,1], np.sqrt(BM_intensity[:,1]) 

		BM_intensity_avg = np.mean(BM_intensity)								# Average Intensity [counts/ tick]
		BM_intensity_avg_err = calculate_rms(BM_intensity, BM_intensity_avg)	# Err on Average Intensity [counts/ tick]

		# Update selected_runs_appended with bm_intensity_avg for corresponding run
		selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), 'BM Average Intensity [per tick]'] = BM_intensity_avg 
		selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), 'BM Average Intensity Err [per tick]'] = BM_intensity_avg_err 

		print(f'{run_number} BM Intensity Avg: {BM_intensity_avg}')
	print(f' ======== BM Intensity Analysis Done. ========')
	
	return selected_runs_appended

# Calculate sum within ROI's and BM intensity average using MNO file 
def analyze_MNO(selected_runs, ROI_list, data_dir, output_summary_file):
	"""
	Function to find counts in different ROI's for all selected runs

	Parameters:
		selected_runs (DataFrame): DataFrame containing selected run summary data.
		ROI_list (array): Contains [Rectangular ROI limits [x0, x1, y0, y1]] [DataFrame column header for ROI sum output] [DataFrame column header for ROI sum err output]
		data_dir (string): Name of directory containing McStas-format data files to be analyzed.
		output_summary_file (string): Name of file containing output summary of analysis.

	Returns:
		selected_runs (DataFrame): DataFram now has extra column containing calculated ROI counts, BM intensity average.
	"""

	# If output_summary_file exists and has all the required data, load it instead of redoing analysis
	if os.path.exists(output_summary_file):
		# Load output_summary_file into a pandas DataFrame
		output_summary = pd.read_csv(output_summary_file)
		load_summary = True
	
		for (ROI, roi_string, roi_err_string) in ROI_list:
			# Check if output_summary has column roi_string
			if roi_string not in output_summary.columns:
				print(f"Column '{roi_string}' not found in output_summary.")
				load_summary = False
	
			# Check if output_summary has column roi_err_string
			if roi_err_string not in output_summary.columns:
				print(f"Column '{roi_err_string}' not found in output_summary.")
				load_summary = False
	
		# Check if output_summary has column 'BM Average Intensity [per tick]'
		if 'BM Average Intensity [per tick]' not in output_summary.columns:
			print("Column 'BM Average Intensity [per tick]' not found in output_summary.")
			load_summary = False
		if 'BM Average Intensity Err [per tick]' not in output_summary.columns:
			print("Column 'BM Average Intensity Err [per tick]' not found in output_summary.")
			load_summary = False
	
		runs = selected_runs['Run #']
		# Check if output_summary['Run #'] has every run that can be found in runs
		if not output_summary['Run #'].isin(runs).all():
			print("Not all runs from selected_runs are present in output_summary['Run #'].")
			load_summary = False

		# If all checks pass
		if load_summary == True:
			print(f'Loading data from {output_summary_file}')
			return output_summary
	
	# If any of the conditions above were false:
	# Redo the analysis

	selected_runs_appended = selected_runs.copy()

	selected_files = []

	# Iterate through selected runs to find corresponding MNO files
	for run_number in selected_runs['Run #']:
		# Construct file path for each run
		file_path = os.path.join(settings.mno_dir, f'MNO_GPSANS_{run_number}.txt')
		# Check if file exists
		if os.path.exists(file_path):
			selected_files.append(file_path)
		else:
			print(f"File not found for run {run_number}: {file_path}")

	print(f' ======== Converting MNO Files ========')

	for file in selected_files:
		# Extract run number from file name
		run_number = (os.path.basename(file).split('_')[2]).split('.')[0]

		processed_outfile = os.path.join(settings.processed_outdir, f'reduced_xbin_{run_number}.txt')

		if not os.path.exists(processed_outfile):
			# Preprocess MNO datafiles
			run_duration = os.system(f"python3 preprocess/extract_duration.py {settings.run_summary_file} {run_number}")
			detIDmap_file = os.path.join(settings.data_dir, 'DetIDmap.csv')

			os.system(f"./MNO_to_intensity.out {file} {detIDmap_file} {processed_outfile} {run_duration}")

	print(f' ======== MNO Conversion Done ========')
	
	print(ROI_list)
	
	print(f' ======== Analyzing Runs ========')

	for file in selected_files:

		# -------------
		# ROI Sums
		# -------------
		roi_count_str = '	 '

		# Find counts within rectangular ROI limits for each MNO file
		for (ROI, roi_string, roi_err_string) in ROI_list:
			# Extract run number from file name
			run_number = (os.path.basename(file).split('_')[2]).split('.')[0]
			processed_outfile = os.path.join(settings.processed_outdir, f'reduced_xbin_{run_number}.txt')
			
			ROI_sum, ROI_sum_err = calculate_roi_sum(processed_outfile, square=ROI, show=False)

			# Update selected_runs_appended with ROI sum for corresponding run
			selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), roi_string] = ROI_sum
			selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), roi_err_string] = ROI_sum_err

			roi_count_str = f'{roi_count_str}{roi_string} Total Counts: {int(ROI_sum)} ± {ROI_sum_err:.2f}\n	 '

		# -------------
		# BM Intensity 
		# -------------
		BM_intensity, B_value, GPSANS_counts = extract_from_MNO(file)	# Intensity per 16667us 
		BM_intensity, BM_intensity_err = BM_intensity[:,1], np.sqrt(BM_intensity[:,1]) 

		BM_intensity_avg, BM_intensity_avg_err = weighted_mean(BM_intensity, BM_intensity_err)	# Average Intensity [counts/ tick]
		#BM_intensity_avg = np.mean(BM_intensity)								# Average Intensity [counts/ tick]
		#BM_intensity_avg_err = calculate_rms(BM_intensity, BM_intensity_avg)	# Err on Average Intensity [counts/ tick]

		# Update selected_runs_appended with bm_intensity_avg for corresponding run
		selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), 'BM Average Intensity [per tick]'] = BM_intensity_avg 
		selected_runs_appended.loc[selected_runs_appended['Run #'] == int(run_number), 'BM Average Intensity Err [per tick]'] = BM_intensity_avg_err 

		print(f'{run_number} BM Intensity Avg: {BM_intensity_avg:.3e} ± {BM_intensity_avg_err:.3e} Counts/ Hour\n{roi_count_str}')

	print(f' ======== Run Analysis Done. ========')

	# Save selected_runs_appended to output_summary_file
	selected_runs_appended.to_csv(output_summary_file)

	return selected_runs_appended
