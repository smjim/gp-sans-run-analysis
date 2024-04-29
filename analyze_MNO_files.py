import mcstasHelper as mc
import numpy as np
import pandas as pd
import os

import settings

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
	Function to select runs from the run summary DataFrame based on given |B| values.

	Parameters:
		selected_runs (DataFrame): DataFrame containing selected run summary data.
		ROI_list (array): Contains [Rectangular ROI limits [x0, x1, y0, y1]] [DataFrame column header for ROI sum output] [DataFrame column header for ROI sum err output]
		data_dir (string): Name of directory containing McStas-format data files to be analyzed.

	Returns:
		selected_runs (DataFrame): DataFram now has extra column containing calculated ROI counts.
	"""

	mno_dir = os.path.join(settings.data_dir, 'MNO_files')
	processed_outdir = os.path.join(settings.data_dir, 'processed_data')

	if not os.path.exists(mno_dir):
		print('Please help us find your MNO files!')
	if not os.path.exists(processed_outdir):
		os.makedirs(processed_outdir)

	selected_runs_appended = selected_runs.copy()

	selected_files = []

	# Iterate through selected runs to find corresponding MNO files
	for run_number in selected_runs['Run #']:
		# Construct file path for each run
		file_path = os.path.join(mno_dir, f'MNO_GPSANS_{run_number}.txt')
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

			processed_outfile = os.path.join(processed_outdir, f'reduced_xbin_{run_number}.txt')
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
