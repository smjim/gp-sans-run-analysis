import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

import settings

from analyze_MNO_files import *

# Colors
dark_blue_hex = '#5928ED'
blue_hex = '#0073E6'
gray_hex = '#111111'

# Font settings
fontsize = 14

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : fontsize}

mpl.rc('font', **font)

def plot_fit(ax, x, m, b, chisq, ndof):
	fit_label = fr'$y = {m:.2e} \cdot x + {b:.2e}$' + '\n' + r'$ \frac{\chi^2}{nDof} = ' + f'{chisq/ndof:.2f}$'
	ax.plot(x, m * x + b, '--', label=fit_label, color='black')

def calculate_chi_square(x, y, y_err, m, b):
	y_fit = m * x + b
	chi_sq = np.sum(((y - y_fit) / y_err) ** 2)  # Calculate chi-square value
	ndof = y.size - 2
	return chi_sq, ndof

def rebin_histogram(hist, new_size):
	old_size = hist.size
	
	if new_size >= old_size:
		raise ValueError("The new size should be smaller than the old size.")
	
	if old_size % new_size != 0:
		rem = old_size % new_size
		print("! New size is not a factor of old size, cutting out last ", rem, " entries")
		hist = hist[:-rem]
	
	factor = old_size // new_size
	
	# Reshape the histogram array to a 2D array
	hist_2d = hist.reshape((new_size, factor))
	
	# Sum the values in each bin of the new histogram
	rebinned_hist = np.sum(hist_2d, axis=1)
	
	return rebinned_hist

# Rebin data_array of the column format [time, value, value_err]
def rebin_run_data(data_array, bin_width):
	if bin_width <=0:	# If bin_width <= 0, dont do the rebinning
		return data_array, None
	
	duration = np.max(data_array[:,0])
	new_size = int(duration/ bin_width)

	rebinned_time = rebin_histogram(data_array[:-1,0], new_size)/ (data_array[:-1,0].size // new_size)	# Average col-wise
	bin_width_t = (data_array[-1,0] - data_array[0,0])/ new_size 

	# sum col wise gives (counts/ tick * ticks/ hour) * bin_width[t] = counts/hour * binwidth
	# divide by bin_width_t to renormalize
	rebinned_data = rebin_histogram(data_array[:,1], new_size)							# Sum col-wise
	rebinned_data_err = np.sqrt(rebin_histogram(np.square(data_array[:,2]), new_size))	# sqrt of sum col-wise of squared values

	# Combine to create rebinned data array
	rebinned_data_array = np.column_stack((rebinned_time, rebinned_data, rebinned_data_err))
	
	return rebinned_data_array, bin_width_t

# rebin the GPSANS histogram
# TODO fix this so that it is bin_width instead of num_bins
def event_time_histogram(event_times, bin_width):
	"""
	Calculate a histogram of a 1D array of event times.

	Parameters:
		event_times (array-like): The 1D array listing event times.
		bin_width (int): The width (t) of bins to use for the histogram. 

	Returns:
		numpy.ndarray: A 2-column array with the columns being "bin_t_midpoint" and "counts_per_bin".
	"""
	
	if bin_width <= 0:
		bin_width = 1
		print('binning GP-SANS event time histogram to 1s')

	# Calculate the range of the event times
	min_time = np.min(event_times)
	max_time = np.max(event_times)

	# Calculate the number of bins based on the bin width
	num_bins = int(np.ceil((max_time - min_time) / bin_width))

	# Create an array to store the counts in each bin
	counts_per_bin = np.zeros(num_bins)

	# Calculate the bin edges
	bin_edges = np.linspace(min_time, min_time + num_bins * bin_width, num_bins + 1)

	# Calculate the bin midpoints
	bin_t_midpoint = (bin_edges[:-1] + bin_edges[1:]) / 2

	# Iterate through event times and count them into bins
	for event_time in event_times:
		bin_index = int((event_time - min_time) // bin_width)
		if bin_index >= num_bins:
			bin_index = num_bins - 1
		counts_per_bin[bin_index] += 1

	return np.column_stack((bin_t_midpoint, counts_per_bin))

def plot_counts_vs_time(selected_runs, roi_list):
	"""
	Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
		roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
	"""

	# Establish colorscheme
	cmap = mpl.colormaps['RdYlGn']

	for roi, roi_str, roi_err_str in roi_list:
		plt.figure(figsize=(10, 6))
		for i, B_value in enumerate(settings.B_select):
			runs_for_B = selected_runs[selected_runs['Field'] == B_value]
			duration = runs_for_B['Duration [Hours]']
			midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

			counts, counts_err = runs_for_B[roi_str], runs_for_B[roi_err_str]

			color = cmap(i / len(settings.B_select))  # Sample color from colormap
			plt.errorbar(midpoint_time, counts, xerr=duration/ 2, yerr=counts_err,
					marker='o', markersize=3, capsize=2, fmt=' ', color=color, label=f'|B| = {B_value}')

		plt.xlabel('Start Time [Hrs]', weight='bold')
		plt.ylabel(f'{roi_str} Total Counts', weight='bold')
		plt.title('Counts vs. Start Time', weight='bold')
		plt.legend(loc='lower right', fontsize=fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_counts_plot.pdf', format='pdf')
		plt.show()

def plot_countrates_vs_time(selected_runs, roi_list, fit=False):
	"""
	Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
		roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
	"""

	# Establish colorscheme
	cmap = mpl.colormaps['RdYlGn']

	for roi, roi_str, roi_err_str in roi_list:
		fig, ax = plt.subplots(figsize=(10, 6))

		# Calculate line of best fit
		if fit:
			duration = selected_runs['Duration [Hours]']
			midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

			countrates, countrates_err = selected_runs[roi_str]/ duration, selected_runs[roi_err_str]/ duration

			m, b = np.polyfit(midpoint_time, countrates, 1)
			chisq, ndof = calculate_chi_square(midpoint_time, countrates, countrates_err, m, b) 

			plot_fit(ax, midpoint_time, m, b, chisq, ndof)	

		for i, B_value in enumerate(settings.B_select):
			runs_for_B = selected_runs[selected_runs['Field'] == B_value]
			duration = runs_for_B['Duration [Hours]']
			midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

			countrates, countrates_err = runs_for_B[roi_str]/ duration, runs_for_B[roi_err_str]/ duration

			color = cmap(i / len(settings.B_select))  # Sample color from colormap
			ax.errorbar(midpoint_time, countrates, xerr=duration/ 2, yerr=countrates_err,
					marker='o', markersize=3, capsize=2, fmt=' ', color=color, label=f'|B| = {B_value}')

		plt.xlabel('Start Time [Hrs]', weight='bold')
		plt.ylabel(roi_str + ' Count Rate [$hours^{-1}$]', weight='bold')
		plt.title('Count Rate vs. Start Time', weight='bold')
		plt.legend(loc='lower right', fontsize=fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_countRate_plot.pdf', format='pdf')
		plt.show()

def plot_BM_intensity_vs_time(selected_runs, hl=None):
	"""
	Function to plot the `Beam Monitor` (ID 2048) Average value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
	"""

	plt.figure(figsize=(10, 6))

	duration = selected_runs['Duration [Hours]']
	midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

	BM_avg_intensity, BM_avg_intensity_err = selected_runs['BM Average Intensity [per tick]'], selected_runs['BM Average Intensity Err [per tick]']

	plt.errorbar(midpoint_time, BM_avg_intensity*60, xerr=duration/ 2, yerr=BM_avg_intensity_err*60,
			fmt=' ', marker='o', markersize=5, capsize=2, color='black', label=f'BM Average Intensity for Selected Runs')

	if hl is not None:
		# Convert 'Run #' to string in selected_runs for consistency
		selected_runs['Run #'] = selected_runs['Run #'].astype(str)
		for run_num in hl:
			run_num_str = str(run_num).strip()
			# Find the index where 'Run #' equals run_num_str
			index = selected_runs[selected_runs['Run #'] == run_num_str].index
			if not index.empty:
				plt.errorbar(midpoint_time[index], BM_avg_intensity[index] * 60, xerr=duration[index]/ 2, yerr=BM_avg_intensity_err[index] * 60,
						fmt=' ', marker='o', markersize=10, capsize=2, alpha=0.6, label=f'Run #{run_num_str}')

	plt.xlabel('Start Time [Hrs]', weight='bold')
	plt.ylabel('Weighted Mean Count Rate [$Bq$]', weight='bold')
	#plt.ylabel('Average Count Rate [$(16667 \mu s)^{-1}$]')
	plt.title('Average BM Count Rate vs. Start Time for Selected Runs', weight='bold')
	plt.legend(fontsize=fontsize*0.8)

	plt.tight_layout()
	plt.savefig(f'{settings.output_dir}BM_avg_intensity_plot.pdf', format='pdf')
	plt.show()

# Plots countrates for selected runs vs midpoint time elapsed of run, with |B| coloring, normalized to BM average intensity
def plot_countrates_vs_time_normalized(selected_runs, roi_list, fit=False, showBM=True):
	"""
	Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
		roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
	"""

	# Establish colorscheme
	cmap = mpl.colormaps['RdYlGn']

	for roi, roi_str, roi_err_str in roi_list:

		width = 10
		if fit: 
			width += 2
		if showBM:
			fig, (ax_bm, ax_counts) = plt.subplots(2, 1, figsize=(width, 8), sharex=True)
		else:
			fig, ax_counts = plt.subplots(figsize=(width, 6))

		# Calculate line of best fit
		if fit:
			duration = selected_runs['Duration [Hours]']
			midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

			BM_avg_intensity, BM_avg_intensity_err = selected_runs['BM Average Intensity [per tick]'], selected_runs['BM Average Intensity Err [per tick]']
			countrates, countrates_err = selected_runs[roi_str]/ duration, selected_runs[roi_err_str]/ duration

			ratio = countrates/ BM_avg_intensity
			ratio_err = np.sqrt((countrates_err/BM_avg_intensity)**2 + (ratio)**2 * (BM_avg_intensity_err/BM_avg_intensity)**2) 

			m, b = np.polyfit(midpoint_time, ratio, 1)
			chisq, ndof = calculate_chi_square(midpoint_time, ratio, ratio_err, m, b) 

			plot_fit(ax_counts, midpoint_time, m, b, chisq, ndof)	

		for i, B_value in enumerate(settings.B_select):
			runs_for_B = selected_runs[selected_runs['Field'] == B_value]
			duration = runs_for_B['Duration [Hours]']
			midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

			# Plot Normalized countrate on top figure
			BM_avg_intensity, BM_avg_intensity_err = runs_for_B['BM Average Intensity [per tick]'], runs_for_B['BM Average Intensity Err [per tick]']
			countrates, countrates_err = runs_for_B[roi_str]/ duration, runs_for_B[roi_err_str]/ duration

			ratio = countrates/ BM_avg_intensity
			ratio_err = np.sqrt((countrates_err/BM_avg_intensity)**2 + (ratio)**2 * (BM_avg_intensity_err/BM_avg_intensity)**2) 

			color = cmap(i / len(settings.B_select))  # Sample color from colormap

			ax_counts.errorbar(midpoint_time, ratio, xerr=duration/ 2, yerr=ratio_err,
					marker='o', markersize=5, capsize=2, fmt=' ', color=color, label=f'|B| = {B_value}')

			# Plot BM Intensity on bottom figure
			if showBM:
				if i == 0:
					bm_label='BM average intensity for runs'
				else:
					bm_label='_nolegend_'
				ax_bm.errorbar(midpoint_time, BM_avg_intensity, xerr=duration/ 2, yerr=BM_avg_intensity_err,
						fmt=' ', marker='o', markersize=5, capsize=2, color='black', label=bm_label)

		plt.xlabel('Start Time [Hrs]', weight='bold')
		ax_counts.set_ylabel('Count Rate/ BM Count Rate', weight='bold')
		if showBM:
			ax_bm.set_ylabel('BM Average Count Rate [$(16667 \mu s)^{-1}$]', weight='bold')
			ax_bm.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0,
				fancybox=True, shadow=True)

		plt.suptitle(f'BM Normalized Count Rate vs. Start Time, {roi_str}', weight='bold')
		ax_counts.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0,
			fancybox=True, shadow=True)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_normalized_countRate_plot.pdf', format='pdf')
		plt.show()


# Plot summary of ROI counts for each |B| value over entire experiment (x: |B|, y: counts, color: ROI def)
def plot_effect_analysis(selected_runs, ROI_list):
	"""
	Function to summarize the observed effect for different |B| values, compared by ROI def 

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
		roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
	"""

	# Establish colorscheme
	cmap = mpl.colormaps['tab10']

	for roi, roi_str, roi_err_str in ROI_list:
		plt.figure(figsize=(10, 6))


		B_list = []
		intensity_list = []

		# For each roi, calculated summary of effect for each |B| value
		for i, B_value in enumerate(settings.B_select):
			runs_for_B = selected_runs[selected_runs['Field'] == B_value]

			# Calculate sum of intensity within ROI for all runs with given |B| value
			intensity_sum = runs_for_B[roi_str].sum()  # Sum of intensity within ROI

			B_list.append(B_value)
			intensity_list.append(intensity_sum)

		intensity_array = np.array(intensity_list)
		B_array = np.array(B_list)

		B_intensity_array = np.column_stack((B_array, intensity_array))

		plt.scatter(B_intensity_array[:,0], B_intensity_array[:,1], label=f'{roi_str}')

		# TODO: plot each roi with a different color, but plot (|B| on x) and (summary fom on y)

#			runs_for_B = selected_runs[selected_runs['Field'] == B_value]
#			duration = runs_for_B['Duration [Hours]']
#			midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2
#
#			counts, counts_err = runs_for_B[roi_str], runs_for_B[roi_err_str]
#
#			color = cmap(i / len(settings.B_select))  # Sample color from colormap
#			plt.errorbar(midpoint_time, counts, xerr=duration/ 2, yerr=counts_err,
#					marker='o', markersize=3, capsize=2, fmt=' ', color=color, label=f'|B| = {B_value}')

#	plt.xlim(0)

	plt.xlabel('|B| Value [T]', weight='bold')
	plt.ylabel('Total Counts for Each ROI', weight='bold')
#	plt.title('|B| effect', weight='bold')
	plt.legend(loc='lower right', fontsize=fontsize*0.8)

	plt.tight_layout()
	plt.savefig(f'{settings.output_dir}B_effect_plot.pdf', format='pdf')
	plt.show()

def plot_run_stats(run_string, rebin=0, showB=False, showGPSANS=False):
	"""
	Function to plot the `Beam Monitor` (ID 2048) values as a function of time for a selected run.

	Parameters:
		run_string (string): Selected run for plotting
		rebin (int): Bin width (seconds) to rebin data to. Leave 0 for no rebinning.
		showB (bool): Show |B| value readback (tube 2048)
		showGPSANS (bool): Show counts in GPSANS monitor.
	"""

	# Get file_path from run number
	file_path = os.path.join(settings.mno_dir, f'MNO_GPSANS_{run_string}.txt')

	# Load data from file_path using functions in extract_from_MNO 
	BM_intensity, B_value, GPSANS_counts = extract_from_MNO(file_path)

	# Rebin counts
	print(BM_intensity)
	print(GPSANS_counts)
	rebinned_BM_intensity, new_binwidth = rebin_run_data(BM_intensity, rebin)
	rebinned_GPSANS = rebinned_BM_intensity #event_time_histogram(GPSANS_counts[:,0], rebin) # TODO fix this 
	#rebinned_GPSANS = event_time_histogram(GPSANS_counts[:,0], rebin) 

	# Show additional plots conditionally
	if showB and not showGPSANS: 
		fig, (ax_counts, ax_b) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
	elif showGPSANS and not showB:
		fig, (ax_counts, ax_gpsans) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
	elif showB and showGPSANS:
		fig, (ax_counts, ax_b, ax_gpsans) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
	else:
		fig, ax_counts = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
	
	# Plot BM Counts
	#ax_counts.errorbar(rebinned_BM_intensity[:,0], rebinned_BM_intensity[:,1], yerr=rebinned_BM_intensity[:,2],
	ax_counts.errorbar(rebinned_BM_intensity[:,0], rebinned_BM_intensity[:,1], yerr=np.sqrt(rebinned_BM_intensity[:,1]),
			alpha=0.8, marker='o', markersize=3, capsize=2, fmt=' ', color=gray_hex, label=f'BM Intensity for Run {run_string}')
	ax_counts.set_ylabel(f'BM Count Rate per {rebin:.3f} s', weight='bold')
	# Plot |B| field
	if showB:
		ax_b.scatter(B_value[:,0], B_value[:,1]/1e4,
				alpha=0.8, color='black', label=f'|B| Value for Run {run_string}')
		ax_b.set_ylabel('|B| Readback Value [T]', weight='bold')
	# Plot GPSANS counts
	if showGPSANS:
		ax_gpsans.errorbar(rebinned_GPSANS[:,0], rebinned_GPSANS[:,1], yerr=np.sqrt(rebinned_GPSANS[:,1]),
				alpha=0.8, marker='o', markersize=3, capsize=2, fmt=' ', color=gray_hex, label=f'GPSANS Counts for Run {run_string}')
		ax_gpsans.set_ylabel(f'GPSANS Counts per {rebin:.3f} s', weight='bold')
	
	# Plot BM mean and mean error
	mean_BM = np.mean(BM_intensity[:,1])
	mean_BM_err = calculate_rms(BM_intensity[:,1], mean_BM) # rms(BM intensity, BM intensity err)

	#mean_rebinned_BM = mean_BM * t_binWidth_BM/ (1/60) 
	mean_rebinned_BM = np.mean(rebinned_BM_intensity[:,1])
	mean_rebinned_BM_err = calculate_rms(rebinned_BM_intensity[:,1], mean_rebinned_BM) 

	ax_counts.axhline(mean_rebinned_BM,
			linewidth=3, linestyle=':', label=f'mean = {mean_BM:.2f}*60 cps', color='black')
	ax_counts.axhspan(mean_rebinned_BM - mean_rebinned_BM_err, mean_rebinned_BM + mean_rebinned_BM_err,
			alpha=0.9, label=f'mean error = {mean_BM_err:.2f}*60 cps', facecolor=blue_hex)

	plt.xlabel('Time [s]', weight='bold')
	plt.suptitle(f'Statistics of Run {run_string}', weight='bold')
	ax_counts.legend(loc='lower right', fontsize=fontsize*0.8)

	plt.tight_layout()
	plt.savefig(f'{settings.output_dir}run_{run_string}_analysis_plot.pdf', format='pdf')
	plt.show()
