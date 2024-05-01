import matplotlib.pyplot as plt
import numpy as np
import matplotlib

import settings

from analyze_MNO_files import *

def plot_fit(ax, x, m, b, chisq):
	fit_label = fr'$y = {m:.2e} \cdot x + {b:.2e}$'+'\n'+fr'$ \chi^2 = {chisq:.2f}$'
	ax.plot(x, m * x + b, '--', label=fit_label, color='black')

def calculate_chi_square(x, y, y_err, m, b):
	y_fit = m * x + b
	chi_sq = np.sum(((y - y_fit) / y_err) ** 2)  # Calculate chi-square value
	return chi_sq


def plot_counts_vs_time(selected_runs, roi_list):
	"""
	Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
		roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
	"""

	# Establish colorscheme
	cmap = matplotlib.colormaps['RdYlGn']

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

		plt.xlabel('Start Time [Hrs]')
		plt.ylabel(f'{roi_str} Total Counts')
		plt.title('Counts vs. Start Time')
		plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

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
	cmap = matplotlib.colormaps['RdYlGn']

	for roi, roi_str, roi_err_str in roi_list:
		fig, ax = plt.subplots(figsize=(10, 6))

		# Calculate line of best fit
		if fit:
			duration = selected_runs['Duration [Hours]']
			midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

			countrates, countrates_err = selected_runs[roi_str]/ duration, selected_runs[roi_err_str]/ duration

			m, b = np.polyfit(midpoint_time, countrates, 1)
			chisq = calculate_chi_square(midpoint_time, countrates, countrates_err, m, b) 

			plot_fit(ax, midpoint_time, m, b, chisq)	

		for i, B_value in enumerate(settings.B_select):
			runs_for_B = selected_runs[selected_runs['Field'] == B_value]
			duration = runs_for_B['Duration [Hours]']
			midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

			countrates, countrates_err = runs_for_B[roi_str]/ duration, runs_for_B[roi_err_str]/ duration

			color = cmap(i / len(settings.B_select))  # Sample color from colormap
			ax.errorbar(midpoint_time, countrates, xerr=duration/ 2, yerr=countrates_err,
					marker='o', markersize=3, capsize=2, fmt=' ', color=color, label=f'|B| = {B_value}')

		plt.xlabel('Start Time [Hrs]')
		plt.ylabel(roi_str + ' Count Rate [$hours^{-1}$]')
		plt.title('Count Rate vs. Start Time')
		plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_countRate_plot.pdf', format='pdf')
		plt.show()

def plot_run_stats(run_string, rebin=1000, showB=False):
	"""
	Function to plot the `Beam Monitor` (ID 2048) values as a function of time for a selected run.

	Parameters:
		run_string (string): Selected run for plotting
	"""

	# Get file_path from run number
	file_path = os.path.join(settings.mno_dir, f'MNO_GPSANS_{run_string}.txt')

	# Load data from file_path using functions in analyze_MNO
	BM_intensity, B_value = B_BM_intensity_MNO(file_path)

	# Rebin and propagate error
	if rebin is not False:
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

		# Rebin counts
		rebinned_counts = rebin_histogram(BM_intensity[:,1], rebin)											# Sum col-wise
		rebinned_time = rebin_histogram(BM_intensity[:-1,0], rebin)/ (BM_intensity[:-1,0].size // rebin)	# Average col-wise
		rebinned_counts_err = np.sqrt(rebin_histogram(np.square(np.sqrt(BM_intensity[:,1])), rebin))		# sqrt of sum col-wise of squared values

		# Rebin |B| Readback 
		#rebinned_counts = rebin_histogram(B_value[:,1], rebin)									# Sum col-wise
		#rebinned_time = rebin_histogram(B_value[:-1,0], rebin)/ (B_value[:-1,0].size // rebin)	# Average col-wise
		#rebinned_counts_err = np.sqrt(rebin_histogram(np.square(np.sqrt(B_value[:,1])), rebin))	# sqrt of sum col-wise of squared values

		plt.figure(figsize=(10, 6))

		plt.errorbar(rebinned_time, rebinned_counts, yerr=rebinned_counts_err,
				marker='o', markersize=3, capsize=2, fmt=' ', color='black', label=f'BM Intensity for Run {run_string}')
	
		plt.xlabel('Time [s]')
		plt.ylabel('BM Count Rate [$16667 \mu s^{-1}$]')
		plt.title('BM Count Rate vs |B| Value')
		plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)
	
		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}BM_countRate_{run_string}_plot.pdf', format='pdf')
		plt.show()

	else:
		 if showB:
		 	print('you dont want to rebin and you want to show the |B| readback at the same time')
#		fig, (ax_counts, ax_b) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
#	
#		ax_counts.errorbar(BM_intensity[:,0], BM_intensity[:,1], yerr=np.sqrt(BM_intensity[:,1]),
#				marker='o', markersize=3, capsize=2, fmt=' ', color='black', label=f'BM Intensity for Run {run_string}')
#		ax_b.scatter(B_value[:,0], B_value[:,1], 
#				color='black', label=f'|B| value for Run {run_string}')
#	
#		plt.xlabel('Time [s]')
#		ax_counts.set_ylabel('BM Count Rate [$16667 \mu s^{-1}$]')
#		ax_counts.set_ylabel('|B| Readback Value [Gauss]')
#		plt.suptitle('BM Count Rate vs |B| Value')
#		ax_counts.legend(loc='lower right', fontsize=settings.fontsize*0.8)
#	
#		plt.tight_layout()
#		plt.savefig(f'{settings.output_dir}_BM_countRate_{run_string}_plot.pdf', format='pdf')
#		plt.show()

def plot_BM_intensity_vs_time(selected_runs):
	"""
	Function to plot the `Beam Monitor` (ID 2048) Average value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
	"""

	plt.figure(figsize=(10, 6))

	duration = selected_runs['Duration [Hours]']
	midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

	BM_avg_intensity, BM_avg_intensity_err = selected_runs['BM Average Intensity'], selected_runs['BM Average Intensity Err']

	plt.errorbar(midpoint_time, BM_avg_intensity, xerr=duration/ 2, yerr=BM_avg_intensity_err,
			fmt=' ', marker='o', markersize=5, capsize=2, color='black', label=f'BM Average Intensity for Selected Runs')

	plt.xlabel('Start Time [Hrs]')
	plt.ylabel('Count Rate [$(16667 \mu s)^{-1}$]')
	plt.title('Count Rate vs. Start Time')
	plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

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
	cmap = matplotlib.colormaps['RdYlGn']

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

			BM_avg_intensity, BM_avg_intensity_err = selected_runs['BM Average Intensity'], selected_runs['BM Average Intensity Err']
			countrates, countrates_err = selected_runs[roi_str]/ duration, selected_runs[roi_err_str]/ duration

			ratio = countrates/ BM_avg_intensity
			ratio_err = np.sqrt((countrates_err/BM_avg_intensity)**2 + (ratio)**2 * (BM_avg_intensity_err/BM_avg_intensity)**2) 

			m, b = np.polyfit(midpoint_time, ratio, 1)
			chisq = calculate_chi_square(midpoint_time, ratio, ratio_err, m, b) 

			plot_fit(ax_counts, midpoint_time, m, b, chisq)	

		for i, B_value in enumerate(settings.B_select):
			runs_for_B = selected_runs[selected_runs['Field'] == B_value]
			duration = runs_for_B['Duration [Hours]']
			midpoint_time = runs_for_B['Start Time Elapsed [Hours]'] + duration/ 2

			# Plot Normalized countrate on top figure
			BM_avg_intensity, BM_avg_intensity_err = runs_for_B['BM Average Intensity'], runs_for_B['BM Average Intensity Err']
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

		plt.xlabel('Start Time [Hrs]')
		ax_counts.set_ylabel('Count Rate/ BM Count Rate', fontsize=settings.fontsize*0.8)
		if showBM:
			ax_bm.set_ylabel('BM Average Count Rate [$hours^{-1}$]', fontsize=settings.fontsize*0.8)
			ax_bm.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0,
				fancybox=True, shadow=True, fontsize=settings.fontsize*0.8)

		plt.suptitle(f'BM Normalized Count Rate vs. Start Time, {roi_str}')
		ax_counts.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0,
			fancybox=True, shadow=True, fontsize=settings.fontsize*0.8)

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
	cmap = matplotlib.colormaps['tab10']

	for roi, roi_str, roi_err_str in roi_list:
		plt.figure(figsize=(10, 6))

		for i, B_value in enumerate(settings.B_select):

			# TODO: for each roi, calculated summary of effect for each |B| value
			print('tmp')

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

		plt.xlabel('|B| Value [T]')
#		plt.ylabel(f'{roi_str} Total Counts')
#		plt.title('|B| effect')
		plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}B_effect_plot.pdf', format='pdf')
		plt.show()

