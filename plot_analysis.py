import matplotlib.pyplot as plt
import numpy as np
import matplotlib

import settings

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

			#total_counts, total_counts_err = runs_for_B['Total Counts'], np.sqrt(runs_for_B['Total Counts'])
			counts, counts_err = runs_for_B[roi_str], runs_for_B[roi_err_str]

			color = cmap(i / len(settings.B_select))  # Sample color from colormap
			#plt.errorbar(midpoint_time, counts, xerr=half_duration, yerr=counts_err,
			#	   capsize=2, ms=3, color=color, fmt=' ', label=f'|B| = {B_value}')
			plt.scatter(midpoint_time, counts,
					s=30, color=color, label=f'|B| = {B_value}')

		plt.xlabel('Start Time [Hrs]')
		plt.ylabel(f'{roi_str} Total Counts')
		plt.title('Counts vs. Start Time')
		plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_counts_plot.pdf', format='pdf')
		plt.show()

def plot_countrates_vs_time(selected_runs, roi_list):
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

			countrates, countrates_err = runs_for_B[roi_str]/ duration, runs_for_B[roi_err_str]/ duration

			color = cmap(i / len(settings.B_select))  # Sample color from colormap
			plt.scatter(midpoint_time, countrates,
					s=30, color=color, label=f'|B| = {B_value}')

		plt.xlabel('Start Time [Hrs]')
		plt.ylabel(roi_str + ' Count Rate [$hours^{-1}$]')
		plt.title('Count Rate vs. Start Time')
		plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_countRate_plot.pdf', format='pdf')
		plt.show()

# TODO make this!
#plot_BM_intensity_run(89847)
def plot_BM_intensity_run(run_string):
	"""
	Function to plot the `Beam Monitor` (ID 2048) values as a function of time for a selected run.

	Parameters:
		run_string (string): Selected run for plotting
	"""

	#filename = 
	# TODO make this!

	plt.figure(figsize=(10, 6))

	duration = selected_runs['Duration [Hours]']
	midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

	BM_avg, BM_avg_err = selected_runs['BM Average Intensity']	

	plt.errorbar(midpoint_time, BM_avg, yerr=BM_avg_err,
			ms=30, color='black', label=f'BM Average Intensity for Selected Runs')

	plt.xlabel('Start Time [Hrs]')
	plt.ylabel(roi_str + ' Count Rate [$hours^{-1}$]')
	plt.title('Count Rate vs. Start Time')
	plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

	plt.tight_layout()
	plt.savefig(f'{settings.output_dir}_BM_avg_intensity_plot.pdf', format='pdf')
	plt.show()

def plot_BM_intensity_vs_time(selected_runs):
	"""
	Function to plot the `Beam Monitor` (ID 2048) Average value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
	"""

	plt.figure(figsize=(10, 6))

	duration = selected_runs['Duration [Hours]']
	midpoint_time = selected_runs['Start Time Elapsed [Hours]'] + duration/ 2

	BM_avg = selected_runs['BM Average Intensity']	
	BM_avg_err = np.sqrt(BM_avg) 


	plt.errorbar(midpoint_time, BM_avg, yerr=BM_avg_err,
			fmt=' ', marker='o', markersize=5, capsize=2, color='black', label=f'BM Average Intensity for Selected Runs')

	plt.xlabel('Start Time [Hrs]')
	plt.ylabel('Count Rate [$(16667 \mu s)^{-1}$]')
	plt.title('Count Rate vs. Start Time')
	plt.legend(loc='lower right', fontsize=settings.fontsize*0.8)

	plt.tight_layout()
	plt.savefig(f'{settings.output_dir}BM_avg_intensity_plot.pdf', format='pdf')
	plt.show()
