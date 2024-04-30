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
			plt.errorbar(midpoint_time, countrates, xerr=duration/ 2, yerr=countrates_err,
					marker='o', markersize=3, capsize=2, fmt=' ', color=color, label=f'|B| = {B_value}')

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

	BM_avg_intensity, BM_avg_intensity_err = runs_for_B['BM Average Intensity'], runs_for_B['BM Average Intensity Err']

	plt.errorbar(midpoint_time, BM_avg_intensity, xerr=duration/ 2, yerr=BM_avg_intensity_err,
			marker='o', markersize=3, capsize=2, fmt=' ', color='black', label=f'BM Average Intensity for Selected Runs')

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
def plot_countrates_vs_time_normalized(selected_runs, roi_list):
	"""
	Function to plot the 'Total Counts' value as a function of 'Start Time' for selected runs.

	Parameters:
		selected_runs (DataFrame): DataFrame containing the selected runs.
		roi_list (array): Array containing ROI definitions, ROI counts column header, ROI counts err column header
	"""

	# Establish colorscheme
	cmap = matplotlib.colormaps['RdYlGn']

	for roi, roi_str, roi_err_str in roi_list:

		fig, (ax_counts, ax_bm) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

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
			ax_bm.errorbar(midpoint_time, BM_avg_intensity, xerr=duration/ 2, yerr=BM_avg_intensity_err,
					fmt=' ', marker='o', markersize=5, capsize=2, color='black')

		plt.xlabel('Start Time [Hrs]')
		ax_counts.set_ylabel(roi_str + ' Count Rate/ BM Count Rate [$hours^{-1}$]')
		ax_bm.set_ylabel('BM Average Count Rate [$hours^{-1}$]')
		plt.suptitle('BM Normalized Count Rate vs. Start Time')

		ax_counts.legend(loc='lower right', fontsize=settings.fontsize*0.8)

		plt.tight_layout()
		plt.savefig(f'{settings.output_dir}{roi_str}_normalized_countRate_plot.pdf', format='pdf')
		plt.show()
