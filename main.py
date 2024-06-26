# James Rogers
# 4-28-2024
# Analyze runs for given |B| all at the same time:
#	output counts per second/ ROI for each ROI
#	output summary plots for each |B| value

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import pandas as pd

import settings

from analyze_MNO_files import * 
from plot_analysis import *

## Define functions
def load_run_summary(file_path):
	run_summary = pd.read_csv(file_path)
	return run_summary

def parse_start_time(run_summary):
	"""
	Function to parse the start time information from the run summary DataFrame.
	
	Parameters:
		run_summary (DataFrame): DataFrame containing the run summary data.
		
	Returns:
		Series: Series containing the parsed start time information.
	"""

	# Parse run summary
	start_time_str = run_summary['Start Time']
	duration_str = run_summary['Duration']
	
	# Extract Start Time Information
	t0_time_str = "2024/01/18 12:27:25"
	t0_time_zone = "EST"
	t0_time = pd.Timestamp(t0_time_str, tz=t0_time_zone)
	start_time = pd.to_datetime(start_time_str, format="%Y/%m/%d %H:%M:%S %Z")
	
	start_time_elapsed = start_time - t0_time
	start_time_elapsed = start_time_elapsed.dt.total_seconds() / 3600

	# Extract Duration Information
	duration_timedelta = pd.to_timedelta(duration_str)
	duration_h = duration_timedelta.dt.total_seconds() / 3600

	return start_time_elapsed, duration_h

def select_runs_by_B(run_summary, B_values):
	"""
	Function to select runs from the run summary DataFrame based on given |B| values.
	
	Parameters:
		run_summary (DataFrame): DataFrame containing the run summary data.
		B_values (list): List of |B| values to select runs for.
		
	Returns:
		DataFrame: DataFrame containing the selected runs.
	"""
	return selected_runs

## Main function
def main():
	# Load run summary
	run_summary = load_run_summary(settings.run_summary_file)
	
	# Parse start time
	run_summary['Start Time Elapsed [Hours]'], run_summary['Duration [Hours]'] = parse_start_time(run_summary)
	
	# -----------------------
	# Filter and select runs 
	# -----------------------

	selected_runs = run_summary.copy()
	selected_runs = selected_runs[selected_runs['Field'].isin(settings.B_select)]				# Select runs by |B| values in B_select list
	selected_runs = selected_runs[selected_runs['Duration [Hours]'] > settings.min_duration]	# Filter runs by minimum duration
	selected_runs = selected_runs[selected_runs['Ramping'] == settings.ramping_included] 		# Choose whether to include ramping runs

	# -----------------------
	# Calculate quantities and store to DataFrame
	# -----------------------

	# Calculate counts/ ROI for each ROI definition
	ROI_list = [[settings.ROI_0, 'ROI 0', 'ROI 0 Err'],
				[settings.ROI_1, 'ROI 1', 'ROI 1 Err'],
				[settings.ROI_2, 'ROI 2', 'ROI 2 Err']]

	# Calculate ROI counts and BM intensity from MNO at the same time
	selected_runs = analyze_MNO(selected_runs, ROI_list, settings.data_dir, settings.output_summary_file)

	# Calculate only ROI counts or BM intensity
	#selected_runs = find_ROI_counts_MNO(selected_runs, ROI_list, settings.data_dir)
	#selected_runs = find_BM_intensity_MNO(selected_runs)
	
	# -----------------------
	# Plot outputs
	# -----------------------

	# Plot counts vs. time for selected runs, with |B| field colored
	#plot_counts_vs_time(selected_runs, ROI_list)
	#plot_countrates_vs_time(selected_runs, ROI_list, fit=True)

	# Plot BM Intensity-corrected countrates
	plot_countrates_vs_time_normalized(selected_runs, ROI_list, fit=True, showBM=False)

	# Plot BM Avg Intensity Variation across selected runs 
	#plot_BM_intensity_vs_time(selected_runs)

	# Plot BM Intensity vs Time for a single run
	#plot_run_stats(89847, rebin=1000, showB=False)

	# Plot summary of ROI counts for each |B| value over entire experiment
	# TODO: complete this function
	#plot_effect_analysis(selected_runs, ROI_list)	

	# TODO: plot counts vs time for each ROI for given |B| value (x: time, y: counts, color: ROI def)
	#plot_countrates_vs_time_B(selected_runs)
	# (not super useful tbh)

if __name__ == "__main__":
	# Plot settings
	matplotlib.rcParams.update({'font.size': settings.fontsize})
	main()

