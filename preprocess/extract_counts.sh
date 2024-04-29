#!/bin/bash

# Define the range of values for run_rum 
first=89821
last=89976

# Define Data directory
data_dir="../../../DATA/"

# Run all runs in range
for ((run_num = $first; run_num <= $last; run_num++));
do
	# Find run duration from ONCAT summary
	duration=$(python3 extract_duration.py ONCAT_summary_HFIR_CG2_IPTS-27957.csv ${run_num})

    # Construct the command with the current value of run_num 
	command="./reduced_xbin_MNO_to_intensityPSD DATA/MNO_GPSANS_${run_num}.txt ${data_dir}DetIDmap.csv ${data_dir}processed_data/reduced_xbin_${run_num}_output.dat ${duration}"

    # Run the command
    echo "Running command ${command}"
    $command

done
