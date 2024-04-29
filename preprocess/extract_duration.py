# James Rogers
# January 22, 2024
# Extract run duration [s] from ONCAT summary file

import csv
import sys

def extract_run_duration(input_filename, input_run_number):
	# Read CSV data from the input file
	with open(input_filename, 'r') as csvfile:
		reader = csv.DictReader(csvfile)

		# Iterate through rows and find the "Duration" for the specified run number
		for row in reader:
			run_number = int(row["Run #"])
			if run_number == input_run_number:
				duration_string = row["Duration"]

				# Duration string is of the form HH:MM:SS
				duration_parts = list(map(int, duration_string.split(':')))
				duration_seconds = duration_parts[0] * 3600 + duration_parts[1] * 60 + duration_parts[2]

				# Print duration (s)
				#print(f"Run {input_run_number} Duration: {duration_seconds} seconds")
				return duration_seconds

		# If the run number is not found
		print(f"Run {input_run_number} not found in the ONCAT summary.")
		return None 

if __name__ == "__main__":
	# Example usage:
	# python script.py input_file.csv 123456
	if len(sys.argv) != 3:
		print("Usage: python script.py input_file.csv run_number")
		sys.exit(1)

	input_filename = sys.argv[1]
	input_run_number = int(sys.argv[2])

	duration = extract_run_duration(input_filename, input_run_number)

	print(duration)
