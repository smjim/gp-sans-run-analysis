## Analysis constants 
# Data paths 
run_summary_file = '../../DATA/ONCAT_summary_HFIR_CG2_IPTS-27957.csv'
data_dir = '../../DATA/processed_data/' # Preprocessed data using scripts from https://github.com/smjim/mcStas_optimization_reduction
output_dir = '../figures-analysis/'

# Define ROI limits 
# Pixel definition [x0 x1 y0 y1]
ROI_0 = [0, 96, 0, 256]     # Full Detector
ROI_1 = [40, 57, 100, 140]  # Rough main region
ROI_2 = [46, 50, 116, 125]  # TODO peak exclusion region(?)

# Selection criteria
min_duration = 1/60 # [Hours]
B_select = [-4.8, -3.6, -2.4, 0, 2.4, 3.6, 4.8] # |B| Field [T]

## Plotting variables
fontsize=18
