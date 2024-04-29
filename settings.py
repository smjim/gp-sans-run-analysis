## Analysis constants 
# Data paths 
run_summary_file = '../../DATA/ONCAT_summary_HFIR_CG2_IPTS-27957.csv'
data_dir = '../../DATA/' # MNO files stored in subdirectory `MNO_files` and processed histogram files will be stored in `processed_data` 
output_dir = './output/'

# Define ROI limits 
# Pixel definition [x0 x1 y0 y1]
ROI_0 = [0, 96, 0, 256]     # Full Detector
ROI_1 = [40, 57, 100, 140]  # Rough main region
ROI_2 = [46, 50, 116, 125]  # TODO peak exclusion region(?)

# Selection criteria
min_duration = 59/60 # [Hours]
B_select = [-4.8, -3.6, -2.4, 0, 2.4, 3.6, 4.8] # |B| Field [T]
ramping_included = False 

## Plotting variables
fontsize=18
