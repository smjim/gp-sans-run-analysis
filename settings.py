## Analysis constants 

# Data paths 
run_summary_file = '../../DATA/ONCAT_summary_HFIR_CG2_IPTS-27957.csv'
output_summary_file = './output/output_summary.csv'

data_dir = '../../DATA/'						# Contains ONCAT Summary, DetIDmap.csv
mno_dir = '../../DATA/MNO_files'				# Contains MNO data files
processed_outdir = '../../DATA/processed_data'	# Contains McStas-format output histograms
output_dir = './output/'

# Define ROI limits: Pixel definition [x0 x1 y0 y1] 
ROI_0 = [0, 96, 0, 256]     # Full Detector
ROI_1 = [40, 57, 100, 140]  # Rough main region
ROI_2 = [46, 50, 116, 125]  # TODO peak exclusion region(?)

# Selection criteria
min_duration = 59/60 # [Hours]
B_select = [0] #[-4.8, -3.6, -2.4, 0, 2.4, 3.6, 4.8] # |B| Field [T]
ramping_included = False 
run_select = None #[89824, 89949, 89950]
