To run the SMiRF protocol, do the following:

1. Generate a list of iRF interactions using irf_on_data.R and .csv file (sample shown in rawdata_noclay.csv)

2. Copy and paste interactions from output of irf_on_data.R into appropriate HAL*.R files (one for 2nd, 3rd, and 4th order interactions; higher order interactions unlikely to run at all with hal9001)

3. Combine the output of HAL*.R files with surface_combiner.py into one single file

4. Run ROSM2.py on combined_HAL_points.csv (output of surface_combiner.py); this will generate 4 ROSMs (Ordinary Least Squares, Lasso, Ridge, and Forward-Backward)

Codes may need to be modified for non-spatial data. First two columns are assumed to be positional and not features.
