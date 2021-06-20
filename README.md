# 2021_ast09
Set up directories as seen in DIRECTORY_TREE.png

images  - Original images. Each galaxy is in a separate sub-directory such as VCC0940
results - result files from sextractor, or images with added stars. Each galaxy is in a separate sub-directory such as VCC0940
pythonscripts - python scripts, as well as the run directory

#commands to run for the processing pipeline
git clone https://github.com/22nicolet/2020_ast07.git
cd 2020_ast07/
mkdir results
cd pythonscripts
python add_stars.py          # Add artifical stars to model subtracted galaxy images 
python completeness_test.py  # Run completeness tests for the images produced by add_stars.py
python CCD.py                # Identify GC candidates, compute surface density, correction factor, number of GCs per annulus, and specific frequency
python sn_plot.py ../results/gcs_info.csv 2008_ACSVCS_galaxy_table.csv 2007_HST_WFPC2.csv  # Plot NGVS result against results from two previous studies. 
