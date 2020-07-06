from   astropy.io import fits
import copy
import csv
from astropy.table import Table
import matplotlib.pyplot as plt
from   matplotlib.colors import LogNorm
import numpy as np
from pydl.pydlutils.spheregroup import spherematch
import os,sys
import statistics
import time

def run_sextractor(catalog_file, image_file1, image_file2, config_file, weight_file1, weight_file2):
    # run sextractor to generate catalog
    command = "sex " + image_file1 + "," + image_file2 + " -c " + config_file + " -WEIGHT_IMAGE " + weight_file1 + "," + weight_file2
    print("running command " + command)
    os.system(command)

    # rename the catalog
    command = "mv output.fits " + catalog_file
    os.system(command)

# precondition : config and parameter files are copied to where the python script is run.
def run_sextractor_for_all():
    # Loop on all galaxies
    for g in galaxy_names:
        # Run sextractor for original images
        base_file_name = image_dir_name + "/" + g + "/" + g + "_"
        image_file1 = base_file_name + "g_" + image_file_suffix + ".fits"           # Always use green filter image for source detection
        weight_file1 = weight_dir_name + "/" + g + "/" + g + "_" + "g_sig.fits"     # Weight file of green filter image
        for f in filter_names:
            image_file2 = base_file_name + f + "_" + image_file_suffix + ".fits"    # File for measurement
            weight_file2 = weight_dir_name + "/" + g + "/" + g + "_" + f + "_sig.fits"        # Corresponding weight file
            catalog_file =  g + "_" + f + "_" + image_file_suffix + "_cat.fits"
            run_sextractor(catalog_file, image_file1, image_file2, config_file_name, weight_file1, weight_file2)

            tmp_dir_name = result_dir_name + "/" + g + "/"
            if not os.path.exists(tmp_dir_name):
                os.mkdir(tmp_dir_name)
            command = "mv " + catalog_file + " " + tmp_dir_name
            os.system(command)

        # Run sextractor for images with added stars
        base_file_name = result_dir_name + "/" + g + "/" + g + "_"
        for mag in magnitudes:
            for suffix in range(0, num_suffix):
                image_file1 = base_file_name + "g_addstar" + str(mag) + "_" + str(suffix) + ".fits"  # Always use green filter image for source detection
                weight_file1 = weight_dir_name + "/" + g + "/" + g + "_" + "g_sig.fits"                 # Weight file of green filter image
                for f in filter_names:
                    image_file2 = base_file_name + f + "_addstar" + str(mag) + "_" + str(suffix) + ".fits"
                    weight_file2 = weight_dir_name + "/" + g + "/" + g + "_" + f + "_sig.fits"

                    catalog_file =  g + "_" + f + "_addstar" + str(mag) + "_" + str(suffix)  + "_cat.fits"
                    run_sextractor(catalog_file, image_file1, image_file2, config_file_name, weight_file1, weight_file2)

                    tmp_dir_name = result_dir_name + "/" + g + "/"
                    if not os.path.exists(tmp_dir_name):
                        os.mkdir(tmp_dir_name)
                    command = "mv " + catalog_file + " " + tmp_dir_name
                    os.system(command)

def plot_completeness_function(x_data, y_data):
    color = "black"
    plt.figure(1)
    plt.scatter(x=x_data, y=y_data, c='r', marker='x', label='GC detection ratio', alpha=0.5)
    plt.ylim(0.0, 1.0)
    plt.legend()
    plt.title('Completeness function', fontsize=14)
    plt.xlabel("Magnitudes")
    plt.ylabel("Detection ratio")
    plt.show()

def completeness_test():
    for g in galaxy_names:
        filter_index = 0
        for f in filter_names:
            num_star_per_mag_bin = [0 for x in range(len(magnitudes))]

            # Get number of detected sources from sextractor produced catalog of original image
            orig_fits_file = result_dir_name + "/" + g + "/" + g + "_" + f + "_modsub2_cat.fits"
            orig_cat = fits.open(orig_fits_file)
            orig_cat.info()
            orig_cat_data = Table(orig_cat[2].data)
            orig_num_detected = len(orig_cat_data)

            # Get number of detected sources from sextractor produced catalog of image with added stars
            mag_index = 0
            for mag in magnitudes:
                for suffix in range(0, num_suffix):
                    tmp_dir_name = result_dir_name + "/" + g + "/"
                    catalog_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(suffix) + "_cat.fits"
                    addstar_cat = fits.open(catalog_file)
                    addstar_cat.info()
                    addstar_cat_data = Table(addstar_cat[2].data)
                    num_star_per_mag_bin[mag_index] += len(addstar_cat_data) - orig_num_detected

                # Average numbers for all catalogs per magnitude bin, divided by true number of added stars
                num_star_per_mag_bin[mag_index] /= (num_suffix * (NUM_GRID_POINTS*NUM_GRID_POINTS))
                mag_index += 1

            filter_index += 1
            plot_completeness_function(magnitudes, num_star_per_mag_bin)

# main
filter_names = ("g")
image_dir_name = "../images"
result_dir_name = "../results"
config_file_name = "ngvs.sex"
galaxy_names = ["VCC0940"]
magnitudes = (23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28)
NUM_STARS_TO_ADD_PER_MAG = 100
INDEX_4 = 5
INDEX_8 = 9
num_suffix = 4
NUM_GRID_POINTS = 30    #number of grid points along x or y

filter_index = 0

weight_dir_name = "../images"
image_file_suffix = "modsub2"

if len(sys.argv) == 3:
    image_dir_name = sys.argv[1]
    image_file_suffix = sys.argv[2]

#run_sextractor_for_all()
completeness_test()