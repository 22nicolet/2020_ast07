from astropy.io import fits
from astropy.table import Table
import os
import sys

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
        base_dir_name = image_dir_name + "/" + g
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

# main
filter_names = ("u", "g", "r", "i", "z")
image_dir_name = "../images"
weight_dir_name = "../images"
image_file_suffix = "modsub2"
result_dir_name = "../results"
config_file_name = "ngvs.sex"
galaxy_names = ["VCC0940"]

if len(sys.argv) == 3:
    image_dir_name = sys.argv[1]
    image_file_suffix = sys.argv[2]

run_sextractor_for_all()
