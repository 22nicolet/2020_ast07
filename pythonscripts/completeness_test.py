from   astropy.io import fits
import copy
import csv
import config as cfg
from astropy.table import Table
import matplotlib.pyplot as plt
from   matplotlib.colors import LogNorm
import numpy as np
from pydl.pydlutils.spheregroup import spherematch
import os,sys
import statistics
import time
import astropy.wcs as wcs
from scipy.optimize import curve_fit


def same_source(x1, y1, x2, y2):
    if abs(x1-x2)<=RADIUS and abs(y1-y2)<=RADIUS:
        return True
    else:
        return False

# Return number points described by (x_data, y_data) within radius of (x_center, y_center)
def num_src_in_radius(x_center, y_center, radius, x_data, y_data):
    radius_square = radius*radius
    count = 0
    for i in range(0, len(x_data)):
        dist_square = (x_data[i] - x_center)*(x_data[i] - x_center) + (y_data[i] - y_center)*(y_data[i] - y_center)
        if dist_square <= radius_square:
            count += 1
    return [count, len(x_data)-count]

# generates command to run SExtractor
def run_sextractor(catalog_file, image_file1, image_file2, config_file, weight_file1, weight_file2):
    # run sextractor to generate catalog
    command = "sex " + image_file1 + "," + image_file2 + " -c " + config_file + " -WEIGHT_IMAGE " + weight_file1 + "," + weight_file2
    print("running command " + command)
    os.system(command)

    # rename the catalog
    command = "mv output.fits " + catalog_file
    os.system(command)

# precondition : config and parameter files are copied to where the python script is run.
# runs SExtractor for all original galaxies and images with artificial stars to generate catalog files
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
        for n in range(0, cfg.NUMBER_NOISES):
            for mag in magnitudes:
                for suffix in range(0, num_suffix):
                    image_file1 = base_file_name + "g_addstar" + str(mag) + "_" + str(suffix) + "_" + str(n) + ".fits"  # Always use green filter image for source detection
                    weight_file1 = weight_dir_name + "/" + g + "/" + g + "_" + "g_sig.fits"                 # Weight file of green filter image
                    for f in filter_names:
                        image_file2 = base_file_name + f + "_addstar" + str(mag) + "_" + str(suffix) + "_" + str(n) + ".fits"
                        weight_file2 = weight_dir_name + "/" + g + "/" + g + "_" + f + "_sig.fits"

                        catalog_file =  g + "_" + f + "_addstar" + str(mag) + "_" + str(suffix)  + "_" + str(n) + "_cat.fits"
                        run_sextractor(catalog_file, image_file1, image_file2, config_file_name, weight_file1, weight_file2)

                        tmp_dir_name = result_dir_name + "/" + g + "/"
                        if not os.path.exists(tmp_dir_name):
                            os.mkdir(tmp_dir_name)
                        command = "mv " + catalog_file + " " + tmp_dir_name
                        os.system(command)

# plots the completeness function
def plot_completeness_function(x_data, y_data, popt, color, mark):
    color = "black"
    #plt.figure(1)
    plt.scatter(x=x_data, y=y_data, c=color, marker=mark, label='GC detection ratio', alpha=0.5)
    plt.ylim(0.0, 1.0)
    plt.xticks(np.arange(22.5, 28.5, 0.5))
    plt.yticks(np.arange(0.0, 1.0, 0.1))
    plt.grid(color = 'grey', linestyle = '-.', linewidth = 1)
    plt.title('Completeness function (NUM_GRID_POINTS = ' + str(cfg.NUM_GRID_POINTS) + ')', fontsize=14)
    plt.xlabel("Magnitudes")
    plt.ylabel("Detection ratio")

    x = np.arange(x_data[0], x_data[-1], 0.1)   #create a new array with magnitudes, each spaced by 0.1
    y = pritchet_formula(x, *popt)
    plt.plot(x, y, 'g--', label='fit: Alpha_f%5.3f, T1_limit = %5.3f, Factor_max = %5.3f' % tuple(popt))

    plt.legend()
    plt.show()

def plot_completeness_function_all(x_data_inside, y_data_inside, x_data_outside, y_data_outside, x_data, y_data, popt, color, mark):
    color = "black"
    #plt.figure(1)
    plt.scatter(x=x_data, y=y_data, c=color, marker=mark, label='GC detection ratio', alpha=0.5)
    plt.ylim(0.0, 1.0)
    plt.xticks(np.arange(22.5, 28.5, 0.5))
    plt.yticks(np.arange(0.0, 1.0, 0.1))
    plt.grid(color = 'grey', linestyle = '-.', linewidth = 1)
    plt.title('Completeness function (NUM_GRID_POINTS = ' + str(cfg.NUM_GRID_POINTS) + ')', fontsize=14)
    plt.xlabel("Magnitudes")
    plt.ylabel("Detection ratio")

    x = np.arange(x_data[0], x_data[-1], 0.1)   #create a new array with magnitudes, each spaced by 0.1
    y = pritchet_formula(x, *popt)
    plt.plot(x, y, 'g--', label='fit: Alpha_f%5.3f, T1_limit = %5.3f, Factor_max = %5.3f' % tuple(popt))

    plt.scatter(x=x_data_inside, y=y_data_inside, c='r', marker='+', label='GC detection ratio', alpha=0.5)
    plt.scatter(x=x_data_outside, y=y_data_outside, c='b', marker='*', label='GC detection ratio', alpha=0.5)

    plt.legend()
    plt.show()

def pritchet_formula(t1, alpha_f, t1_limit, factor_max):
    return 0.5 * (1 - (alpha_f*(t1 - t1_limit)/np.sqrt(1 + alpha_f * alpha_f * (t1 - t1_limit) * (t1 - t1_limit)))) * factor_max

# finds detection ratio:
# number of artificial stars detected = number of stars detected in the images with added psf - number detected in the original image
# average the difference among the four steps in each magnitude
# divide difference by total number of artificial stars added to get detection ratio
def completeness_test():
    for g in galaxy_names:
        filter_index = 0
        for f in filter_names:
            start = time.time()
            num_star_per_mag_bin = [0 for x in range(len(magnitudes))]
            num_star_per_mag_bin_inside = [0 for x in range(len(magnitudes))]
            num_star_per_mag_bin_outside = [0 for x in range(len(magnitudes))]

            # Get number of detected sources from sextractor produced catalog of original image
            orig_fits_file = result_dir_name + "/" + g + "/" + g + "_" + f + "_modsub2_cat.fits"
            orig_cat = fits.open(orig_fits_file)
            #orig_cat.info()
            orig_cat_data = Table(orig_cat[2].data)
            orig_num_detected = len(orig_cat_data)

            for n in range(0, cfg.NUMBER_NOISES):
                # Get number of detected sources from sextractor produced catalog of image with added stars
                mag_index = 0

                for mag in magnitudes:
                    truth_inside_total = 0
                    truth_outside_total = 0
                    for suffix in range(0, num_suffix):
                        tmp_dir_name = result_dir_name + "/" + g + "/"
                        catalog_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(suffix) + "_" + str(n) + "_cat.fits"
                        image_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(suffix) + "_" + str(n) + ".fits"
                        w = wcs.WCS(image_file)
                        x_center, y_center = w.all_world2pix(cfg.RA_GALAXY_CENTER[0], cfg.DEC_GALAXY_CENTER[0], 1)
                        addstar_cat = fits.open(catalog_file)
                        #addstar_cat.info()
                        addstar_cat_data = Table(addstar_cat[2].data)
                        x_cat = addstar_cat_data['X_IMAGE']
                        y_cat = addstar_cat_data['Y_IMAGE']
                        ra_cat_col = addstar_cat_data['ALPHA_J2000']
                        dec_cat_col = addstar_cat_data['DELTA_J2000']
                        ra_cat = np.array(ra_cat_col)
                        dec_cat = np.array(dec_cat_col)
                        addstar_cat.close()

                        # read truth csv file
                        x_truth = []
                        y_truth = []
                        truth_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(suffix) + "_" + str(n) + "_truth.csv"
                        with open(truth_file) as csv_file:
                            csv_reader = csv.reader(csv_file, delimiter=',')
                            for row in csv_reader:
                                x_truth.append(float(row[0]))
                                y_truth.append(float(row[1]))

                        [truth_inside, truth_outside] = num_src_in_radius(x_center, y_center, cfg.CR_RADIUS, x_truth, y_truth)
                        truth_inside_total += truth_inside
                        truth_outside_total += truth_outside
                        #translates x y coordinates into ra dec (wcs world coordinates)
                        ra_truth, dec_truth = w.all_pix2world(x_truth, y_truth, 1)

                        # spherematch on each inserted artificial GC (truth) to check if it's detected.
                        # If yes, add to a set for quick look up.
                        coords = ()
                        coords = spherematch(ra_cat, dec_cat, ra_truth, dec_truth, matchlength=0.0001, chunksize=0.0005)
                        num_detected = len(coords[0])

                        x_detected = []
                        y_detected = []

                        for index in coords[1]:
                            x_detected.append(x_truth[index])
                            y_detected.append(y_truth[index])
                        [detected_inside, detected_outside] = num_src_in_radius(x_center, y_center, cfg.CR_RADIUS, x_detected, y_detected)

                        num_star_per_mag_bin_inside[mag_index] += detected_inside
                        num_star_per_mag_bin_outside[mag_index] += detected_outside
                        num_star_per_mag_bin[mag_index] += num_detected
                        #print(str(mag) + ", " + str(suffix) + ", " + str(num_detected))
                    # Average numbers for all catalogs per magnitude bin, divided by true number of added stars
                    num_star_per_mag_bin[mag_index] /= (1.0 * num_suffix * (cfg.NUM_GRID_POINTS*cfg.NUM_GRID_POINTS))
                    num_star_per_mag_bin_outside[mag_index] /= 1.0 * truth_outside_total
                    num_star_per_mag_bin_inside[mag_index] /= 1.0 * truth_inside_total
                    mag_index += 1

            end = time.time()
            print('Elapsed time of completeness test = ' + str(end - start) + ' seconds')
            popt, pcov = curve_fit(pritchet_formula, magnitudes, num_star_per_mag_bin)
            filter_index += 1

            completeness_csv_file = result_dir_name + "/" + g + "/" + g + "_" + f + "_completeness.csv"
            with open(completeness_csv_file, 'w') as csvout:
                csv_writer = csv.writer(csvout)
                csv_writer.writerow(magnitudes)
                csv_writer.writerow(num_star_per_mag_bin)
                csv_writer.writerow(['alpha_f', popt[0]])
                csv_writer.writerow(['t1_limit', popt[1]])
            plt.figure(1)
            #plot_completeness_function(magnitudes, num_star_per_mag_bin, popt, 'r', 'x')
            #plot_completeness_function(magnitudes, num_star_per_mag_bin_inside, popt, 'g', '+')
            #plot_completeness_function(magnitudes, num_star_per_mag_bin_outside, popt, 'b', '*')
            plot_completeness_function_all(magnitudes, num_star_per_mag_bin_inside, \
                                       magnitudes, num_star_per_mag_bin_outside, \
                                       magnitudes, num_star_per_mag_bin, \
                                       popt, 'b', '*')


# main
filter_names = ("g")
image_dir_name = "../../../SIP2020/images"
result_dir_name = "../../../SIP2020/results"
config_file_name = "ngvs.sex"
galaxy_names = ["VCC0940"]
magnitudes = (23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28)
NUM_STARS_TO_ADD_PER_MAG = 100
INDEX_4 = 5
INDEX_8 = 9
num_suffix = 4

RADIUS = 1  # If x1 and x2 are less than RADIUS pixels away, consider them to be the same source

filter_index = 0

weight_dir_name = "../../../SIP2020/images"
image_file_suffix = "modsub2"

if len(sys.argv) == 3:
    image_dir_name = sys.argv[1]
    image_file_suffix = sys.argv[2]

start = time.time()
#run_sextractor_for_all()
end = time.time()
print('Elapsed time = '+str(end-start)+' seconds')
completeness_test()
