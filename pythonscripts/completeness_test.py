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
import time
import astropy.wcs as wcs
from scipy.optimize import curve_fit

def same_source(x1, y1, x2, y2):
    if abs(x1-x2)<=RADIUS and abs(y1-y2)<=RADIUS:
        return True
    else:
        return False

# Return number points described by (x_data, y_data) within radius of (x_center, y_center)
def num_src_in_ring(x_center, y_center, radius_inner, radius_outer, x_data, y_data):
    radius_inner_square = radius_inner * radius_inner
    radius_outer_square = radius_outer * radius_outer
    count = 0
    for i in range(0, len(x_data)):
        dist_square = (x_data[i] - x_center)*(x_data[i] - x_center) + (y_data[i] - y_center)*(y_data[i] - y_center)
        if dist_square <= radius_outer_square and dist_square >= radius_inner_square:
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
    start = time.time()
    for g in galaxy_names:
        # Run SExtractor for original images
        base_file_name = image_dir_name + "/" + g + "/" + g + "_"
        image_file1 = base_file_name + "g_" + image_file_suffix + "0.fits"           # Always use green filter image for source detection
        suffix = 9
        while not os.path.exists(image_file1):
            image_file1 = base_file_name + "g_" + image_file_suffix + str(suffix) + ".fits"
            suffix -= 1
            if suffix < 1:
                break
        weight_file1 = weight_dir_name + "/" + g + "/" + g + "_" + "g_sig.fits"     # Weight file of green filter image
        for f in filter_names:
            image_file2 = base_file_name + f + "_" + image_file_suffix + "0.fits"   # File for measurement
            suffix = 9
            while not os.path.exists(image_file2):
                image_file2 = base_file_name + f + "_" + image_file_suffix + str(suffix) + ".fits"
                suffix -= 1
                if suffix < 1:
                    break
            weight_file2 = weight_dir_name + "/" + g + "/" + g + "_" + f + "_sig.fits"  # Corresponding weight file
            catalog_file =  g + "_" + f + "_" + image_file_suffix + "_cat.fits"         # Drop suffix after modsub to make it easy for completeness_test
            run_sextractor(catalog_file, image_file1, image_file2, config_file_name, weight_file1, weight_file2)

            tmp_dir_name = result_dir_name + "/" + g + "/"
            if not os.path.exists(tmp_dir_name):
                os.mkdir(tmp_dir_name)
            command = "mv " + catalog_file + " " + tmp_dir_name
            os.system(command)

        # Run SExtractor for images with added stars
        base_file_name = result_dir_name + "/" + g + "/" + g + "_"
        for n in range(0, cfg.NUMBER_NOISES):
            for mag in cfg.MAGNITUDES:
                for step in range(0, num_steps):
                    image_file1 = base_file_name + "g_addstar" + str(mag) + "_" + str(step) + "_" + str(n) + ".fits"  # Always use green filter image for source detection
                    weight_file1 = weight_dir_name + "/" + g + "/" + g + "_" + "g_sig.fits"                 # Weight file of green filter image
                    for f in filter_names:
                        image_file2 = base_file_name + f + "_addstar" + str(mag) + "_" + str(step) + "_" + str(n) + ".fits"
                        weight_file2 = weight_dir_name + "/" + g + "/" + g + "_" + f + "_sig.fits"

                        catalog_file =  g + "_" + f + "_addstar" + str(mag) + "_" + str(step)  + "_" + str(n) + "_cat.fits"
                        run_sextractor(catalog_file, image_file1, image_file2, config_file_name, weight_file1, weight_file2)

                        tmp_dir_name = result_dir_name + "/" + g + "/"
                        if not os.path.exists(tmp_dir_name):
                            os.mkdir(tmp_dir_name)
                        command = "mv " + catalog_file + " " + tmp_dir_name
                        os.system(command)
    end = time.time()
    print('Elapsed time of SExtractor runs = '+str(end-start)+' seconds')

# plots the completeness function
def plot_completeness_function_all(galaxy, x_data_inside, y_data_inside, x_data_outside, y_data_outside, \
                                   poptinside, poptoutside, color, mark, annulus):
    color = "black"
    plt.ylim(0.0, 1.0)
    plt.xticks(np.arange(22.5, 28.5, 0.5))
    plt.yticks(np.arange(0.0, 1.0, 0.1))
#    plt.grid(color = 'grey', linestyle = '-.', linewidth = 1)
    plt.axhline(y=0.5, c='grey', ls=":", alpha=0.5)

    plt.title('Galaxy ' + galaxy + ' completeness function\n' \
              '(annulus=' + annulus +' Re)', fontsize=14)
    plt.xlabel("Magnitudes")
    plt.ylabel("Detection ratio")

    plt.scatter(x=x_data_inside, y=y_data_inside, c='r', marker='+', label='GC detection ratio in annulus',\
                alpha=0.5)
#    plt.scatter(x=x_data_outside, y=y_data_outside, c='b', marker='*', label='GC detection ratio out of region',\
#                alpha=0.5)

    x = np.arange(x_data_inside[0], x_data_inside[-1], 0.1)   #create a new array with magnitudes, each spaced by 0.1
    yinside = pritchet_formula(x, *poptinside)
#    youtside = pritchet_formula(x, *poptoutside)
    plt.plot(x, yinside, 'r--', label='Pritchet formula fit:\nAlpha_f=%5.3f\nT1_limit=%5.3f\nFactor_max=%5.3f' % tuple(poptinside))
#    plt.plot(x, youtside, 'b:', label='Alpha_f=%5.3f\nT1_limit=%5.3f\nFactor_max=%5.3f' % tuple(poptoutside))

    plt.legend()

def pritchet_formula(t1, alpha_f, t1_limit, factor_max):
    return 0.5 * (1 - (alpha_f*(t1 - t1_limit)/np.sqrt(1 + alpha_f * alpha_f * (t1 - t1_limit) * (t1 - t1_limit)))) * factor_max

# finds detection ratio:
# number of artificial stars detected = number of stars detected in the images with added psf - number detected in the original image
# average the difference among the four steps in each magnitude
# divide difference by total number of artificial stars added to get detection ratio
def completeness_test():
    start = time.time()
    galaxy_index = 0
    for g in galaxy_names:
        print("Completeness testing on Galaxy " + g)
        list = [8.0] + cfg.CR_RADIUS
        radius_in_pixels = np.array(list) * cfg.GALAXY_EFFECTIVE_RADII[galaxy_index] / cfg.ARCSEC_PER_PIXEL

        #Get the size of the image
        image_file = result_dir_name + "/" + g + "/" + g + "_g_addstar24.0_0_0.fits"
        image = fits.open(image_file)
        image_data = image[0].data
        image_size = len(image_data[1])

        x_center = -1
        y_center = -1
        filter_index = 0
        for f in filter_names:

            # Read SExtractor catalog, match it with truth to get list of detected artificial GCs.
            x_detected = [] # x_detected[noise_index][mag_index][step] is a list of x coordinates
            y_detected = [] # y_detected[noise_index][mag_index][step] is a list of y coordinates
            for noise_index in range(0, cfg.NUMBER_NOISES):
                # Get number of detected sources from sextractor produced catalog of image with added stars
                x_detected_noise = []
                y_decteded_noise = []
                mag_index = 0
                for mag in cfg.MAGNITUDES:
                    x_detected_mag = []
                    y_decteded_mag = []
                    for step in range(0, num_steps):
                        # read SExtractor catalog into ra_cat, dec_cat
                        x_detected_step = []
                        y_detected_step = []
                        tmp_dir_name = result_dir_name + "/" + g + "/"
                        catalog_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(step) + "_" + str(noise_index) + "_cat.fits"
                        image_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(step) + "_" + str(noise_index) + ".fits"
                        w = wcs.WCS(image_file)
                        x_center, y_center = w.all_world2pix(cfg.GALAXY_CENTER_RA[galaxy_index], cfg.GALAXY_CENTER_DEC[galaxy_index], 1)
                        radius_in_pixels[0] = max(max(x_center, image_size - 1 - x_center), max(y_center, image_size - 1 - y_center)) - 100
                        addstar_cat = fits.open(catalog_file)
                        addstar_cat_data = Table(addstar_cat[2].data)
                        x_cat = addstar_cat_data['X_IMAGE']
                        y_cat = addstar_cat_data['Y_IMAGE']
                        ra_cat_col = addstar_cat_data['ALPHA_J2000']
                        dec_cat_col = addstar_cat_data['DELTA_J2000']
                        ra_cat = np.array(ra_cat_col)
                        dec_cat = np.array(dec_cat_col)
                        addstar_cat.close()

                        # read truth csv file into ra_truth, dec_truth
                        x_truth = []
                        y_truth = []
                        truth_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(step) + "_" + str(noise_index) + "_truth.csv"
                        with open(truth_file) as csv_file:
                            csv_reader = csv.reader(csv_file, delimiter=',')
                            for row in csv_reader:
                                x_truth.append(float(row[0]))
                                y_truth.append(float(row[1]))
                        #translates x y coordinates in pixel into ra dec (wcs world coordinates)
                        ra_truth, dec_truth = w.all_pix2world(x_truth, y_truth, 1)

                        # call spherematch to get a list of detected artificial GCs
                        # coords[0] is indices matched in ra_cat/dec_cat.
                        # coords[1] is indices matched in ra_truth/dec_truth
                        coords = ()
                        coords = spherematch(ra_cat, dec_cat, ra_truth, dec_truth, matchlength=0.0001, chunksize=0.0005)

                        for index in coords[1]:
                            x_detected_step.append(x_truth[index])
                            y_detected_step.append(y_truth[index])
                        x_detected_mag.append(x_detected_step)
                        y_decteded_mag.append(y_detected_step)
                        # for step
                    x_detected_noise.append(x_detected_mag)
                    y_decteded_noise.append(y_decteded_mag)
                    mag_index += 1
                    # for mag
                x_detected.append(x_detected_noise)
                y_detected.append(y_decteded_noise)
                # for noise_index

            # Calculate detection ratio for each annulus.  For each annulus, accumulate number of
            # truth and detected artificial GCs for each magnitude bin across all noise_index and num_steps
            for i in range (0, len(radius_in_pixels) - 1):
                detected_star_per_mag_bin_inside = np.zeros((len(cfg.MAGNITUDES)), dtype=float)
                detected_star_per_mag_bin_outside = np.zeros((len(cfg.MAGNITUDES)), dtype=float)
                added_star_per_mag_bin_inside = np.zeros((len(cfg.MAGNITUDES)), dtype=float)
                added_star_per_mag_bin_outside = np.zeros((len(cfg.MAGNITUDES)), dtype=float)

                for noise_index in range(0, cfg.NUMBER_NOISES):
                    mag_index = 0
                    for mag in cfg.MAGNITUDES:
                        for step in range(0, num_steps):
                            tmp_dir_name = result_dir_name + "/" + g + "/"
                            x_truth = []
                            y_truth = []
                            truth_file = tmp_dir_name + g + "_" + f + "_addstar" + str(mag) + "_" + str(step) + "_" + str(noise_index) + "_truth.csv"
                            with open(truth_file) as csv_file:
                                csv_reader = csv.reader(csv_file, delimiter=',')
                                for row in csv_reader:
                                    x_truth.append(float(row[0]))
                                    y_truth.append(float(row[1]))

                            [truth_inside, truth_outside] = num_src_in_ring(x_center, y_center, radius_in_pixels[i + 1],
                                                                radius_in_pixels[i], x_truth, y_truth)
                            added_star_per_mag_bin_inside[mag_index] += 1.0 * truth_inside
                            added_star_per_mag_bin_outside[mag_index] += 1.0 * truth_outside

                            [detected_inside, detected_outside] = num_src_in_ring(x_center, y_center, radius_in_pixels[i + 1], \
                                radius_in_pixels[i], x_detected[noise_index][mag_index][step], y_detected[noise_index][mag_index][step])
                            detected_star_per_mag_bin_inside[mag_index] += 1.0 * detected_inside
                            detected_star_per_mag_bin_outside[mag_index] += 1.0 * detected_outside
                            # for step
                        mag_index += 1
                        # for mag
                    # for noise_index

                # For each magnitude bin, across all noise_index and num_steps, calculate detection ratio:
                # detected_star_per_mag_bin_inside / added_star_per_mag_bin_inside
                ratio_detected_star_per_mag_bin_outside = detected_star_per_mag_bin_outside / added_star_per_mag_bin_outside
                ratio_detected_star_per_mag_bin_inside = detected_star_per_mag_bin_inside / added_star_per_mag_bin_inside

                poptinside, pcovinside = curve_fit(pritchet_formula, cfg.MAGNITUDES, ratio_detected_star_per_mag_bin_inside, maxfev=8000)
                poptoutside, pcovoutside = curve_fit(pritchet_formula, cfg.MAGNITUDES, ratio_detected_star_per_mag_bin_outside, maxfev=8000)

                completeness_csv_file = result_dir_name + "/" + g + "/" + g + "_" + f + "_radius_" + str(list[i + 1]) + "-" + str(list[i]) + "_completeness.csv"
                completeness_fig_file = result_dir_name + "/" + g + "/" + g + "_" + f + "_radius_" + str(list[i + 1]) + "-" + str(list[i]) + "_completeness.png"
                with open(completeness_csv_file, 'w') as csvout:
                    csv_writer = csv.writer(csvout)
                    csv_writer.writerow("inside region")
                    tmparray = np.array(cfg.MAGNITUDES)
                    tmprow = tmparray.tolist()
                    tmprow.insert(0, 'magnitude')
                    csv_writer.writerow(tmprow)
                    tmprow = detected_star_per_mag_bin_inside.tolist()
                    tmprow.insert(0, 'detected')
                    csv_writer.writerow(tmprow)
                    tmprow = added_star_per_mag_bin_inside.tolist()
                    tmprow.insert(0, 'added')
                    csv_writer.writerow(tmprow)
                    csv_writer.writerow(['alpha_f', poptinside[0]])
                    csv_writer.writerow(['t1_limit', poptinside[1]])
                    csv_writer.writerow(['factor_max', poptinside[2]])
                    csv_writer.writerow("outside region")
                    tmparray = np.array(cfg.MAGNITUDES)
                    tmprow = tmparray.tolist()
                    tmprow.insert(0, 'magnitude')
                    csv_writer.writerow(tmprow)
                    tmprow = detected_star_per_mag_bin_outside.tolist()
                    tmprow.insert(0, 'detected')
                    csv_writer.writerow(tmprow)
                    tmprow = added_star_per_mag_bin_outside.tolist()
                    tmprow.insert(0, 'added')
                    csv_writer.writerow(tmprow)
                    csv_writer.writerow(['alpha_f', poptoutside[0]])
                    csv_writer.writerow(['t1_limit', poptoutside[1]])
                    csv_writer.writerow(['factor_max', poptoutside[2]])
                plt.figure(cfg.CR_RADIUS[i])
                plot_completeness_function_all(g, cfg.MAGNITUDES, ratio_detected_star_per_mag_bin_inside, \
                                           cfg.MAGNITUDES, ratio_detected_star_per_mag_bin_outside, \
                                           poptinside, poptoutside, 'b', '*', \
                                           str(list[i+1]) + "-" + str(list[i]))
                plt.savefig(completeness_fig_file)
                plt.close('all')
                # for i
            filter_index += 1
            # for f
        galaxy_index += 1
        # for g
    end = time.time()
    print('Elapsed time of completeness test = ' + str(end - start) + ' seconds')

# main
cfg.load_sample_gal()
filter_names = cfg.FILTER_NAMES
image_dir_name = cfg.IMAGE_DIR_NAME
result_dir_name = cfg.RESULT_DIR_NAME
weight_dir_name = cfg.WEIGHT_DIR_NAME
config_file_name = "ngvs.sex"
galaxy_names = cfg.GALAXY_NAMES
num_steps = 4
RADIUS = 1  # If x1 and x2 are less than RADIUS pixels away, consider them to be the same source
image_file_suffix = "modsub"

if len(sys.argv) == 4:
    image_dir_name = sys.argv[1]
    result_dir_name = sys.argv[2]
    weight_dir_name = sys.argv[3]

run_sextractor_for_all()
completeness_test()
