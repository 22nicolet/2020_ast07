#adds artificial stars to the original image
from astropy.io import fits
import config as cfg
import copy
import csv
import matplotlib.pyplot as plt
import numpy as np
import os,sys
import time

# @precondition psf data is from a psf file, sigma_array has no negatives
# adds gaussian noise (mean = pixel value, sigma = sqrt(pixel value/gain)) to each pixel of the psf
def add_noise_to_psf(filter_index):
    mean_array = np.copy(psf_data_array[filter_index])
    sigma_array = np.copy(psf_data_array[filter_index])
    sigma_array = np.sqrt(sigma_array / int(gain))
    psf_mag_noise = np.random.normal(mean_array, sigma_array)

    #print(psf_data[filter_index])
    #plt.figure(figsize=[5,5])
    #plt.imshow(psf_data[filter_index],origin='lower',cmap='gray_r')
    #plt.show()
    #print(psf_data)
    #plt.figure(figsize=[5,5])
    #plt.imshow(psf_mag_noise[filter_index],origin='lower',cmap='gray_r')
    #plt.show()
    #print(psf_mag_noise)

    return psf_mag_noise

# adds in NUM_GRID_POINTS * NUM_GRID_POINTS stars into the image, stars are added to a square grid around the center of the PSF
def add_star_one_step(orig_image, image_data_origin, grid_space, x_step_size, y_step_size, ratio, mag, suffix, galaxy_name, filter_index, filter_name, noise_index, psf_data_noise):
    image_data_addstar = copy.deepcopy(image_data_origin)
    for i in range (0, cfg.NUM_GRID_POINTS):
        for j in range (0, cfg.NUM_GRID_POINTS):
            position_x = int(grid_space * (i + 1) + x_step_size)
            position_y = int(grid_space * (j + 1) + y_step_size)

            # PSF image is 31x31 pixels, calculate x_min/max, y_min/max so PSF's center is at position_x and position_y
            x_min = position_x - 15
            x_max = position_x + 16
            y_min = position_y - 15
            y_max = position_y + 16

            # Overlay the fake star on top of the image background
            image_data_addstar[x_min:x_max, y_min:y_max] = image_data_addstar[x_min:x_max, y_min:y_max] + psf_data_noise / ratio

    #plt.figure(figsize=[8, 8])
    #plt.imshow(image_data_addstar, origin='lower', cmap='gray_r', norm=LogNorm())
    #plt.title('After add stars on the image', fontsize=14)
    #plt.show()

    # Write out image with added stars
    new_image_file = result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_" + filter_name + "_addstar" + str(mag) + "_" + suffix + "_" + str(noise_index) + ".fits"
    hdu = fits.PrimaryHDU(image_data_addstar, orig_image[0].header)
    hdu.writeto(new_image_file, overwrite=True)

# nested loop to add stars into images (galaxy, filter, number of different noises, magnitudes)
# each magnitude bin has four steps spaced by half the grid_space (size of the image/(NUM_GRID_POINTS + 1)), the steps form a square (shifted left, right, diagonally)
def add_stars():
    for g in galaxy_names:
        base_dir_name = image_dir_name + "/" + g
        base_file_name = image_dir_name + "/" + g + "/" + g + "_"
        filter_index = 0
        for f in filter_names:
            orig_image_file = base_file_name + f + "_modsub2.fits"  # File for measurement
            orig_image = fits.open(orig_image_file)
            image_data_origin = orig_image[0].data
            size = len(image_data_origin[1])                        # size is the dimension of the square image
            grid_space = size/(cfg.NUM_GRID_POINTS + 1)                 # the number of pixels between two grid points
            for n in range(0, cfg.NUMBER_NOISES):
                # Add Gaussian noises to each pixel of the PSF
                psf_data_noise = add_noise_to_psf(filter_index)

                # Scale to some magnitude
                # Equations: m1-m2=-2.5log(F1/F2), F1/F2=10**(-0.4*(m1-m2))
                psf_mag = -2.5 * np.log10(sum(sum(psf_data_noise))) + 30
                print('galaxy='+g+' filter='+f+' noise_idx='+str(n)+' psf_mag='+str(psf_mag))

                for mag in magnitude_ranges:  # generate fake stars of different magnitude
                    ratio = 10 ** (-0.4 * (psf_mag - mag))  # calculate the flux ratio

                    add_star_one_step(orig_image, image_data_origin, grid_space, 0, 0, ratio, mag, "0", g, filter_index, f, n, psf_data_noise)
                    add_star_one_step(orig_image, image_data_origin, grid_space, grid_space/2, 0, ratio, mag, "1", g, filter_index, f, n, psf_data_noise)
                    add_star_one_step(orig_image, image_data_origin, grid_space, 0, grid_space/2, ratio, mag, "2", g, filter_index, f, n, psf_data_noise)
                    add_star_one_step(orig_image, image_data_origin, grid_space, grid_space/2,  grid_space/2, ratio, mag, "3", g, filter_index, f, n, psf_data_noise)

        #            plt.figure(figsize=[8, 8])
        #            plt.imshow(image_data_origin, origin='lower', cmap='gray_r', norm=LogNorm())
        #            plt.title('Galaxy model subtracted image', fontsize=14)
        #            plt.show()

        #            plt.figure(figsize=[8, 8])
        #            plt.imshow(image_data_addstar, origin='lower', cmap='gray_r', norm=LogNorm())
        #            plt.title('After add stars on the image', fontsize=14)
        #            plt.show()


#main
filter_names = ("g")
image_dir_name = "../images"
result_dir_name = "../results"
galaxy_names = ["VCC0940"]
magnitude_ranges = (23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28)
NUM_STARS_TO_ADD_PER_MAG = 100

# Downloading PSF files
psf_data = [0 for x in range(len(filter_names))]
psf_data_array = [0 for x in range(len(filter_names))]
filter_index = 0

# Load effective gain table from csv file from
# http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/ngvs/docs/gain.html
gain = 1
gain_table = {}
with open('ngvs_gain_table.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        gain_table[row[0]] = row[1]


for f in filter_names:
    image_key = 'NGVS-1+0.l.' + f
    image = 'image=' + image_key + '.Mg002.fits'   # You may change the fieldname, band, x&y position etc.
    x = '17366'
    y = '19409'
    psf_URL = 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/community/ngvs/NGVSpsf.pl?'+image+'&x='+x+'&y='+y
    print(psf_URL)
    os.system('wget -O psf.fits %s' %psf_URL)           # Use wget command to download the file
    time.sleep(1)                                       # Give wget some time to get the file from web

    #reads in psf, changes the negative data numbers to 0, stores in a numpy array for further processing
    hdulist = fits.open('psf.fits')
    psf_data[filter_index] = hdulist[0].data
    np_array = np.array(psf_data[filter_index]) # Convert list into a numpy array
    np_array_sign = np.sign(np_array)  # negative elements become -1, 0 stays 0, positive elements become 1
    np_array *= (np_array > 0)                  # converts negative elements to -0
    np_array *= np_array_sign                   # replaces -0 with 0 to prevent error later on
    psf_data_array[filter_index] = np.copy(np_array)
    filter_index += 1

gain = gain_table[image_key]
start = time.time()
add_stars()
end = time.time()
print('Elapsed time = '+str(end-start)+' seconds')
