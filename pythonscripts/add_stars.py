#!/usr/bin/env python
# coding: utf-8

from   astropy.io import fits
import copy
import csv
import matplotlib.pyplot as plt
from   matplotlib.colors import LogNorm
import numpy as np
import os,sys
import random
import time
def add_star_one_step(orig_image, image_data_origin, grid_space, x_step_size, y_step_size, ratio, mag, suffix, galaxy_name, filter_index, filter_name):
    image_data_addstar = copy.deepcopy(image_data_origin)
    for i in range (0, NUM_GRID_POINTS):
        for j in range (0, NUM_GRID_POINTS):
            position_x = int(grid_space * (i + 1) + x_step_size)
            position_y = int(grid_space * (j + 1) + y_step_size)

            # PSF image is 31x31 pixels, calculate x_min/max, y_min/max so PSF's center is at position_x and position_y
            x_min = position_x - 15
            x_max = position_x + 16
            y_min = position_y - 15
            y_max = position_y + 16

            # Overlay the fake star on top of the image background
            image_data_addstar[x_min:x_max, y_min:y_max] = image_data_addstar[x_min:x_max, y_min:y_max] + psf_data[filter_index] / ratio

    #plt.figure(figsize=[8, 8])
    #plt.imshow(image_data_addstar, origin='lower', cmap='gray_r', norm=LogNorm())
    #plt.title('After add stars on the image', fontsize=14)
    #plt.show()

    # Write out image with added stars
    new_image_file = result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_" + filter_name + "_addstar" + str(mag) + "_" + suffix + ".fits"
    hdu = fits.PrimaryHDU(image_data_addstar, orig_image[0].header)
    hdu.writeto(new_image_file, overwrite=True)

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
            grid_space = size/(NUM_GRID_POINTS + 1)                 # the number of pixels between two grid points
            for mag in magnitude_ranges:  # generate fake stars of different magnitude
                ratio = 10 ** (-0.4 * (psf_mag[filter_index] - mag))  # calculate the flux ratio

                add_star_one_step(orig_image, image_data_origin, grid_space, 0, 0, ratio, mag, "0", g, filter_index, f)
                add_star_one_step(orig_image, image_data_origin, grid_space, grid_space/2, 0, ratio, mag, "1", g, filter_index, f)
                add_star_one_step(orig_image, image_data_origin, grid_space, 0, grid_space/2, ratio, mag, "2", g, filter_index, f)
                add_star_one_step(orig_image, image_data_origin, grid_space, grid_space/2,  grid_space/2, ratio, mag, "3", g, filter_index, f)

    #            plt.figure(figsize=[8, 8])
    #            plt.imshow(image_data_origin, origin='lower', cmap='gray_r', norm=LogNorm())
    #            plt.title('Galaxy model subtracted image', fontsize=14)
    #            plt.show()

    #            plt.figure(figsize=[8, 8])
    #            plt.imshow(image_data_addstar, origin='lower', cmap='gray_r', norm=LogNorm())
    #            plt.title('After add stars on the image', fontsize=14)
    #            plt.show()



filter_names = ("g")
image_dir_name = "../images"
result_dir_name = "../results"
galaxy_names = ["VCC0940"]
magnitude_ranges = (23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28)
NUM_STARS_TO_ADD_PER_MAG = 100
NUM_GRID_POINTS = 30    #number of grid points along x or y

# Downloading PSF files
psf_data = [0 for x in range(len(filter_names))]
psf_mag  = [0 for x in range(len(filter_names))]
filter_index = 0
for f in filter_names:
    image='image=NGVS-1+0.l.' + f + '.Mg002.fits'   # You may change the fieldname, band, x&y position etc.
    x='17366'
    y='19409'
    psf_URL='https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/community/ngvs/NGVSpsf.pl?'+image+'&x='+x+'&y='+y
    print(psf_URL)
    os.system('wget -O psf.fits %s' %psf_URL)       # Use wget command to download the file
    time.sleep(1)                                 # Give wget some time to get the file from web

    hdulist=fits.open('psf.fits')
    psf_data[filter_index]=hdulist[0].data
    #print(psf_data[filter_index])
    #plt.figure(figsize=[3,3])
    #plt.imshow(psf_data[filter_index],origin='lower',cmap='gray_r')
    #plt.show()

    # Scale to some magnitude
    # Equations: m1-m2=-2.5log(F1/F2), F1/F2=10**(-0.4*(m1-m2))
    psf_mag[filter_index]=-2.5*np.log10(sum(sum(psf_data[filter_index])))+30
    print(psf_mag[filter_index])

    filter_index += 1

add_stars()
