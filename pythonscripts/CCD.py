import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
import astropy.wcs as wcs
import math
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
###%matplotlib inline
from pydl.pydlutils.spheregroup import spherematch
import copy
import csv
import os,sys
import time
from scipy.optimize import curve_fit

import config as cfg

#Here's what each section of the code does:

#run_sextractor and run_sextractor_for_all are yours

#create_table:
#     1. opens the sextractor catalogs for the u,g,i,z bands
#     2. converts each catalog into an astropy-readable Table() object
#     3. combines the data from each Table() object into a 'Source' table with all of the data for each SExtractor source
#     4. writes two copies of this source table for later manipulation: "copycat" will be manipulated and "clustercat" will eventually contain the data for all GCs
#     5. closes all of the opened SExtractor catalogs to save RAM

#create_table_for_all
#     This code iterates over every galaxy given to create a "Source" table for each. 

#analyze_galaxy
#This is the big function that does all of the heavy lifting!
#---------Initial Logistics Stuff--------------
#     1. opens the Source Tables that are made in the create_table function as astropy Table() objects
#     2. creates the 'ignore' list. This will allow me to flag all of the sources that are *not* GCs later on
#     3. creates a list of all of the coordinates of each light source using Spherematch. 
#----------Aperture Correction-------------------
#     1. compares the brightness of each source for every color filter, for both 8px and 4px apertures. This step eliminates sources that are flagged by NGVS or are not point sources.
#     2. determines the median difference in brightness and applies this value to every source that SExtractor picked up
#----------Point Sources-------------------
#     1. determines the concentration index (I4-I8) and flags every source outside this range in IC with the [IGNORE] parameter.
#     2. make a list of all of the g_mag for each source and [IGNORE] the ones that are dimmer than g=24.5
#----------Extinction Correction-------------------
#     Subtracts the extinction values from NGVS for each color filter.
#----------Graphing-------------------
#     1. creates two lists: the G-I difference for each source and the I-Z difference. This is mostly just a test and my attempt to create a GI-IZ color-color diagram
#     2. sets global fonts to ones that look pretty for the graphs later
#     3. Here's where the two different "Source" tables come into play. The "copy_cat" table is the dataset I will use to graph the color-color diagrams. The "s_cat" will enevtually be shortened to a list of all GC candidates *after* graphing.
#     4. calculate the differences of U-G and G-I for the axes of the color-color diagram (and I-Z for testing). This is done in the 8px aperture as I've found that graph looks better.
#     5. uses matplotlib.path to construct the polygon that is given to us in the Lim et al. 2017 paper. 
#     6. checks to see of the (U-G, G-I) coordinates for each light source falls within this polygon. Matplotlib.path is a little weird, so the output of the Path() function is an array filled with boolean values for each source. 
#     7. makes/updates 4 lists for graphing: it makes 2 lists for the x and y coordinates of all GC candidates, and removes the coordinates from the list with all other point sources.
#     8. plots a graph of g-mag vs. concentration index with all point sources and the boundaries on IC
#     9. plots the U-G, G-I color-color diagram. both the point source and CCD plots are saved to file, which you will have to change.
#----------List of GCs-------------------
#     1. Checks to see which points are inside the Lim et al. polygon and flags to [IGNORE] all other point sources
#     2. appends the [ignore] list to the "Source" table as a column
#     3. removes all rows where the [IGNORE] flag is true
#     4. creates a region file for DS9 with the coordinates of every GC candidate and saves that to file
#     5. writes the table of GC candidates to file (as a FITS table)
#----------More Graphing-------------------
#     This is my attempt to create a G-I, I-Z color-color diagram, with limited results.
#----------Tie it Up-------------------
#     Closes all open files to save RAM


#analyze_all_galaxies
#     This iterates through every galaxy and performs the analyze_galaxy function

#do_everything
#     This function runs all other functions: runs SExtractor, then Creates Tables, then Analyzes Galaxies

# Return number points described by (x_data, y_data) within radius of (x_center, y_center)
def get_points_in_ring(x_center, y_center, radius_inner, radius_outer, x_data, y_data):
    radius_inner_square = radius_inner * radius_inner
    radius_outer_square = radius_outer * radius_outer
    indices = []
    for i in range(0, len(x_data)):
        dist_square = (x_data[i] - x_center)*(x_data[i] - x_center) + (y_data[i] - y_center)*(y_data[i] - y_center)
        if dist_square <= radius_outer_square and dist_square >= radius_inner_square:
            indices.append(i)
    return [indices, len(x_data)-len(indices)]

#cd to dir where python script is run
#os.chdir("./")
#Credit to Nicole for the Source Extractor code!
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
        image_file1 = base_file_name + "g_" + image_file_suffix + "0.fits"          # Always use green filter image for source detection
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

            weight_file2 = weight_dir_name + "/" + g + "/" + g + "_" + f + "_sig.fits"        # Corresponding weight file
            catalog_file =  g + "_" + f + "_" + image_file_suffix + "_cat.fits"
            run_sextractor(catalog_file, image_file1, image_file2, config_file_name, weight_file1, weight_file2)

            tmp_dir_name = result_dir_name + "/" + g + "/"
            if not os.path.exists(tmp_dir_name):
                os.mkdir(tmp_dir_name)
            command = "mv " + catalog_file + " " + tmp_dir_name
            os.system(command)
                    
def create_table(u_file, g_file, i_file, z_file, galaxy_name):
    INDEX_4 = 1
    INDEX_8 = 5
    #open the FITS data filesg (your SE catalog and the FITS catalog)
    u_cat = fits.open(u_file)
    g_cat = fits.open(g_file)
    i_cat = fits.open(i_file)
    z_cat = fits.open(z_file)

    #Convert the FITS files into astropy data tables
    u_data = Table(u_cat[1].data)
    g_data = Table(g_cat[1].data)
    i_data = Table(i_cat[1].data)
    z_data = Table(z_cat[1].data)

    #Create a new 'Master' table with all of the data.
    #Format:
    #NUMBER, RA, DEC, U4, U8, G4, G8, I4, I8
    numbers = np.array(u_data['NUMBER'], dtype='f')
    source_table = Table([numbers, u_data['ALPHA_J2000'], u_data['DELTA_J2000'],
                          u_data['MAG_APER'][:,INDEX_4], u_data['MAG_APER'][:,INDEX_8],
                          g_data['MAG_APER'][:,INDEX_4], g_data['MAG_APER'][:,INDEX_8],
                          i_data['MAG_APER'][:,INDEX_4], i_data['MAG_APER'][:,INDEX_8],
                          z_data['MAG_APER'][:,INDEX_4], z_data['MAG_APER'][:,INDEX_8],
                          g_data['FLAGS']],
                         names=('NUMBER','RA','DEC','U4', 'U8','G4','G8', 'I4','I8', 'Z4', 'Z8', 'FLAGS'))
    source_table.write(result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_copycat.fits",
                       format='fits', overwrite=True)
    source_table.write(result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_clustercat.fits",
                       format='fits', overwrite=True)
    #Close the Files
    u_cat.close()
    g_cat.close()
    i_cat.close()
    z_cat.close()

def create_table_for_all():
    for g in galaxy_names:
        u_file = result_dir_name + "/" + g + "/" + g + "_u_" + image_file_suffix + "_cat.fits"
        g_file = result_dir_name + "/" + g + "/" + g + "_g_" + image_file_suffix + "_cat.fits"
        i_file = result_dir_name + "/" + g + "/" + g + "_i_" + image_file_suffix + "_cat.fits"
        z_file = result_dir_name + "/" + g + "/" + g + "_z_" + image_file_suffix + "_cat.fits"
        
        print("Creating table for " + g)
        create_table(u_file, g_file, i_file, z_file, g)
        
def analyze_galaxy(sources_file, cat_file, copy_file, galaxy_name, galaxy_index):
    """
    #This code 
    #** opens a table (from above)
    #** calibrates magnitudes based on an input catalog
    #** corrects for extinction
    #** sorts for point sources using concentration indices
    #** creates two color-color diagrams: UGI and GIZ
    #** creates a point-source diagram
    #** selects GC candidates within a certain phase-space on the UGI CCD
    #** outputs CCD png files
    #** outputs a FITS table with a list of GC candidates and their information
    #** outputs a region file with the GCs
    """
    ##################################################
    #Opening Tables###################################
    #print("Loading tables...")
    #start = time.time()
    
    #open the FITS data files (your SE source table and the NGVS "true" FITS catalog)
    cat_cat = fits.open(cat_file)
    s_cat = fits.open(sources_file)
    copy_cat = fits.open(copy_file)
    
    #Convert the FITS files into astropy data tables
    cat_data = Table(cat_cat[1].data)
    s_data = Table(s_cat[1].data)
    copy_dat = Table(copy_cat[1].data)
    
    #Ignore list -- this will be turned into a column later that will tell me when to remove sources
    ignore = [False] * len(s_data)

    #List of magnitude from catalog -- will be attached to s_data and used for GCLF later
    cat_mag = [0.0] * len(s_data)

    #end = time.time()
    #print("Time taken: " + str(end - start))
    ##################################################
    #Coordinate Matching##############################
    #print("Matching coordinates...")
    #start = time.time()
    
    #Define the tuple "coords" as the output of the spherematch function (which matches your coordinates). 
    #coords[0] = Matching SExtractor coordinate indices (the index of each matched point in the original SExtractor table)
    #coords[1] = The indices of each of the matched objects in the FITS catalog
    #coords[2] = How close the two matched coordinates are in degrees
    #Since I'm using the same detection image in SE (g), the objects it finds are the same across filters
    coords = ()
    #Column names are in capital letters sometime, but it's not guaranteed
    ra_name = 'RA'
    dec_name = 'DEC'
    gmag_name = 'GMAG'
    ic_name = 'IC'
    flags_name = 'FLAGS'
    names_in_capital = 1
    column_names = cat_data.columns
    if not ra_name in column_names:
        ra_name = 'ra'
        dec_name = 'dec'
        gmag_name = 'gmag'
        ic_name = 'iC'
        flags_name = 'flags'
        names_in_capital = 0

    coords = spherematch(s_data['RA'], s_data['DEC'], cat_data[ra_name], cat_data[dec_name], matchlength=0.0001, chunksize=0.0005)
    for i in range(len(coords[0])):
        cat_mag[coords[0][i]] = cat_data[gmag_name][coords[1][i]]

    #end = time.time()
    #print("Time taken: " + str(end - start))
    ##################################################
    #Offsets##########################################
    #print("Correcting magnitudes...")
    #start = time.time()
    for f in filter_names:
        #8px Magnitude unpacking and correcting

        #List all of the magnitudes with an 8px aperture
        f_mag_8 = []
        cat_mag_f = []
        for i in range(len(coords[0])):
            if cat_data[flags_name][coords[1][i]][1] == 0 and cfg.CI_LOWER_BOUND < cat_data[ic_name][coords[1][i]] < cfg.CI_UPPER_BOUND \
                    and s_data['FLAGS'][coords[0][i]] == 0:
                if names_in_capital:
                    cat_mag_f.append(cat_data[f.upper() + "MAG"][coords[1][i]])
                else:
                    cat_mag_f.append(cat_data[f + "mag"][coords[1][i]])
                f_mag_8.append(s_data[f.upper() + "8"][coords[0][i]])
                
        #Calculate the median magnitude offset
        diff_8f = list(map(lambda x, y: x - y, cat_mag_f, f_mag_8))
        offset_8f = np.median(diff_8f)
            
        ###########################################################
        #Now, do it all again for the 4px magnitudes

        #List all of the magnitudes with an 4px aperture
        f_mag_4 = []
        for i in range(len(coords[0])):
            if cat_data[flags_name][coords[1][i]][1] == 0 and cfg.CI_LOWER_BOUND < cat_data[ic_name][coords[1][i]] < cfg.CI_UPPER_BOUND \
                    and s_data['FLAGS'][coords[0][i]] == 0:
                f_mag_4.append(s_data[f.upper() + "4"][coords[0][i]])

        #Calculate the median magnitude offset 
        diff_4f = list(map(lambda x, y: x - y, cat_mag_f, f_mag_4))
        offset_4f = np.median(diff_4f)
        #Correct the tables 
        s_data[f.upper() + "8"] = s_data[f.upper() + "8"] + offset_8f
        s_data[f.upper() + "4"] = s_data[f.upper() + "4"] + offset_4f
        copy_dat[f.upper() + "8"] = copy_dat[f.upper() + "8"] + offset_8f
        copy_dat[f.upper() + "4"] = copy_dat[f.upper() + "4"] + offset_4f
    ###########################################################
    # Extinction Correction
    ###########################################################
    # Format: (Au, Ag, Ai, Az)
    for f in filter_names:
        s_data[f.upper() + "8"] = s_data[f.upper() + "8"] - cfg.extinction_magnitudes[galaxy_name][
            filter_names.index(f)]
        s_data[f.upper() + "4"] = s_data[f.upper() + "4"] - cfg.extinction_magnitudes[galaxy_name][
            filter_names.index(f)]

        copy_dat[f.upper() + "8"] = copy_dat[f.upper() + "8"] - cfg.extinction_magnitudes[galaxy_name][
            filter_names.index(f)]
        copy_dat[f.upper() + "4"] = copy_dat[f.upper() + "4"] - cfg.extinction_magnitudes[galaxy_name][
            filter_names.index(f)]

    ##################################################
    #Concentration Indices
    ##################################################
    #ONLY NEED I FILTER

    #Calculate concentration indices, find point source, and append their indices to a list.
    #Also remove all non-ps's from the s_data table
    point_indices = []
    ci_list = []
    for i in range(len(s_data)):
        ci = s_data['I4'][i] - s_data['I8'][i]
        if cfg.CI_LOWER_BOUND < ci < cfg.CI_UPPER_BOUND:
            point_indices.append(i)
        ci_list.append(ci)
    #Ignore points with a bad I4 - I8 (CI)
    bad_ci_count = 0
    bad_mag_count = 0
    bad_ci_mag_count = 0
    for i in range(len(s_data)):
        ci = s_data['I4'][i] - s_data['I8'][i]
        if ci > cfg.CI_UPPER_BOUND or ci < cfg.CI_LOWER_BOUND or s_data['G4'][i] > cfg.MAG_THRESHOLD:
            bad_ci_mag_count += 1
            if ci > cfg.CI_UPPER_BOUND or ci < cfg.CI_LOWER_BOUND:
                bad_ci_count += 1
            if s_data['G4'][i] > cfg.MAG_THRESHOLD:
                bad_mag_count += 1
            ignore[i] = True
            
    gmag = []
    for i in range(len(s_data)):
        gmag.append(s_data['G4'][i])

    ###########################################################
    #Throw out all of the 'bad' data points in the table. 
    ###########################################################
    #Ignore all points with flags or a bad IC
    #For source i in NGVS catalog that doesn't qualify and if it
    #matched with a source in SExtractor catalog, identify the index
    #of i in the matched list coords[1] and save it in coordinate_index
    #Get the corresponding SExtractor source index and mark ignore[s_index]
    for i in range(len(cat_data)):
        if cat_data[flags_name][i][1] != 0 or cfg.CI_LOWER_BOUND > cat_data[ic_name][i] \
                or cat_data[ic_name][i] > cfg.CI_UPPER_BOUND or s_data['FLAGS'][i] != 0:
            if i in list(coords[1]):
                coordinate_index = list(coords[1]).index(i)
                s_index = coords[0][coordinate_index]
                ignore[s_index] = True
                
    #Ignore all points with a G magnitude dimmer than cfg.MAG_THRESHOLD
    for i in range(len(s_data)):
        if s_data['G8'][i] > cfg.MAG_THRESHOLD:
            ignore[i] = True
        elif s_data['G4'][i] > cfg.MAG_THRESHOLD:
            ignore[i] = True
        else:
            pass

    #end = time.time()
    #print("Time taken: " + str(end - start))
    ###########################################################
    #Graphing (yay!)
    ###########################################################
    #print("Graphing CCDs and PS graph...")
    #start = time.time()
    
    #All CCD GIZ points prior to fitting a polygon
    g8 = np.array(s_data['G8'])
    i8 = np.array(s_data['I8'])
    z8 = np.array(s_data['Z8'])
    gi_b = g8 - i8 #Used in G-I, I-Z graph
    iz_b = i8 - z8 #Used in G-I, I-Z graph

    #Set global fonts
    font = {'family' : 'Times New Roman'}
    rc('font', **font)  # pass in the font dict as kwargs
    
    #Use data from the duplicate s_data table to graph all of the points before making a gc catalog
    #Extract rows in copy_dat where IGNORE is False
    ignorecol_copy = Column(ignore, name='IGNORE')
    copy_dat.add_column(ignorecol_copy)
    tmp_copy_dat = Table(dtype=copy_dat.dtype)
    for i in range(len(copy_dat)):
        if copy_dat['IGNORE'][i] != True:
            tmp_copy_dat.add_row(copy_dat[i])
    copy_dat = tmp_copy_dat

    #Calculate subtractions for the graphs, and make lists for matplotlib to use
    #u-g (8px)
    u_g_8 = []
    for i in range(len(copy_dat)):
        diff = copy_dat['U8'][i] - copy_dat['G8'][i]
        u_g_8.append(diff)
    #g-i (8px)
    g_i_8 = []
    for i in range(len(copy_dat)):
        diff = copy_dat['G8'][i] - copy_dat['I8'][i]
        g_i_8.append(diff)
    ###################################
    #i-z (8px)
    i_z_8 = []
    for i in range(len(copy_dat)):
        diff = copy_dat['I8'][i] - copy_dat['Z8'][i]
        i_z_8.append(diff)

    #Select points inside the Lim+2017 polygon using matplotlib.path
    verts = [
        (0.8, 0.5),  # far left, top
        (1.1, 0.5),  # left, top
        (1.9, 1.0),  # right, top
        (1.9, 1.2),  # far right, bottom
        (1.5, 1.2),  # right, bottom
        (0.8, 0.8),  # left, bottom
        (0.8, 0.5)   # ignored
    ]
    
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    path = Path(verts, codes)
    
    #Create a list of x,y coord tuples to check if inside polygon
    coords8 = []
    for i in range(len(u_g_8)):
        coords8.append((u_g_8[i], g_i_8[i]))

    #Check if coords are inside the polygon. This outputs an array for each tuple: [[True/False], ...]
    #If inside polygon, remove the coords from the u_g_x and g_i_x lists
    tf_array8 = path.contains_points(coords8)

    gc_8_x = []
    gc_8_y = []
    gc_index = []
   
    for i in range(len(coords8)):
        if tf_array8[i] == True:
            gc_8_x.append(coords8[i][0])
            gc_8_y.append(coords8[i][1])
            u_g_8.remove(coords8[i][0])
            g_i_8.remove(coords8[i][1])
        else:pass
    
    #U-G, G-I Graph (8px)
    fig2 = plt.figure(2, figsize = (4, 4), dpi=200)
    ccd8 = fig2.add_subplot(111)
    #Plot points and polygon
    ccd8.scatter(x=u_g_8, y=g_i_8, c='black', s=0.5)
    ccd8.scatter(x=gc_8_x, y=gc_8_y, c='r', s=0.5)
    patch8 = patches.PathPatch(path, ec='blue', fill=False, facecolor=None, lw=0.5, ls='--')
    ccd8.add_patch(patch8)
    #Set tick marks
    ccd8.xaxis.set_minor_locator(AutoMinorLocator())
    ccd8.yaxis.set_minor_locator(AutoMinorLocator())
    ccd8.tick_params(axis="x", direction="in")
    ccd8.tick_params(which="minor", axis="x", direction="in")
    ccd8.tick_params(axis="y", direction="in")
    ccd8.tick_params(which="minor", axis="y", direction="in")
    ccd8.tick_params(bottom=True, top=True, left=True, right=True)
    ccd8.tick_params(which="minor", bottom=True, top=True, left=True, right=True)
    #Axis labels and limits
    plt.xlabel(r'$(u^* - g\prime)_0$')
    plt.ylabel(r'$(g\prime - i\prime)_0$')
    plt.title(galaxy_name + " Color-Color Diagram (8px Aperture)")
    plt.ylim(-0.5, 2)
    plt.xlim(0, 2.5)
    plt.gca().invert_yaxis()
    #Save and show plot
    fig2.savefig(result_dir_name + "/" + galaxy_name + "/CCD_8px.png", dpi=300, bbox_inches='tight')
    plt.close('all')

    #Point Source graph
    fig3 = plt.figure(3, figsize = (5, 5), dpi=200)
    ps = fig3.add_subplot(111)
    ps.scatter(x=ci_list, y=gmag, c='black', s=0.5)
    plt.axvline(x=cfg.CI_LOWER_BOUND, c='grey', ls=":", alpha=0.5)
    plt.axvline(x=cfg.CI_UPPER_BOUND, c='grey', ls=":", alpha=0.5)
    ps.xaxis.set_minor_locator(AutoMinorLocator())
    ps.yaxis.set_minor_locator(AutoMinorLocator())
    ps.tick_params(axis="x", direction="in")
    ps.tick_params(which="minor", axis="x", direction="in")
    ps.tick_params(axis="y", direction="in")
    ps.tick_params(which="minor", axis="y", direction="in")
    ps.tick_params(bottom=True, top=True, left=True, right=True)
    ps.tick_params(which="minor", bottom=True, top=True, left=True, right=True)
    plt.xlabel(r"$i_4 - i_8$ [mag]")
    plt.ylabel("g [mag]")
    plt.title(galaxy_name + " Point Sources")
    plt.ylim(17, 28)
    plt.xlim(-.4, .4)
    plt.gca().invert_yaxis()
    fig3.savefig(result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "pointsources.png", dpi=200, bbox_inches='tight')
    plt.close('all')

    #end = time.time()
    #print("Time taken: " + str(end - start))
    ###################################################################
    #Finally, the most important step, return a table of GC candidates!
    ###################################################################
    #print("Saving GC catalog...")
    #start = time.time()
    
    for i in range(len(s_data)):
        diffx = s_data['U8'][i] - s_data['G8'][i]
        diffy = s_data['G8'][i] - s_data['I8'][i]
        if path.contains_points([(diffx, diffy)]) == False:
            ignore[i] = True
            
    ignorecol = Column(ignore, name='IGNORE')
    s_data.add_column(ignorecol)
    #Save G magnitude from catalog for GCLF. Use G8 instead.
    #magcol = Column(cat_mag, name='CATALOG_GMAG')
    #s_data.add_column(magcol)
    tmp_s_data = Table(dtype=s_data.dtype)
    for i in range(len(s_data)):
        if s_data['IGNORE'][i] != True:
            tmp_s_data.add_row(s_data[i])
    s_data = tmp_s_data

    #Create a region file for DS9
    r = open(result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_gcs.reg", "w")
    for i in range(len(s_data['NUMBER'])):
        r.write("J2000; box " + str(s_data['RA'][i]) + "d " + str(s_data['DEC'][i]) + "d " + "2\" " + "2\" " + "45 " + "# " + "color=cyan" + "\n")
    #Draw two circles in region file: GCs inside the inner circle is counted as GCs of the galaxy.
    #The outer circle and the inner circle form an annulus for background GC surface density calculation.
    base_file_name = image_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_g_" + image_file_suffix
    image_file = base_file_name + "0.fits"
    suffix = 9
    while not os.path.exists(image_file):
        image_file = base_file_name + str(suffix) + ".fits"
        suffix -= 1
        if suffix < 1:
            break
    image = fits.open(image_file)
    image_data = image[0].data
    image_size = len(image_data[1])
    w = wcs.WCS(image_file)
    x_center, y_center = w.all_world2pix(cfg.GALAXY_CENTER_RA[galaxy_index], cfg.GALAXY_CENTER_DEC[galaxy_index], 1)
    re_pixel = cfg.GALAXY_EFFECTIVE_RADII[galaxy_index] / cfg.ARCSEC_PER_PIXEL
    annuli_radius = []
    for factor in cfg.GALAXY_RE_FACTOR:
        tmp_radius = factor * re_pixel
        annuli_radius.append(tmp_radius)
        r.write("IMAGE; circle " + str(x_center) + ' ' + str(y_center) + ' ' + str(tmp_radius) + ' # color=green\n')
    outer_radius = max(max(x_center, image_size-1-x_center), max(y_center, image_size-1-y_center)) - 100
    annuli_radius.append(outer_radius)
    r.write("IMAGE; circle " + str(x_center) + ' ' + str(y_center) + ' ' + str(outer_radius) + ' # color=green\n')
    r.close()
    
    #Write table of GC candidates to file        
    s_data.write(result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_gc_cat.fits",
                 format='fits', overwrite=True)
    #s_data.show_in_browser(jsviewer=True)
    
    #end = time.time()
    #print("Time taken: " + str(end - start))

    ######################################################################
    #Graphing pt. 2
    #This is my attempt to create a G-I, I-Z color-color diagram, with limited results.
    ######################################################################
    # gi = []
    # iz = []
    # for i in range(len(s_data)):
    #     diff_gi = s_data['G8'][i] - s_data['I8'][i]
    #     diff_iz = s_data['I8'][i] - s_data['Z8'][i]
    #     gi.append(diff_gi)
    #     iz.append(diff_iz)
    #
    # #G-I, I-Z Graph (8px)
    # figgi = plt.figure(2, figsize = (4, 4), dpi=200)
    # ccdgi = figgi.add_subplot(111)
    # #Plot points and polygon
    # ccdgi.scatter(x=gi_b, y=iz_b, c='black', s=0.5)
    # ccdgi.scatter(x=gi, y=iz, c='r', s=0.5)
    # #Set tick marks
    # ccdgi.xaxis.set_minor_locator(AutoMinorLocator())
    # ccdgi.yaxis.set_minor_locator(AutoMinorLocator())
    # ccdgi.tick_params(axis="x", direction="in")
    # ccdgi.tick_params(which="minor", axis="x", direction="in")
    # ccdgi.tick_params(axis="y", direction="in")
    # ccdgi.tick_params(which="minor", axis="y", direction="in")
    # ccdgi.tick_params(bottom=True, top=True, left=True, right=True)
    # ccdgi.tick_params(which="minor", bottom=True, top=True, left=True, right=True)
    # #Axis labels and limits
    # plt.xlabel(r'$(g - i)_0$')
    # plt.ylabel(r'$(i - z)_0$')
    # plt.title(galaxy_name + " Color-Color Diagram (8px Aperture)")
    # plt.ylim(-0.5, 2)
    # plt.xlim(0, 2.5)
    # plt.gca().invert_yaxis()
    # #Save and show plot
    # figgi.savefig(result_dir_name + "CCD_8px_2.png", dpi=300, bbox_inches='tight')
    # plt.show()

    #Close the Files
    s_cat.close()
    cat_cat.close()
    copy_cat.close()

#Plot the surface density of GCs as a function of radius
def plot_surface_density(galaxy_name, galaxy_index, annuli_radius, annuli_radius_arcsec, surface_density, error, radius_arcsec_values, sersic_values):
    plt.figure(galaxy_index)
    color = "black"

    modified_surface_density = np.array(surface_density)
    modified_error = np.array(error)

    plt.yscale('log', nonposy='clip')
    plt.title('Galaxy ' + galaxy_name + '\nSurface Density of GCs', fontsize=14)
    plt.xlabel('Radius   (arcsec)')
    plt.ylabel('Surface Density  (#GCs/(arcmin^2)')

    plt.scatter(x=annuli_radius_arcsec, y=modified_surface_density, c='r', marker='*', alpha=0.5)
    plt.axvline(x=annuli_radius[1] * cfg.ARCSEC_PER_PIXEL,
                label='Radius (1.5 x Re,galaxy)'.format(annuli_radius[1] * cfg.ARCSEC_PER_PIXEL),
                color='green', linestyle='--', alpha=0.3)
    plt.errorbar(x=annuli_radius_arcsec, y=modified_surface_density, yerr=modified_error, color='blue', capsize=4, ls='none')

    #Convert radius_values from pixel to arcsec
    plt.plot(radius_arcsec_values, sersic_values, 'm--', label='Sersic profile', alpha=0.5)
    plt.legend()

    plt.savefig(result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_surface_density.png", dpi=300, bbox_inches='tight')
    plt.close('all')

#Plot Pritchet formular and Gaussian distribution
def plot_pritchet_gaussiandef(galaxy_name, r1, r2, x_data, y1_data, y2_data):
    color = "black"
    plt.ylim(0.0, 1.0)
    plt.xticks(np.arange(18, 28.5, 0.5))
    plt.yticks(np.arange(0.0, 1.0, 0.1))
#    plt.grid(color = 'grey', linestyle = '-.', linewidth = 1)
    plt.title('Galaxy ' + galaxy_name + ' Pritchet formula and Expected GCLF function\n'\
              '(annulus=' + str(r1) + '-' + str(r2) +' Re)', fontsize=14)
    plt.xlabel("Magnitudes")
#    plt.ylabel("Detection ratio")

    plt.scatter(x=x_data, y=y1_data, c='r', s=1, label='Completeness function', alpha=0.5)
    plt.scatter(x=x_data, y=y2_data, c='b', s=1, label='Expected GCLF function', alpha=0.5)
    multiplication_result = y1_data * y2_data
    #Create a step function, where the value is 1.0 when magnitude <= 24.5 and
    #0.0 when magnitude > 24.5. Multiply the step function with multiplication_result
    #This effectively set completeness function value to 0.0 where magnitude > 24.5
    step_function = np.where(x_data<=24.5, 1.0, 0.0)
    multiplication_result *= step_function

    plt.plot(x_data, multiplication_result, c='darkgreen', label='Product of completeness and expected GCLF function', alpha=0.5)
    plt.legend(loc='lower left', bbox_to_anchor=(0, -0.4))
    plt.subplots_adjust(bottom=0.25)

#Given Prichet curve described by alpha_f, t1_limit, foctor_max, calculate value of given magnitude
def pritchet_formula(mag, alpha_f, t1_limit, factor_max):
    return 0.5 * (1 - (alpha_f*(mag - t1_limit)/np.sqrt(1 + alpha_f * alpha_f * (mag - t1_limit) * (mag - t1_limit)))) * factor_max

def gaussian_dist(mag, mu, sigma):
    return (1.0/(sigma*cfg.SQRT_2_PI))*np.exp(-0.5*(mag-mu)*(mag-mu)/(sigma*sigma))

#Sersic profile N(R) to fit the GC radial density profile
#R - radius, Re - Effective radius, Ne - Density of GC at Re, Bg - Density of GC in background
def sersic_profile(R, Re, Ne, Bg):
    n = 1.0             #Sersic index
    bn = 2 * n - 0.324
    return Ne * np.exp(-1.0*bn*(pow(R/Re, 1.0/n) - 1.0)) + Bg

#Calculate correction of the surface density in given radial bin for completeness.
#Integrating the completeness curve multiplied by the expected GCLF with
#mu_g = 23.9 mag, sigma_g = 0.8 mag.  Jordan et al. (2006), (2007)
#divided by the integrated expected GCLF
def calculate_correction(galaxy_name, filter_name, factors, start_end_indices):
    #Read saved completeness data from r1 to r2.  This could be a range
    #that consists of multiple annuli. Add up all detected GCs and added GCs
    #across the annuli and curve_fit Pritchet's formula to come up with
    #alpha_f, t1_limit, factor_max
    detected = np.zeros(len(cfg.MAGNITUDES), dtype=float)
    added = np.zeros(len(cfg.MAGNITUDES), dtype=float)
    #Completeness is measured at 8.0 * Re_galaxy. Completeness result file is named so.
    #But last element in factors list is not necessarily at 8.0 * Re_galaxy.  Hack the
    #the last element so files can be named correctly.
    factors_8 = list(factors)
    factors_8[-1] = 8.0
    for i in range(start_end_indices[0], start_end_indices[1]+1):
        if i == 0:
            r1 = 0
        else:
            r1 = factors_8[i-1]
        r2 = factors_8[i]
        completeness_csv_file = result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_" + filter_name \
                                + "_radius_" + str(r1) + "-" + str(r2) + "_completeness.csv"
        with open(completeness_csv_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            completeness_header = next(csv_reader)
            completeness_magnitudes = next(csv_reader)
            completeness_detected = next(csv_reader)
            completeness_added = next(csv_reader)
            tmp_detected = np.array(completeness_detected[1:], dtype=float)
            tmp_added = np.array(completeness_added[1:], dtype=float)
        detected += tmp_detected
        added += tmp_added
    completeness_ratio = (1.0 * detected) / added
    poptinside, pcovinside = curve_fit(pritchet_formula, cfg.MAGNITUDES, completeness_ratio, maxfev=8000)
    alpha_f = poptinside[0]
    t1_limit = poptinside[1]
    factor_max = poptinside[2]

    #Take sufficiently large range in magnitude of +/- 5 sigma_g from the mean mu_g
    #Sample in finely spaced grid (n=1000) and use Simpson's rule to approximate
    #the integration of multiplication of completeness and expected GCLF
    mag_lower = cfg.MU_G - cfg.SIGMA_G * 5
    mag_upper = cfg.MU_G + cfg.SIGMA_G * 5
    mag_delta = (mag_upper - mag_lower) / cfg.NUM_MAG_DELTA
    coefficient = 0
    numerator = 0.0
    denominator = 0.0
    magnitude = mag_lower
    #Store sample magnitude values and corresponding Pritchet formular and Gaussian
    #distribution for plotting
    magnitude_values = np.zeros(cfg.NUM_MAG_DELTA+1, dtype=float)    #Store values of sample magnitudes
    pritchet_values = np.zeros(cfg.NUM_MAG_DELTA+1, dtype=float)     #Store values of Pritchet formular
    gaussian_values = np.zeros(cfg.NUM_MAG_DELTA+1, dtype=float)     #Store values of Gaussian distribution
    for i in range(0, cfg.NUM_MAG_DELTA+1):
        if i == 0 or i == cfg.NUM_MAG_DELTA:
            coefficient = 1.0
        else:
            if i % 2 == 1:
                coefficient = 4
            else:
                coefficient = 2
        tmp_p = pritchet_formula(magnitude, alpha_f, t1_limit, factor_max)
        pritchet_values[i] = tmp_p
        tmp_g = gaussian_dist(magnitude, cfg.MU_G, cfg.SIGMA_G)
        gaussian_values[i] = tmp_g
        magnitude_values[i] = magnitude
        #Sources fainter than 24.5 magnitude (cfg.MAG_THRESHOLD) are NOT considered
        #as GCs. So set completeness function value to 0 for magnitude > 24.5
        if magnitude <= 24.5:       #This is as if a step function is multiplied to it.
            numerator += coefficient * pritchet_values[i] * gaussian_values[i]
        denominator += coefficient * gaussian_values[i]
        magnitude += mag_delta
    #Simpson's rule is applied to both numerator and denominator,
    #The factor mag_delta/3 is canceled out.
    correction = numerator / denominator

    #Plot pritchet formular and gaussian distribution in the same chart.
    pritchet_gclf_png_file = result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_" + filter_name \
                    + "_radius_" + str(r1) + "-" + str(r2) + "_pritchet_gclf.png"
    plt.figure(1)
    plot_pritchet_gaussiandef(galaxy_name, r1, r2, magnitude_values, pritchet_values, gaussian_values)
    plt.savefig(pritchet_gclf_png_file)
    plt.close('all')

    return correction

#Integrate Sersic function over radius to get total number of GCs
#From radius 0 to 10 * Re, calculate number of GCs by summing up
#the number of GCs on each of the infinitely thin annuli. The number
#of GC of such annulus is sersic_function(r) * 2 * Pi * r * dr
def integrate_sersic_function(radius_arcsec_stop, Re_sersic, Ne_sersic, Bg_sersic):
    #Use Simpson's rule to approximate the integration of 2 * Pi * r * sersic_function(r)
    radius_arcsec_delta = (radius_arcsec_stop - 0.0) / cfg.NUM_MAG_DELTA
    num_gcs = 0.0
    radius_arcsec = 0.0
    coefficient = 0
    for i in range(0, cfg.NUM_MAG_DELTA+1):
        if i == 0 or i == cfg.NUM_MAG_DELTA:
            coefficient = 1
        else:
            if i % 2 == 1:
                coefficient = 4
            else:
                coefficient = 2
        tmp_s = 2.0 * math.pi * radius_arcsec * sersic_profile(radius_arcsec, Re_sersic, Ne_sersic, Bg_sersic)
        tmp_s /= 3600.00    #surface density number was in 1/arcmin^2.  Convert to 1/arcsec^2
        num_gcs += coefficient * tmp_s
        radius_arcsec += radius_arcsec_delta
    num_gcs *= (radius_arcsec_delta / 3.0)
    return num_gcs

#Some of the radial bins have 0 GC in them. Combine these bins with
#neighboring bins where there are GCs.
#new_radial_bins - a list of pairs, each pair is the start and end
#index of the original radial bins.
#If there is no GC in any of the bins, return without combining
#Any bins with zero GC will be combined with the first bin with GC
#to its right (bigger index). If there are bins without GC but aren't
#bins with GC to their right, combine them with the first bin with
#GC to their left (smaller index).
def combine_radial_bins(raw_num_gcs_in_annuli):
    new_radial_bins = []
    #Find the index of the rightmost bin with GC
    index_of_last_nonzero_bin = -1
    for i in reversed(range(len(raw_num_gcs_in_annuli))):
        if raw_num_gcs_in_annuli[i] != 0:
            index_of_last_nonzero_bin = i
            break

    if index_of_last_nonzero_bin == -1:
        # If there is no GC in any of the bins, return without combining
        for i in range(len(raw_num_gcs_in_annuli)):
            new_radial_bins.append([i, i])
    else:
        index = 0
        while index <= index_of_last_nonzero_bin:
            #Any bins with zero GC will be combined with the
            #first bin with GC to its right (bigger index).
            start_index = index
            while raw_num_gcs_in_annuli[index] == 0 and index <= index_of_last_nonzero_bin:
                index += 1
            end_index = index
            new_radial_bins.append([start_index, end_index])
            index += 1
        #If index_of_last_nonzero_bin is not the last bin, then any bin with
        #index bigger than it contains 0 GCs.  Add them to the last new_radial_bins
        if index_of_last_nonzero_bin < (len(raw_num_gcs_in_annuli) - 1):
            new_radial_bins[-1][-1] = len(raw_num_gcs_in_annuli) - 1
    return new_radial_bins

#Count number of GCs per radial bin per galaxy (Default to green filter)
#Calculate surface density and error and use completeness test data to
#calculate correction factor to surface density and error. Fit the data to
#fit Sersic profile and integrate it across radius of galaxy to predict
#total number of GCs for the galaxy and its GC system effective radius
def calculate_num_gcs(galaxy_name, galaxy_index, gcs_observed_table, surface_density_info_table,
                      error_info_table, correction_info_table, gcs_info_table):
    filter_name = 'g'
    gc_cat_file = result_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_gc_cat.fits"
    gc_cat = fits.open(gc_cat_file)
    s_data = Table(gc_cat[1].data)

    #Get image size in pixel, galaxy center coordinates in pixel, effective radius
    #in pixel, coordinates of candidate GCs in pixel
    base_file_name = image_dir_name + "/" + galaxy_name + "/" + galaxy_name + "_g_" + image_file_suffix
    image_file = base_file_name + "0.fits"
    suffix = 9
    while not os.path.exists(image_file):
        image_file = base_file_name + str(suffix) + ".fits"
        suffix -= 1
        if suffix < 1:
            break
    image = fits.open(image_file)
    image_data = image[0].data
    image_size = len(image_data[1])
    w = wcs.WCS(image_file)
    x_gcs = np.zeros(len(s_data))   #x coordinate of GCs
    y_gcs = np.zeros(len(s_data))   #y coordinate of GCs
    for i in range(len(s_data)):
        x_gcs[i], y_gcs[i] = w.all_world2pix(s_data['RA'][i], s_data['DEC'][i], 1)

    #Get galaxy center coordinates in pixel.  Get array of radius of a set of annuli
    x_center, y_center = w.all_world2pix(cfg.GALAXY_CENTER_RA[galaxy_index], cfg.GALAXY_CENTER_DEC[galaxy_index], 1)
    re_pixel = cfg.GALAXY_EFFECTIVE_RADII[galaxy_index] / cfg.ARCSEC_PER_PIXEL
    annuli_radius = []                      #In unit of pixel
    annuli_area = []                        #In unit of arcmin^2 base on Prof. Peng's suggestion
    previous_area_of_circle_arcmin = 0.0
    for factor in cfg.GALAXY_RE_FACTOR:
        tmp_radius = factor * re_pixel
        annuli_radius.append(tmp_radius)
        area_of_circle_arcmin = (math.pi * tmp_radius * cfg.ARCSEC_PER_PIXEL * tmp_radius * cfg.ARCSEC_PER_PIXEL) / (60.0 * 60.0)
        annuli_area.append(area_of_circle_arcmin - previous_area_of_circle_arcmin)
        previous_area_of_circle_arcmin = area_of_circle_arcmin
    outer_radius = max(max(x_center, image_size-1-x_center), max(y_center, image_size-1-y_center)) - 100 #Draw an outside circle that's 100 pixel smaller than image size
    annuli_radius.append(outer_radius)
    area_of_circle_arcmin = (math.pi * outer_radius * cfg.ARCSEC_PER_PIXEL * outer_radius * cfg.ARCSEC_PER_PIXEL) / (60.0 * 60.0)
    annuli_area.append(area_of_circle_arcmin - previous_area_of_circle_arcmin)

    #Count number of observed GC candidates inside each of the annuli
    #The annuli are 0 to 0.75Re, 0.75Re to 1.5Re, 1.5Re to 2.25Re,
    #2.25Re to 3.0Re, 3.0Re to 4.5Re, 4.5Re to 8.0Re.
    #These factors are defined in list cfg.GALAXY_RE_FACTOR
    tmparray = np.array(cfg.GALAXY_RE_FACTOR)
    factors = tmparray.tolist()
    factors.append(outer_radius/re_pixel)   #Calculate how many times of Re,Galaxy the outside circle is
    raw_num_gcs_in_annuli = []
    prev_radius = 0.0
    for i in range(len(factors)):
        [gcs_in_annulus, dummy] = get_points_in_ring(x_center, y_center, prev_radius, annuli_radius[i], x_gcs, y_gcs)
        raw_num_gcs_in_annuli.append(len(gcs_in_annulus))
        prev_radius = annuli_radius[i]

    #Some of the radial bins don't have GC. Merge these bins with
    #neighboring bins when calculating surface density and error.
    #Completeness test result also needs to be re-calculated base
    #on the new radial bins.
    merged_radial_bins = combine_radial_bins(raw_num_gcs_in_annuli)
    merged_annuli_radius = []
    num_gcs_in_merged_annuli = [galaxy_name]    #Number of GCs in each annulus
    surface_density = [galaxy_name]             #Number of GCs/area in each annulus
    error = [galaxy_name]                       #Error = sqrt(Number of GCs)/area in each annulus
    corrections = [galaxy_name]                 #correction of surface_density: integrating completeness and GCLF divided by integrating GCLF
    for start_end_indices in merged_radial_bins:
        num_gcs_in_annuli = 0
        merged_annuli_area = 0
        #Radial bins from start_end_indices[0] to start_end_indices[1]
        #are merged because there is only one bin has GCs. Total number
        #of GCs in these bins is the same as the number of GCs in that
        #one bin.  However, area is the sum of all annulus in these bins.
        #Also, use that one bin's radius for the surface density v.s.
        #radius chart.
        for i in range(start_end_indices[0], start_end_indices[1]+1):
            num_gcs_in_annuli += raw_num_gcs_in_annuli[i]
            merged_annuli_area += annuli_area[i]
            if raw_num_gcs_in_annuli[i] > 0:
                merged_annuli_radius.append(annuli_radius[i])
        num_gcs_in_merged_annuli.append(num_gcs_in_annuli)
        surface_density.append(num_gcs_in_annuli*1.0 / merged_annuli_area)
        error.append(math.sqrt(num_gcs_in_annuli * 1.0) / merged_annuli_area)
        correction = calculate_correction(galaxy_name, filter_name, factors, start_end_indices)
        corrections.append(correction)

    tmp_raw_num_gcs_in_annuli = [galaxy_name]
    tmp_raw_num_gcs_in_annuli += raw_num_gcs_in_annuli
    gcs_observed_table.append(tmp_raw_num_gcs_in_annuli)
    #Make correction to surface density and error
    for i in range(1, len(surface_density)):
        surface_density[i] /= corrections[i]
        error[i] /= corrections[i]
    surface_density_info_table.append(surface_density)
    error_info_table.append(error)
    correction_info_table.append(corrections)

    #curve_fit surface density and error data to Sersic profile. p0 is initial
    #value for the parameters, setting it to reasonable value rather than
    #default value of 1 helps curve_fit to converge.
    Re_sersic = 1.0
    Ne_sersic = 1.0
    Bg_sersic = 1.0
    # The radii in merged_annuli_radius are the radii of the outer side of each of the annuli
    # If we calculate surface density using these values, it would make the numbers look bigger
    # than it should be. One way to estimate is use radius that evenly split the area of an
    # annulus.
    tmp_annuli_radii = []
    r1 = 0.0
    for r2 in merged_annuli_radius:
        tmp_annuli_radii.append(math.sqrt((r1 * r1 + r2 * r2) / 2.0))
        r1 = r2
    annuli_radius_arcsec = np.array(tmp_annuli_radii)
    annuli_radius_arcsec = annuli_radius_arcsec * cfg.ARCSEC_PER_PIXEL
    if len(merged_annuli_radius) >= 4:
        popt, pcov = curve_fit(sersic_profile, annuli_radius_arcsec, surface_density[1:], sigma=error[1:],
                               absolute_sigma=True, p0=[cfg.GALAXY_EFFECTIVE_RADII[galaxy_index]*1.5, 1.0, 0.1])
        Re_sersic = popt[0]
        Ne_sersic = popt[1]
        Bg_sersic = popt[2]

    #Calculate surface density in finer interval of radius for plotting
    #Radius interval is 1/1000 (cfg.NUM_MAG_DELTA) of annuli_radius[-1] * 2
    radius_arcsec_stop = annuli_radius[-1] * cfg.ARCSEC_PER_PIXEL * 1.2
    radius_arcsec_delta = (radius_arcsec_stop - 0.0) / cfg.NUM_MAG_DELTA
    radius_arcsec_values = np.zeros(cfg.NUM_MAG_DELTA+1, dtype=float)  #Store values of sample radius
    sersic_values = np.zeros(cfg.NUM_MAG_DELTA+1, dtype=float)  #Store values of Sersic profile value
    radius_arcsec = 0.0
    for i in range(0, cfg.NUM_MAG_DELTA+1):
        tmp_s = sersic_profile(radius_arcsec, Re_sersic, Ne_sersic, Bg_sersic)
        sersic_values[i] = tmp_s
        radius_arcsec_values[i] = radius_arcsec
        radius_arcsec += radius_arcsec_delta

    #Plot surface density of GCs as a function of radius
    plot_surface_density(galaxy_name, galaxy_index, annuli_radius, annuli_radius_arcsec,
                         surface_density[1:], error[1:], radius_arcsec_values, sersic_values)

    #Integrate Sersic function over radius to get total number of GCs
    #A value of 1.0 for any Sersic parameter means the curve_fit failed.
    #In case of curve_fit failed, take the number of GC candidates within 1.5 R_e,
    #assume a background as defined by the outermost radial bin, background subtract,
    #then double the number to get an estimate of the total number
    if Re_sersic == 1.0 or Ne_sersic == 1.0 or Bg_sersic == 1.0:
        num_gcs = 0
        annuli_area_1_5_Re = 0
        #Count observed GCs in annuli within 1.5*Re_galaxy
        for i in range(len(factors)):
            if factors[i] > 1.5:
                break
            else:
                annuli_area_1_5_Re += annuli_area[i]
            num_gcs += tmp_raw_num_gcs_in_annuli[i+1]
        #Calculate GC surface density in the outermost annulus
        background_surface_density = tmp_raw_num_gcs_in_annuli[-1]/annuli_area[-1]
        #Calculate number of background GC within 1.5*Re
        num_gcs -= background_surface_density * annuli_area_1_5_Re
        num_gcs = 2 * max(0, num_gcs)
    else:
        num_gcs = integrate_sersic_function(radius_arcsec_stop, Re_sersic, Ne_sersic, Bg_sersic)
    num_gcs = int(num_gcs+0.5)

    #Calculate GC specific frequency.  S_N = num_gcs * 10^(0.4*(M_V+15))
    #M_V is the absolute V magnitude of the galaxy, M_V = V - 31.09 mag
    #where V is the observed V' magnitude, and 31.09 is called distance modulus, corresponds to 16.5Mpc distance
    #V is calculated by assuming V = g - 0.3 mag, where g is the observed g' magnitude
    V = cfg.GALAXY_GMAG[galaxy_index] - 0.3
    M_V = V - 31.09
    S_N = 1.0 * num_gcs * pow(10, 0.4 * (M_V + 15.0))

    sersic_result = [galaxy_name, num_gcs, Re_sersic, cfg.GALAXY_EFFECTIVE_RADII[galaxy_index], M_V, S_N, Bg_sersic]
    gcs_info_table.append(sersic_result)

#analyze_all_galaxies does two things:
#analyze_galaxy - gets observed GCs using color-color diagram
#calculate_num_gcs - uses observed GCs and result of completeness test to predict total number of GCs for galaxy
def analyze_all_galaxies():
    csv_header = ['Galaxy Name', '#GCs (total)']
    previous_factor = 0
    for factor in cfg.GALAXY_RE_FACTOR:
        tmpstr = '#GCs (' + str(previous_factor) + ' to ' + str(factor) + ' Re)'
        csv_header.append(tmpstr)
        previous_factor = factor
    tmpstr = '#GCs (beyond ' + str(previous_factor) + ' Re)'
    csv_header.append(tmpstr)
    gcs_header = ['Galaxy Name', '#GCs (Integraton of Sersic)', 'Effective Radius (GC System)',
                  'Effective Radius (Galaxy)', 'M_V', 'GC Specific Frequency', 'GC Background']

    gcs_observed_table = []             #Number of GCs for each annulus
    gcs_observed_table.append(csv_header)
    del csv_header[1]
    surface_density_info_table = []     #Number_of_GCs/area for each annulus
    surface_density_info_table.append(csv_header)
    error_info_table = []               #Number_of_GCs/area for each annulus
    error_info_table.append(csv_header)
    correction_info_table = []          #correction for surface density for each annulus
    correction_info_table.append(csv_header)
    gcs_info_table = []                 #Number of GCs and GCS effective radius (Sersic)
    gcs_info_table.append(gcs_header)

    galaxy_index = 0
    num_gcs_list = []
    for g in galaxy_names:
        start = time.time()
        print("Analyzing galaxy " + g)
        cat_file = catalog_dir_name + "/" + g + "/" + g + "_xD.fits"  #Suffix could be xD or NGVS
        if not os.path.exists(cat_file):
            cat_file = catalog_dir_name + "/" + g + "/" + g + "_NGVS.fits"  # Suffix could be xD or NGVS
        sources_file = result_dir_name + "/" + g + "/" + g + "_copycat.fits"
        copy_file = result_dir_name + "/" + g + "/" + g + "_clustercat.fits"
        
        analyze_galaxy(sources_file, cat_file, copy_file, g, galaxy_index)
        calculate_num_gcs(g, galaxy_index, gcs_observed_table,
                                        surface_density_info_table, error_info_table,
                                        correction_info_table, gcs_info_table)
        galaxy_index += 1
        end = time.time()
        print("Time taken: " + str(end - start))

    gcs_observed_table_file = result_dir_name + '/gcs_observed.csv'
    with open(gcs_observed_table_file, 'w') as csvout:
        csv_writer = csv.writer(csvout)
        for row in gcs_observed_table:
            csv_writer.writerow(row)
    surface_density_info_table_file = result_dir_name + '/surface_density_info.csv'
    with open(surface_density_info_table_file, 'w') as csvout:
        csv_writer = csv.writer(csvout)
        for surface_density in surface_density_info_table:
            csv_writer.writerow(surface_density)
    error_info_table_file = result_dir_name + '/surface_density_error_info.csv'
    with open(error_info_table_file, 'w') as csvout:
        csv_writer = csv.writer(csvout)
        for error in error_info_table:
            csv_writer.writerow(error)
    correction_info_table_file = result_dir_name + '/surface_density_correction_info.csv'
    with open(correction_info_table_file, 'w') as csvout:
        csv_writer = csv.writer(csvout)
        for correction in correction_info_table:
            csv_writer.writerow(correction)
    gcs_info_table_file = result_dir_name + '/gcs_info.csv'
    with open(gcs_info_table_file, 'w') as csvout:
        csv_writer = csv.writer(csvout)
        for row in gcs_info_table:
            csv_writer.writerow(row)

def do_everything():
    tot_start = time.time()
    start = time.time()
#    run_sextractor_for_all()
    end = time.time()
    print("Time taken running SExtractor: " + str(end - start) + " seconds.")

    start = time.time()
#    create_table_for_all()
    end = time.time()
    print("Time taken creating table: " + str(end - start) + " seconds.")

    start = time.time()
    analyze_all_galaxies()
    end = time.time()
    print("Time taken identifying GCs: " + str(end - start) + " seconds.")

    tot_end = time.time()
    print("Time taken in total: " + str(tot_end - tot_start) + " seconds.")

# main
cfg.load_sample_gal()
filter_names = ("u", "g", "i", "z") #Only UGIZ capability
image_dir_name = cfg.IMAGE_DIR_NAME
catalog_dir_name = cfg.IMAGE_DIR_NAME   # Catalog file of a galaxy is in its image directory
result_dir_name = cfg.RESULT_DIR_NAME
config_file_name = "ngvs_1.0.sex"
galaxy_names = cfg.GALAXY_NAMES
weight_dir_name = cfg.WEIGHT_DIR_NAME
image_file_suffix = "modsub"

do_everything()
