# Read two CSV files, both contains SN (specific frequency) and MV(galaxy luminosity in V band) information.
# The 1st file is our own data from NGVS survey, the 2nd file is from other published paper.
from astropy.io import fits
import config as cfg
import copy
import csv
import matplotlib.pyplot as plt
import numpy as np
import os,sys

def plot_sn_mv(gcs_mv, gcs_sn, other_gcs_mv, other_gcs_sn, other_gcs_sn_error, other_survey_name, ymin, ymax, png_file_name):
    plt.figure(0)
    color = "black"

    plt.title('GC Specific Frequency against Galaxy Luminosity', fontsize=14)
    plt.xlabel(r'$M_V$')
    plt.ylabel(r'$S_N$')
    plt.ylim(ymin, ymax)

    plt.errorbar(x=other_gcs_mv, y=other_gcs_sn, yerr=other_gcs_sn_error, color='grey', capsize=0, ls='none', linewidth=1)
    plt.scatter(x=other_gcs_mv, y=other_gcs_sn, c='b', marker='.', alpha=0.5, label=other_survey_name)
    plt.scatter(x=gcs_mv, y=gcs_sn, c='r', marker='*', alpha=0.5, label='NGVS')

    plt.legend(loc='upper left')
    plt.savefig(png_file_name, dpi=300, bbox_inches='tight')
    plt.close('all')

# In gcs_info_file, information for each galaxy is a row (except the header row)
# GC's MV and SN are in column mv_col_index and sn_col_index
def read_gcs_info(gcs_info_file, mv_col_index, sn_col_index):
    gcs_mv = []
    gcs_sn = []
    gcs_sn_error = []
    with open(gcs_info_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)    # Skip the header row in csv file
        for row in csv_reader:
            gcs_mv.append(float(row[mv_col_index]))
            index = row[sn_col_index].find('+-')
            if index != -1:
                gcs_sn.append(float(row[sn_col_index][0:index]))
                gcs_sn_error.append(float(row[sn_col_index][index+2:]))
            else:
                gcs_sn.append(float(row[sn_col_index]))
                gcs_sn_error.append(0.0)
    return gcs_mv, gcs_sn, gcs_sn_error

if len(sys.argv) == 4:
    gcs_info_file = sys.argv[1]
    other_gcs_info_file = sys.argv[2]

    gcs_mv, gcs_sn, gcs_sn_error = read_gcs_info(gcs_info_file, 4, 5)
    # ACSVCS (Peng 2008)
    other_gcs_mv, other_gcs_sn, other_gcs_sn_error = read_gcs_info(other_gcs_info_file, 2, 9)
    plot_sn_mv(gcs_mv, gcs_sn, other_gcs_mv, other_gcs_sn, other_gcs_sn_error, 'Peng et al. (2008); ACSVCS', -28.0, 70.0, 'sn_mv_acsvcs.png')

    # HST_WFPC2 (Miller & Lotz 2007)
    other_gcs_info_file = sys.argv[3]
    other_gcs_mv, other_gcs_sn, other_gcs_sn_error = read_gcs_info(other_gcs_info_file, 1, 3)
    plot_sn_mv(gcs_mv, gcs_sn, other_gcs_mv, other_gcs_sn, other_gcs_sn_error, 'Miller & Lotz (2007); WFPC2', -28.0, 70.0, 'sn_mv_hst_wfpc2.png')
else:
    print('Usage: python sn_plot.py <NGVS GC info csv file1> <ACSVCS GC info csv file2> <HST WFPC2 GC info csv file2>')
