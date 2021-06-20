import csv

NUMBER_NOISES = 1       #number of different noises to add to the psf
FILTER_NAMES = ['g']
FILTER_NAMES_UGI = ['g', 'u', 'i']
FILTER_MAGS = [23.0, 24.0, 22.3]
PIXEL_BETWEEN_GCS = 38
MAG_THRESHOLD=24.5
MAGNITUDES = (23.0, 23.5, 24.0, 24.5, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5, 28.0)
NUM_STEPS = 4
ARCSEC_PER_PIXEL = 0.187
IMAGE_DIR_NAME = "/Volumes/SanDiskMac/images"
RESULT_DIR_NAME = "/Volumes/SanDiskMac/results"
WEIGHT_DIR_NAME = "/Volumes/SanDiskMac/images"
SIGMA_G = 0.8
MU_G = 23.9
NUM_MAG_DELTA = 1000    #number of delta along magnitude axis for Simpson' rule approximation
SQRT_2_PI = 2.5066283
GALAXY_RE_FACTOR = [0.75, 1.50, 2.25, 3.00, 4.50]
CR_RADIUS = [4.50, 3.00, 2.25, 1.50, 0.75, 0]         #central region radius in arcseconds
CI_UPPER_BOUND = 0.15
CI_LOWER_BOUND = -0.1

GALAXY_NAMES = []
GALAXY_CENTER_RA = []
GALAXY_CENTER_DEC = []
GALAXY_GMAG = []
GALAXY_EFFECTIVE_RADII = []
GALAXY_FIELD_NAMES = []
GALAXY_PSF_X = []
GALAXY_PSF_Y = []
extinction_magnitudes = {}

# Incomplete data
# VCC0093, there is no modsub image for g band
# VCC0162, there is no modsub image for u band needed by CCD.py
# VCC0373, error while reading VCC0373_z_sig.fits
# VCC0047, curve_fit failed
# VCC0170, curve_fit failed
# VCC0180, plotting failed
# VCC0263, plotting failed
# VCC0346, there is no modsub image for g band
# VCC0349, plotting failed
# VCC0384, plotting failed
# VCC0389, plotting failed
# VCC0143, 0 to 0.75 Re result file missing
# VCC0428, prichet fitting failed.
# VCC0641, curve_fit failed
# VCC0657, there is no modsub image for g band
# VCC0678, curve_fit failed
# VCC0697, curve_fit failed
# VCC0708, curve_fit failed
# VCC0713, create_table failed, index out of range
# VCC0733, curve_fit failed
# VCC0764, curve_fit failed
# VCC0811, prichet fitting failed.
# VCC0821, prichet fitting failed.
# VCC0862, prichet fitting failed.
# VCC0890, prichet fitting failed.
# VCC0892, prichet fitting failed.
# VCC0898, creating table failed.
# VCC0917, creating table failed.
# VCC0925, curve_fit failed
# VCC0975, curve_fit failed
# VCC0977, prichet fitting failed.
# VCC0997, curve_fit failed
# VCC1042, curve_fit failed
# VCC1055, prichet fitting failed.
# VCC1059, curve_fit failed
# VCC1081, prichet fitting failed.
# VCC1078, prichet fitting failed.
# VCC1100, prichet fitting failed.
# VCC1123, curve_fit failed
# VCC1141, prichet fitting failed.
# VCC1143, curve_fit failed
# VCC1148, curve_fit failed
# VCC1159, add_stars failed
# VCC1175, curve_fit failed
# VCC1188, prichet fitting failed.
# VCC1189, curve_fit failed
# VCC1192, prichet fitting failed
# VCC1236, curve_fit failed
# VCC1239, curve_fit failed
# VCC1249, curve_fit failed
# VCC1297, prichet fitting failed.
# VCC1301, create_table failed
# VCC1313, prichet fitting failed.
# VCC1340, curve_fit failed
# VCC1369, create_table failed
# VCC1317, create_table failed
# VCC1325, curve_fit failed
# VCC1343, curve_fit failed
# VCC1354, curve_fit failed
# VCC1355, curve_fit failed
# VCC1356, curve_fit failed
# VCC1389, there is no modsub image for z band
# VCC1403, there is no modsub image for g band
# VCC1425, prichet fitting failed.
# VCC1435, create_table failed
# VCC1436, curve_fit failed
# VCC1437, there is no modsub image for i band
# VCC1450, there is no modsub image for g band
# VCC1466, curve_fit failed
# VCC1472, curve_fit failed
# VCC1486, prichet fitting failed.
# VCC1501, curve_fit failed
# VCC1512, create_table failed
# VCC1533, there is no modsub image for g band
# VCC1534, there is no modsub image for g band
# VCC1541, prichet fitting failed.
# VCC1580, curve_fit failed
# VCC1589, curve_fit failed
# VCC1591, curve_fit failed
# VCC1612, prichet fitting failed.
# VCC1623, curve_fit failed
# VCC1627, curve_fit failed
# VCC1668, there is no modsub image for g band
# VCC1669, there is no modsub image for g band
# VCC1675, create_table failed
# VCC1681, prichet fitting failed.
# VCC1704, create_table failed
# VCC1722, curve_fit failed
# VCC1750, curve_fit failed
# VCC1755, create_table failed
# VCC1780, prichet fitting failed.
# VCC1809, prichet fitting failed.
# VCC1820, curve_fit failed
# VCC1821, curve_fit failed
# VCC1829, curve_fit failed
# VCC1831, there is no modsub image for g band
# VCC1852, curve_fit failed
# VCC1853, curve_fit failed
# VCC1858, there is no modsub image for g band
# VCC1860, curve_fit failed
# VCC1873, curve_fit failed
# VCC1880, prichet fitting failed.
# VCC1926, curve_fit failed
# VCC1955, create_table failed
# VCC2010, prichet fitting failed.
# VCC2015, curve_fit failed
# VCC2029, curve_fit failed
# VCC2073, curve_fit failed
# VCC2083, there is no modsub image for g band
# VCC2086, curve_fit failed
GALAXY_TO_EXCLUDE = {'VCC0047', 'VCC0093', 'VCC0143', 'VCC0158', 'VCC0162', 'VCC0197', 'VCC0207', 'VCC0213', 'VCC0322',
                     'VCC0332', 'VCC0334', 'VCC0360', 'VCC0361', 'VCC0373', 'VCC0170', 'VCC0180', 'VCC0263', 'VCC0346', 'VCC0349',
                     'VCC0384', 'VCC0389', 'VCC0428', 'VCC0453', 'VCC0476', 'VCC0511', 'VCC0538', 'VCC0565', 'VCC0566',
                     'VCC0641', 'VCC0657', 'VCC0678', 'VCC0697', 'VCC0708', 'VCC0713', 'VCC0725', 'VCC0733', 'VCC0764',
                     'VCC0811', 'VCC0821', 'VCC0862', 'VCC0863', 'VCC0890', 'VCC0892', 'VCC0898', 'VCC0917', 'VCC0925',
                     'VCC0975', 'VCC0977', 'VCC0997', 'VCC1042', 'VCC1055', 'VCC1059', 'VCC1078', 'VCC1081', 'VCC1100',
                     'VCC1123', 'VCC1141', 'VCC1143', 'VCC1148', 'VCC1159', 'VCC1175', 'VCC1188', 'VCC1189', 'VCC1192',
                     'VCC1199', 'VCC1236', 'VCC1239', 'VCC1249', 'VCC1297', 'VCC1301', 'VCC1313', 'VCC1317', 'VCC1325',
                     'VCC1340', 'VCC1343', 'VCC1354', 'VCC1355', 'VCC1356', 'VCC1369', 'VCC1389', 'VCC1403', 'VCC1425',
                     'VCC1435', 'VCC1436', 'VCC1437', 'VCC1450', 'VCC1466', 'VCC1472', 'VCC1474', 'VCC1478', 'VCC1486',
                     'VCC1500', 'VCC1501', 'VCC1512', 'VCC1533', 'VCC1534', 'VCC1541', 'VCC1580', 'VCC1589', 'VCC1591',
                     'VCC1597', 'VCC1612', 'VCC1623', 'VCC1627', 'VCC1668', 'VCC1669', 'VCC1675',
                     'VCC1681', 'VCC1704', 'VCC1722', 'VCC1750', 'VCC1755', 'VCC1780', 'VCC1809', 'VCC1820', 'VCC1821',
                     'VCC1829', 'VCC1831', 'VCC1846', 'VCC1852', 'VCC1853', 'VCC1858', 'VCC1860', 'VCC1873', 'VCC1880',
                     'VCC1926', 'VCC1955', 'VCC2010', 'VCC2015', 'VCC2029', 'VCC2045', 'VCC2073', 'VCC2083', 'VCC2086'}

# Load the sample_gal.csv for basic information of galaxies such as magnitudes,
# half-light radii, positions, etc, which is used for the galaxy-light fitting.
# GALAXY_NAMES column 0
# GALAXY_CENTER_RA column 1
# GALAXY_CENTER_DEC column 2
# GALAXY_GMAG column 3
# GALAXY_EFFECTIVE_RADII column 4
# GALAXY_FIELD_NAMES column 5
# GALAXY_PSF_X column 6
# GALAXY_PSF_Y column 7
# extinction_magnitudes columns 8 to 12
def load_sample_gal():
    with open('sample_gal.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if row[0] in GALAXY_TO_EXCLUDE:
                continue
            GALAXY_NAMES.append(row[0])
            GALAXY_CENTER_RA.append(float(row[1]))
            GALAXY_CENTER_DEC.append(float(row[2]))
            GALAXY_GMAG.append(float(row[3]))
            GALAXY_EFFECTIVE_RADII.append(float(row[4]))
            GALAXY_FIELD_NAMES.append(row[5])
            GALAXY_PSF_X.append(int(row[6]))
            GALAXY_PSF_Y.append(int(row[7]))
            extinction_mags = []
            extinction_mags.append(float(row[8]))  # u
            extinction_mags.append(float(row[9]))  # g
            #extinction_mags.append(float(row[10])) # r, skip it.
            extinction_mags.append(float(row[11])) # i
            extinction_mags.append(float(row[12])) # z
            extinction_magnitudes[row[0]] = extinction_mags
