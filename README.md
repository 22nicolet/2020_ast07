# 2020_ast07
Set up directories as of the following:
.
├── images
│   └── VCC0940
├── pythonscripts
│   ├── add_stars.py
│   ├── completeness_test.py
│   ├── default.nnw
│   ├── gauss_2.5_5x5.conv
│   ├── ngvs.param
│   ├── ngvs.sex
│   ├── run_sextractor.py
└── results
    └── VCC0940

images  - Original images. Each galaxy is in a separate sub-directory such as VCC9040
results - result files from sextractor, or images with added stars. Each galaxy is in a separate sub-directory such as VCC9040
pythonscripts - python scripts, as well as the run directory

#commands to run to generate a completness function plot
python3 add_stars.py
python3 completeness_test.py
