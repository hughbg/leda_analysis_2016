#!/usr/bin/env python2.7

"""
# 00_download_data.py

Downloads required from zenodo repository:
    Danny C. Price et. al. (2017).
    Large-aperture Experiment to Detect the Dark Ages -- February 2016 [Data set].
    Zenodo. http://doi.org/10.5281/zenodo.400690

"""

import os

ZENODO_URL = 'https://zenodo.org/record/400690/files'
data_files = ["outriggers_2016-01-27_14H10M25S.h5"]

print("Downloading calibration data...")
os.system("wget %s/cal_data.zip" % ZENODO_URL)
print("unzipping...")
os.system("unzip -o cal_data.zip")
os.system("rm cal_data.zip")

if not os.path.exists("data"):
    print("Creating data directory...")
    os.mkdir("data")

for filename in data_files:
    print("Downloading %s..." % filename)
    os.system("wget %s/%s" %(ZENODO_URL, filename))
    os.system("mv %s data/" % filename)