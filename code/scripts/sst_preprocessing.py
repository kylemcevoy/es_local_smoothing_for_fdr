#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 12:07:18 2024

@author: Kyle McEvoy

data pre-processing of SST
sst.mnmean.nc file contains the NOAA Extended Reconstructed SST V5
obtained from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
downloaded 05/15/23 covering dates 1854-01-01 to 2023-04-01
see paper for full citation.

this file loads the SST dataset, removes the monthly mean to get anomalies,
subsets to observations with latitudes in [-60, 60] and month = january,
and writes out to a .csv file.
"""

import xarray as xr

sst = xr.open_mfdataset('../../data/sst.mnmean.nc')
sst = sst.drop_vars('time_bnds')
sst = sst.astype('float64')

sst.load()

sst_anom = sst.groupby('time.month') - sst.groupby('time.month').mean('time')

sst_jan_mid = (sst_anom.
               sel(lat = slice(60, -60)).
               sel(time = sst['time.month'] == 1))

(sst_jan_mid.to_dataframe()
 .dropna()
 .reset_index()
 .to_csv("../../data/sst_midlat_anom_jan.csv"))
 
 