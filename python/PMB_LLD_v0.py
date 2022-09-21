# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 10:19:51 2022
    This is the call function of PMB-LLD Method for users. It applies to single location with its
    SST timeseries (one-dimension vector).
    
    Input part is for users to input the start_yr, end_yr of the whole timeseries, as well 
    as the moving baseline period.
    
    Theta and climatology detrending part is to detrend linear trend of the 31-yr fixed SST
    for theta and climatology calculation.
    
    SST detrending part is to detrend linear trend of moving 31-yr SST for MHW detecting.
    
    MHW detecting part firstly calculates the theta and climatology of the fixed 31-yr and
    then detect MHWs using detrended SST.
    
    Ver.0
    update date: Jun 7, 2022
@author: Cindy
"""

import numpy as np
from datetime import date
import scipy.signal as scis
import cal_theta_clim as ctc
import detect_mhws as detm

# import matplotlib.pyplot as plt

# --------------------------------------- Input -----------------------------------------
#  Please input the start year, end year of the whole timeseries ...
#  (including both theta calculating and MHW detetcing) and the moving baseline period

start_yr  = 1901
end_yr    = 2100
mov_base  = 31 # moving baseline period
mov_base2 = int((mov_base-1)/2) # half of the moving baseline period


# ----------------------------- Theta and climatology detrending -----------------------------
# This part is for a fixed 31-year baseline period for the detrending of seasonally varying ...
# theta and climatology

t_theta = np.arange(date(start_yr,1,1).toordinal(),date(start_yr+31-1,12,31).toordinal()+1)
sst_theta = np.zeros(len(t_theta))
sst_theta[0] = 0 # Initial condition
a = 0.85 # autoregressive parameter
for i in range(1,len(t_theta)):
    sst_theta[i] = a*sst_theta[i-1] + 0.75*np.random.randn() + 0.5*np.cos(t_theta[i]*2*np.pi/365.25)
sst_theta = sst_theta - sst_theta.min() + 5.
sst_theta = scis.detrend(sst_theta)


# ------------------------------------- SST detrending --------------------------------------
# This part is for a moving 31-year baseline period for the detrending of SST used for ...
# MHW detetcing
t_all = np.arange(date(start_yr+mov_base,1,1).toordinal(),date(end_yr,12,31).toordinal()+1)
sst = np.zeros(len(t_all))
sst[0] = 0 # Initial condition
a = 0.85 # autoregressive parameter
for i in range(1,len(t_all)):
    sst[i] = a*sst[i-1] + 0.75*np.random.randn() + 0.5*np.cos(t_all[i]*2*np.pi/365.25)
sst = sst - sst.min() + 5.

for i in range(start_yr+mov_base+mov_base2, end_yr-mov_base2+1):
    if i == start_yr+mov_base+mov_base2:
        sst_all = sst[date(i-mov_base2, 1, 1).toordinal()-date(start_yr+mov_base,1,1).toordinal()+1:
                      date(i+mov_base2,12,31).toordinal()-date(start_yr+mov_base,1,1).toordinal()+2]
        sst_all = scis.detrend(sst_all)
        sst_now = sst_all[date(i, 1, 1).toordinal()-date(i-mov_base2,1,1).toordinal()+1:
                                         date(i,12,31).toordinal()-date(i-mov_base2,1,1).toordinal()+2]
    else:
        sst_all = sst[date(i-mov_base2, 1, 1).toordinal()-date(start_yr+mov_base,1,1).toordinal()+1:
                      date(i+mov_base2,12,31).toordinal()-date(start_yr+mov_base,1,1).toordinal()+2]
        sst_all = scis.detrend(sst_all)   
        sst_now = np.append(sst_now, sst_all[date(i, 1, 1).toordinal()-date(i-mov_base2,1,1).toordinal()+1:
                                     date(i,12,31).toordinal()-date(i-mov_base2,1,1).toordinal()+2])


# --------------------------------------- MHW detecting --------------------------------------- 
# Calculate theta and climatology using year 1901-1931 (31yr)
[thres,seas] = ctc.cal_theta(t_theta, sst_theta)

# Detect MHWs using year 1947-2085

# ========== First repmat thres and seas to 1947-2085 timeseries (the same as sst_all) ==============
t_now = np.arange(date(start_yr+mov_base+mov_base2,1,1).toordinal(),date(end_yr-mov_base2,12,31).toordinal()+1)
T = len(t_now)
year = np.zeros((T))
month = np.zeros((T))
day = np.zeros((T))
doy = np.zeros((T))

for i in range(T):
    year[i] = date.fromordinal(t_all[i]).year
    month[i] = date.fromordinal(t_all[i]).month
    day[i] = date.fromordinal(t_all[i]).day

# Leap-year baseline for defining day-of-year values
year_leapYear = 2012 # This year was a leap-year and therefore doy in range of 1 to 366
t_leapYear = np.arange(date(year_leapYear, 1, 1).toordinal(),date(year_leapYear, 12, 31).toordinal()+1)

month_leapYear = np.zeros((len(t_leapYear)))
day_leapYear = np.zeros((len(t_leapYear)))
doy_leapYear = np.zeros((len(t_leapYear)))
for tt in range(len(t_leapYear)):
    month_leapYear[tt] = date.fromordinal(t_leapYear[tt]).month
    day_leapYear[tt] = date.fromordinal(t_leapYear[tt]).day
    doy_leapYear[tt] = t_leapYear[tt] - date(date.fromordinal(t_leapYear[tt]).year,1,1).toordinal() + 1
# Calculate day-of-year values
# The length of doy is 1947-2085
for tt in range(T):
    doy[tt] = doy_leapYear[(month_leapYear == month[tt]) * (day_leapYear == day[tt])]

clim = {}
clim['thresh'] = np.NaN*np.zeros(T)
clim['seas'] = np.NaN*np.zeros(T)
clim['thresh'] = thres[doy.astype(int)-1]
clim['seas'] = seas[doy.astype(int)-1]
# ==================================================================================================

# Then detecting
mhw = detm.detect_mhw(t_now, sst_now, clim)
# ------------------------------------------- Ending --------------------------------------

# -------------- Testing -----------
# plt.plot(thres)
# plt.plot(seas)
# plt.show()

