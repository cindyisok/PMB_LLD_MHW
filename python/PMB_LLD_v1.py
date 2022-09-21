# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 20:09:25 2022

    This is the call function of PMB-LLD Method for users. It applies to both gridded data
    and single location with the SST timeseries (one-dimension vector).
    
    Input part is for users to input the start_yr, end_yr of the whole timeseries, as well 
    as the moving baseline period.
    
    Theta and climatology detrending part is to detrend linear trend of the 31-yr fixed SST
    for theta and climatology calculation.
    
    SST detrending part is to detrend linear trend of moving 31-yr SST for MHW detecting.
    
    MHW detecting part firstly calculates the theta and climatology of the fixed 31-yr and
    then detect MHWs using detrended SST.

    Ver.1
    
@author: Cindy
"""

import numpy as np
from datetime import date
import scipy.signal as scis
import scipy.io as scio
import cal_theta_clim as ctc
from detect_mhws import blockAverage
import detect_mhws as detm
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

#
# --------------------------------------- Input -----------------------------------------
#  Please input the start year, end year of the whole timeseries ...
#  (including both theta calculating and MHW detetcing) and the moving baseline period
#
start_yr  = 1901
end_yr    = 2100
mov_base  = 31 # moving baseline period
mov_base2 = int((mov_base-1)/2) # half of the moving baseline period

#
# ----------------------------- Theta and climatology detrending -----------------------------
# This part is for a fixed 31-year baseline period for the detrending of seasonally varying ...
# theta and climatology
#
t_theta = np.arange(date(start_yr,1,1).toordinal(),date(start_yr+31-1,12,31).toordinal()+1)
data = scio.loadmat('/home/sundi/CESM/SST_NA_LR/SST_year_%04d.mat'%start_yr)
lon = data['lon']
lat = data['lat']
sst = data['SST0']


# read other years
for iy in range(start_yr+1,start_yr+30+1):
    filename = '/home/sundi/CESM/SST_NA_LR/SST_year_%04d.mat'%iy
    data = scio.loadmat(filename)
    sst0 = data['SST0']
    if (np.mod(iy,4)==0 and np.mod(iy,100)!=0) or np.mod(iy,400)==0:
        sst_leap_day = np.nanmean(sst0[:,:,58:60],2)
        sst_leap_day = sst_leap_day.reshape([sst_leap_day.shape[0],sst_leap_day.shape[1],1])
        sst0 = np.concatenate((sst0[:,:,0:59],sst_leap_day,sst0[:,:,59:]),axis=2)
    sst = np.concatenate((sst,sst0),2)

nx = sst.shape[0]
ny = sst.shape[1]

thres_all = np.full([nx,ny,366],np.nan)
seas_all = np.full([nx,ny,366],np.nan)
for ix in range(nx):
    for iy in range(ny):
        sst_theta = sst[ix,iy,:]
        if ~np.isnan(sst_theta[0]):
            sst_theta = scis.detrend(sst_theta)
            
            # ===========================================================
            # Calculate theta and climatology using year 1901-1931 (31yr)
            [thres,seas] = ctc.cal_theta(t_theta, sst_theta)
            # ===========================================================
            
            thres_all[ix,iy,:] = thres
            seas_all[ix,iy,:] = seas

# np.save('m90.npy',thres_all)
# np.save('mclim.npy',seas_all)

#
# ------------------------------------- SST detrending --------------------------------------
# This part is for a moving 31-year baseline period for the detrending of SST used for ...
# MHW detetcing
#
t_all = np.arange(date(start_yr+mov_base,1,1).toordinal(),date(end_yr,12,31).toordinal()+1)

data = scio.loadmat('G:/SST_NA_LR/SST_year_%04d.mat'%(start_yr+mov_base))
lon = data['lon']
lat = data['lat']
sst = data['SST0']
nx = sst.shape[0]
ny = sst.shape[1]
if (np.mod(start_yr+mov_base,4)==0 and np.mod(start_yr+mov_base,100)!=0) or np.mod(start_yr+mov_base,400)==0:
    sst_leap_day = np.nanmean(sst[:,:,58:60],2)
    sst_leap_day = sst_leap_day.reshape([sst_leap_day.shape[0],sst_leap_day.shape[1],1])
    sst = np.concatenate((sst[:,:,0:59],sst_leap_day,sst[:,:,59:]),axis=2)

for iy in range(start_yr+mov_base+1,end_yr+1):
    filename = 'G:/SST_NA_LR/SST_year_%04d.mat'%iy
    data = scio.loadmat(filename)
    sst0 = data['SST0']
    if (np.mod(iy,4)==0 and np.mod(iy,100)!=0) or np.mod(iy,400)==0:
        sst_leap_day = np.nanmean(sst0[:,:,58:60],2)
        sst_leap_day = sst_leap_day.reshape([sst_leap_day.shape[0],sst_leap_day.shape[1],1])
        sst0 = np.concatenate((sst0[:,:,0:59],sst_leap_day,sst0[:,:,59:]),axis=2)
    sst = np.concatenate((sst,sst0),2)
    
sst_now = np.full([nx,ny,date(end_yr-mov_base2,12,31).toordinal()-date(start_yr+mov_base+mov_base2,1,1).toordinal()+1],np.nan)   
 
for i in range(start_yr+mov_base+mov_base2, end_yr-mov_base2+1):
    sst_all = sst[:,:,date(i-mov_base2, 1, 1).toordinal()-date(start_yr+mov_base,1,1).toordinal():
                      date(i+mov_base2,12,31).toordinal()-date(start_yr+mov_base,1,1).toordinal()+1]
    for ix in range(nx):
        for iy in range(ny):
            sst_day = sst_all[ix,iy,:]
            if ~np.isnan(sst_day[0]):
                sst_day = scis.detrend(sst_day)
                sst_now[ix,iy,date(i, 1, 1).toordinal()-date(start_yr+mov_base2+mov_base,1,1).toordinal():
                        date(i,12,31).toordinal()-date(start_yr+mov_base2+mov_base,1,1).toordinal()+1] = sst_day[date(i, 1, 1).toordinal()-date(i-mov_base2,1,1).toordinal():
                                                                                                          date(i, 12, 31).toordinal()-date(i-mov_base2,1,1).toordinal()+1]
    # print(i) # for test
     
     
# np.save('sst_now.npy',sst_now)
# thres_all = np.load('m90.npy')
# seas_all = np.load('mclim.npy')

#
# --------------------------------------- MHW detecting --------------------------------------- 
#
# ========== First repmat thres and seas to 1947-2085 timeseries (the same as sst_all) ==============
#
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

#
# Then MHW detecting
#    
clim = {}
clim['thresh'] = np.NaN*np.zeros(T)
clim['seas'] = np.NaN*np.zeros(T)

mean_cum = np.full([nx,ny],np.nan)
trend_cum = np.full([nx,ny],np.nan)
dtrend_cum = np.full([nx,ny],np.nan)

for ix in range(nx):
    for iy in range(ny): 
        clim['thresh'] = thres_all[ix,iy,doy.astype(int)-1]
        clim['seas'] = seas_all[ix,iy,doy.astype(int)-1]
        
        # ==================================================
        # Detect MHWs
        mhw = detm.detect_mhw(t_now, sst_now[ix,iy,:], clim)
        # ==================================================
        
        # mean and trend calculation
        mhwBlock = blockAverage(t_now, mhw, temp=sst_now[ix,iy,:])
        mean, trend, dtrend = detm.meanTrend(mhwBlock)
        
        # for plotting
        mean_cum[ix,iy] = mean['intensity_cumulative']
        trend_cum[ix,iy] = trend['intensity_cumulative']
        dtrend_cum[ix,iy] = dtrend['intensity_cumulative']
        
    # print(ix) # for test

# ------------------------------- ending -----------------------------------------

#
# ---------------------------------- plot ---------------------------------------
#
    
p_cum = np.abs(trend_cum) - dtrend_cum

fig = plt.figure(figsize=(10, 10))
map = Basemap(
    projection='merc',
    llcrnrlon=120,
    llcrnrlat=20,
    urcrnrlon=240,
    urcrnrlat=70,
    lat_0=20.,lon_0=120.
    # resolution='i'
    )
# Lat, Lon = map(*np.meshgrid(lon0, lat0))
trend_cum[p_cum<0] = np.nan # whether significant in 95% confidence level

map.drawcoastlines(linewidth=0.5)
map.drawmapboundary(fill_color='w')
parallels = np.arange(20, 70, 16)
meridians = np.arange(120, 240.1, 30)
map.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=20)
map.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=20)
map.fillcontinents(color='gray')
x,y = map(lon,lat)
cs = map.pcolormesh(x, y, trend_cum*100, vmin = -150, vmax = 150, cmap='jet')
cb = map.colorbar(cs, extend='both')
cb.set_label('℃·day/century',fontsize=15)
plt.title('k$_{I}$ from 1947-2085 : PMB-LLD' ,fontsize=20)
plt.show()

# ==================================================================================
#

                                                                                                                     