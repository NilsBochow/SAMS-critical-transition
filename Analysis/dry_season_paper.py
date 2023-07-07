import numpy as np 
import pandas as pd 
import sys
from scipy import stats
from skimage.morphology import binary_opening
from scipy.ndimage import uniform_filter1d
import netCDF4 as nc4

### anomaly, pentad sum, 30d climatological 
""" 
For this script, there are some files needed: 
precip_all_years_lat_lon : A file with spatial daily precipitation for all years. For the paper, we use the freely available ERA5 precipitation.
climatological_mean: A file with the climatological mean of the precipitation or can be calculated from the precip_all_years_lat_lon variable
"""



"""
Some short functions used for calculation. A variable precip_nc4 is needed to infer the latitude/longitude and time variable.  
"""

def load_lat_lon(): 
    lat = precip_nc4.variables['lat'][:]
    lon = precip_nc4.variables['lon'][:]
    la = len(lat)
    lo = len(lon)
    n = la * lo
    return lo, la, n, lat, lon

def load_tlen(): 
    time =  precip_nc4.variables["time"][:]
    return time.shape[0]

def lin_reg(DSL_yearly): 
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(np.arange(41), DSL_yearly) 
    return slope, intercept

def rolling_mean_3D(window, array, axis): 
    hW = window//2
    L = array.shape[axis]-window+1   
    indexer = [slice(None) for _ in range(array.ndim)]
    indexer[axis] = slice(hW,hW+L)

    return uniform_filter1d(array, window, axis=axis)[tuple(indexer)]

def rolling_sum_3D(a, n) :
    ret = np.cumsum(a, axis=0, dtype=float)
    ret[ n:] = ret[ n:] - ret[ :-n]
    return ret[ n - 1:]

def rainfall_anomaly(precip): 
    anomaly = precip - np.mean(np.reshape(climatological_mean, (la*lo))[area_indices])
    return anomaly

def filter_wet_spells(WSL, n): 
    return np.uint8(binary_opening(WSL, np.ones((n), np.uint8)))



def calculate_DSL_anomaly(year): 
    """
    Calculates the rainfall anomaly and the beginning and 
    retreat of the wet season based on the method of Liebmann et al. 2012. 
    """
    precip_one_year = precip_all_years_lat_lon[year]
    precip_mean = mean_area_precip(precip_one_year, lat_min=lat_min, lat_max=lat_max, lon_min=lon_min, lon_max=lon_max)
    pads_number = 370-precip_mean.shape[0]
    
    rolling_precip = np.nanmean(np.pad(precip_mean, pads_number, mode = "constant", constant_values = np.nan)[pads_number:].reshape(-1,5), axis = 1)
    anomaly = rainfall_anomaly(rolling_precip)
    rainfall_anomaly_cum = np.cumsum(anomaly, axis = 0)
    onset = np.argmax(rainfall_anomaly_cum, axis=0)
    mask = np.less.outer(np.arange(len(rainfall_anomaly_cum)), onset)  
    rainfall_anomaly_cum[mask] = np.nan  
    retreat = np.nanargmin(rainfall_anomaly_cum, axis=0) 
    DSL = retreat - onset

    return DSL, onset, retreat 

def longestZeroSeqLength(a):
    """
    Calculates the longest sequence of zeros in a given array.
    """
    chg = np.abs(np.diff(np.equal(a, 0).view(np.int8), prepend=[0], append=[0]))
    rng = np.where(chg == 1)[0]
    if rng.size == 0: return 0   
    rng = rng.reshape(-1, 2)
    return np.subtract(rng[:,1], rng[:,0]).max(), rng[np.subtract(rng[:,1], rng[:,0]).argmax(),:][0], rng[np.subtract(rng[:,1], rng[:,0]).argmax(),:][1]

def determine_DSL(WSL):
    condition = (WSL==0)
    DSL = np.count_nonzero(condition, axis = 0)

    return DSL

tlen = load_tlen()
lo, la, n, lat, lon = load_lat_lon()

def coordinate_indices_from_ra(lat, lon, lat_max, lat_min, lon_max, lon_min):
    """ 
    Function to find indices of the chosen area.
    """
    la = len(lat)
    lo = len(lon)
    n = la * lo
    indices = np.arange(n).reshape((la, lo))
    lat_max = lat[np.argmin(np.abs(lat - lat_max))]
    lat_min = lat[np.argmin(np.abs(lat - lat_min))]
    lon_max = lon[np.argmin(np.abs(lon - lon_max))]
    lon_min = lon[np.argmin(np.abs(lon - lon_min))]
    ra_indices = indices[np.where(lat == lat_min)[0][0]  : np.where(lat == lat_max)[0][0] + 1, np.where(lon == lon_min)[0][0] : np.where(lon == lon_max)[0][0] + 1 ]
    return np.array(np.unique(ra_indices), dtype = 'int')

def define_tdjf(lat, lon, lat_max, lat_min, lon_max, lon_min): 
    tdjf = coordinate_indices_from_ra(lat, lon, lat_max, lat_min, lon_max, lon_min)
    return tdjf

"""
South Amazonia area.
"""
lon_min= -70
lon_max = -50 
lat_min = -15 
lat_max = -5
area_indices = define_tdjf(lat, lon, lat_max, lat_min, lon_max, lon_min)

def mean_area_precip(array, lat_min, lat_max, lon_min, lon_max): 
    array = np.reshape(array, (-1,la*lo))
    array_mean = np.mean(array[:, area_indices], axis= 1)
    return array_mean




def calculate_DSL_pentad(year): 
    """
    Calculates DSL with adoption of the method of Marengo et al. 2001.
    """
    pentads=8
    precip_one_year = precip_all_years_lat_lon[year]
    precip_mean = mean_area_precip(precip_one_year, lat_min=lat_min, lat_max=lat_max, lon_min=lon_min, lon_max=lon_max)
    ## padding array because of leap years ##
    pads_number = 370-precip_mean.shape[0]
    rolling_precip = np.nanmean(np.pad(precip_mean, pads_number, mode = "constant", constant_values = np.nan)[pads_number:].reshape(-1,5), axis = 1)
    
    WSL =np.zeros((rolling_precip.shape[0]))
    
    # Values above the climatoligical mean correspond to 1 and we calculate the rolling sum over 8 pentads for further calculations. 
    # This corresponds to the number of consecutive "wet" pentads.
    WSL[rolling_precip>np.mean(np.reshape(climatological_mean, (la*lo))[area_indices])] = 1
    WSL_rolling_sum = rolling_sum_3D(WSL, n=pentads)

    
    # Part where onset and retreat date are calculated. 
    retreat_date = determine_wet_season_retreat(WSL_rolling_sum)
    onset_date = determine_WSL_onset(WSL_rolling_sum, retreat_date)

    # DSL = onset - retreat date of wet season 
    DSL_length = onset_date - retreat_date 

    return DSL_length, retreat_date, onset_date


def determine_wet_season_retreat(WSL_rolling_sum): 
    """ 
    Calculation of yearly wet-season retreat date in pentads. 
    Initial assumption is to find the date with 6 preceeding "wet" pentads and 2 subsequent "dry" pentads. For further information see Methods or Marengo et al. 2001.
    """

    preceeding_wet_days = (WSL_rolling_sum >= 6) 
    subsequent_dry_days = (WSL_rolling_sum <= 2) 

    # This part determines the onset date #
    retreat_date = np.argmax(np.logical_and((preceeding_wet_days), np.roll(subsequent_dry_days, -8, axis = 0)), axis = 0) +8 

    condition = (np.less(retreat_date, 9))
    i = 1
    while np.any(condition): 
        """
        Relaxiation of condition if 6-2 fails. 
        """
        preceeding_wet_days = (WSL_rolling_sum >= 6-i) 
        subsequent_dry_days = (WSL_rolling_sum <= 2+i) 
        retreat_date = np.argmax(np.logical_and((preceeding_wet_days), np.roll(subsequent_dry_days, -8, axis = 0)), axis = 0) +8 
        condition = (np.less(retreat_date, 9))
        i +=1 
        if i > 2: 
            break
    return retreat_date


def determine_WSL_onset(WSL_rolling_sum, retreat_date): 
    """ 
    Calculation of yearly wet-season onset date in pentads. 
    Initial assumption is to find the date with 2 preceeding "dry" pentads and 6 subsequent "wet" pentads. For further information see Methods or Marengo et al. 2001.
    """

    preceeding_dry_days = (WSL_rolling_sum <= 2) 
    subsequent_wet_days = (WSL_rolling_sum >= 6) 
    
    # This part determines the onset date #
    logical_composition = np.logical_and((preceeding_dry_days), np.roll(subsequent_wet_days, -8, axis = 0))
    mask = np.less.outer(np.arange(len(logical_composition)), retreat_date)  
    logical_composition[mask] = False  
    onset_date = np.nanargmax(logical_composition, axis = 0)  + 8 #+8 since we take rolling sum over 8 pentads


    i = 1
    condition = (np.less(onset_date,retreat_date))
    while np.any(condition): 
        """
        Relaxiation of condition if 6-2 fails. 
        """
        preceeding_dry_days = (WSL_rolling_sum <= 2+i) 
        subsequent_wet_days = (WSL_rolling_sum >= 6-i)
        logical_composition = np.logical_and((preceeding_dry_days), np.roll(subsequent_wet_days, -8, axis = 0))
        mask = np.less.outer(np.arange(len(logical_composition)), retreat_date)  
        logical_composition[mask] = False
        onset_date = np.nanargmax(logical_composition, axis = 0)  + 8
        i +=1
        condition = (np.less(onset_date,retreat_date))
        if i > 2: 
            break
    
    
    return onset_date

def calculate_DSL_30d(year): 
    """ 
    Calculation of DSL via Method 1 in Methods.
    """

    window = 30
    precip_one_year = precip_all_years_lat_lon[year]
    precip_mean = mean_area_precip(precip_one_year, lat_min=lat_min, lat_max=lat_max, lon_min=lon_min, lon_max=lon_max)
    rolling_precip = rolling_mean_3D(window, precip_mean, 0)

    WSL =np.zeros((rolling_precip.shape[0]))
    WSL[rolling_precip>np.mean(np.reshape(climatological_mean, (la*lo))[area_indices])] = 1

    WSL = filter_wet_spells(WSL, 10)

    return WSL

def calculate_all_methods(): 
    """ 
    Store results of all methods in one array.
    """
    DSL_lat_lon = np.zeros((3, 41))
    DSL_lat_lon_onset = np.zeros((3, 41))
    WSL_lat_lon_onset = np.zeros((3, 41))

    slope = np.zeros((3, 2))
    slope_onset_DSL = np.zeros((3, 2))
    slope_onset_WSL = np.zeros((3, 2))

    for year in range(41):
        WSL = calculate_DSL_30d(year)
        DSL_lat_lon[0, year], DSL_lat_lon_onset[0, year], WSL_lat_lon_onset[0, year] = np.apply_along_axis(longestZeroSeqLength, 0, WSL)
        
        DSL_lat_lon[1, year], DSL_lat_lon_onset[1, year], WSL_lat_lon_onset[1, year] = calculate_DSL_anomaly(year)
        
        DSL_lat_lon[2, year], DSL_lat_lon_onset[2, year], WSL_lat_lon_onset[2, year] = calculate_DSL_pentad(year)

    for i in range(3): 
        slope[i, :] = lin_reg(DSL_lat_lon[i, :]) 
        slope_onset_DSL[i, :] = lin_reg(DSL_lat_lon_onset[i, :])    
        slope_onset_WSL[i, :] = lin_reg(WSL_lat_lon_onset[i, :])    

    return WSL, DSL_lat_lon, DSL_lat_lon_onset, WSL_lat_lon_onset, slope, slope_onset_DSL, slope_onset_WSL

WSL, DSL, DSL_onset, WSL_onset, slope, slope_onset_DSL, slope_onset_WSL =  calculate_all_methods()