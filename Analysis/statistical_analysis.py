import numpy as np 
import pandas as pd 
from scipy import stats
from matplotlib import pyplot as plt
from statsmodels.tsa.seasonal import STL
from scipy.optimize import curve_fit
import netCDF4 as nc4 
import os 

"""
Script for the statistical analysis and significance tests. 
"""
directory = dir_path = os.path.dirname(os.path.realpath(__file__)) + "/statistics_example/"

def calc_kendall_tau(array):
    """
    Simple function to calculate the Kendall tau of an ordered array.
    """
    tau, p_value = stats.kendalltau(np.arange(array.size), array)
    return tau, p_value


def kendall_tau_array(array1): 
    """
    Simple function to apply calc_kendall_tau() over whole array.
    """
    array1 = array1[240::, :, :]
    tau = np.zeros((array1.shape[1], array1.shape[2]))
    p_value = np.zeros((array1.shape[1], array1.shape[2]))
    for i in range(array1.shape[1]):
        for j in range(array1.shape[2]):
            tau[i, j], p_value[i, j] = calc_kendall_tau(array1[:,i,j])
    return tau

def fourrier_surrogates(ts, ns):
    """
    Generating Fourier surrogates for significance testing. See methods for details. 
    Takes the orignal time series ts and the number of surrogates ns as input. Returns the surrogate time series new_ts. 
    """
    ts_fourier  = np.fft.rfft(ts)
    random_phases = np.exp(np.random.uniform(0, 2 * np.pi, (ns, ts.shape[0] // 2 + 1)) * 1.0j)
    ts_fourier_new = ts_fourier * random_phases
    new_ts = np.real(np.fft.irfft(ts_fourier_new))
    return new_ts

def kendall_tau_test(ts, ns, trend, mode = 'fourier'):
    """
    Calculates the Kendall tau of each surrogate time series and compares it to the original kenddal tau. 
    Takes the original time series ts as input, the number of surrogates ns and the Kendall tau of the original time series. 
    Returns the p-value.
    """
    tlen = ts.shape[0]
    if mode == 'fourier':
        tsf = ts - ts.mean()
        nts = fourrier_surrogates(tsf, ns)
    stat = np.zeros(ns)
    tlen = nts.shape[1]
    p_value = np.zeros(ns) #not needed
    for i in range(ns):

        popt1 = calc_kendall_tau(nts[i])
        stat[i], p_value[i] = popt1

    p = 1 - stats.percentileofscore(stat, trend) / 100.
    return p   


def EWS_ar1(array):
    """
    Function to calculate the temporal autocorrelation at lag-1 for every element in a spatial array. The input array is expected to have the shape (time, lon/lat, lat/lon). 
    Can be used for 1D time series as well.
    """
    window = 20*12 #rolling window length in months
    ar1 = np.zeros_like(array) 
    for i in range(array.shape[1]):
        for j in range(array.shape[2]):
            if np.ma.is_masked(array[:, i, j]):
                continue #if masked continue, e.g., over ocean 
            else: 
                tr = int(10 * 12 + 1) #Trend smoother length in months
                array_stl = STL(array[:, i, j], period = int(12), seasonal=int(13), trend = tr, robust = 'true').fit() #STL decomposition
                array_trend = array_stl.trend
                array_season = array_stl.seasonal
                array_ds = array_stl.resid 
                ar1[:, i, j] = pd.Series(array_ds).rolling(window).apply(lambda x: x.autocorr(), raw=False).values #calculate the autocorrelation at lag-1 for every site of the detrended and deseasoned array
    return ar1



def EWS_variance(array):
    """
    Function to calculate the temporal variance for every element in a spatial array. The input array is expected to have the shape (time, lon/lat, lat/lon).
    Can be used for 1D time series as well.
    """
    window = 20*12 #rolling window length in months
    variance = np.zeros_like(array) 
    for i in range(array.shape[1]):
        for j in range(array.shape[2]):
            if np.ma.is_masked(array[:, i, j]):
                continue #if masked continue, e.g., over ocean 
            else: 
                tr = int(10 * 12 + 1) #Trend smoother length in months
                array_stl = STL(array[:, i, j], period = int(12), seasonal=int(13), trend = tr, robust = 'true').fit() #STL decomposition
                array_trend = array_stl.trend
                array_season = array_stl.seasonal
                array_ds = array_stl.resid 
                variance[:, i, j] = pd.Series(array_ds).rolling(window=window).var().values #calculate the variance for every site of the detrended and deseasoned array
    return variance

def fit_trend_soil(array):
    """
    Function to fit a linear function to a spatial array at every site. The input array is expected to have the shape (time, lon/lat, lat/lon). 
    """
    def lin(x, p0, p1):
        return p1 + p0 * x
    popt = np.zeros((2, array.shape[1], array.shape[2])) 
    for i in range(array.shape[1]):
        for j in range(array.shape[2]): 
            popt1, pcov1 = curve_fit(lin, np.arange(40), array[:, i, j], absolute_sigma = True, p0 =[-2, 800], maxfev = 1000) 
            popt[:, i, j] = popt1
    return popt


### Example calculation of Variance/Autocorrelation, Trend and significance 
gpcc_precip = nc4.Dataset(directory + "GPCC/precip.monitor.mon.total.1x1.v2020_SouthAmerica.nc")["precip"][0:471,:,:] 
#Source:  GPCC operated by DWD under the auspices of the World Meteorological Organization (WMO)

gpcp_precip = nc4.Dataset(directory + "GPCP/subset_SouthAmerica.nc")["precip"][:,:,:] 
#Source: Adler, R.F., G.J. Huffman, A. Chang, R. Ferraro, P. Xie, J. Janowiak, B. Rudolf, U. Schneider, S. Curtis, D. Bolvin, A. Gruber, J. Susskind, and P. Arkin, 2003: The Version 2 Global Precipitation Climatology Project (GPCP) Monthly Precipitation Analysis (1979-Present). J. Hydrometeor., 4,1147-1167.

tlen = gpcp_precip.shape[0]
variance = EWS_variance(gpcp_precip)
ar1 = EWS_ar1(gpcp_precip)
np.savez_compressed(directory + "GPCP/variance_gpcp_SouthAmerica.npz", data=variance.data, mask=variance.mask)
np.savez_compressed(directory + "GPCP/ar1_gpcp_SouthAmerica.npz", data=ar1.data, mask=ar1.mask)
popt_ar1 = kendall_tau_array(ar1)
popt_variance = kendall_tau_array(variance)



np.save(directory + "GPCP/popt_variance.npy", popt_variance)
np.save(directory + "GPCP/popt_ar1.npy", popt_ar1)


def calculate_p_value(array, popt, array_string): 
    array_p_value=np.zeros((array.shape[1], array.shape[2]))
    for i in range(array.shape[1]):
        for j in range(array.shape[2]):
            if np.isnan(array[240::, i, j]).any():
                pass
            elif (not np.any(array[240::, i, j])):
                pass
            else: 
                print(i,j)
                array_p_value[i, j] = kendall_tau_test(array[240::, i, j], 10000, popt[i, j])
    np.save(directory + "GPCP/"+ array_string + "_p_value.npy", array_p_value)

calculate_p_value(ar1, popt_ar1, "ar")