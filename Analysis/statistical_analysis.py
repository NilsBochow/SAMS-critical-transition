import numpy as np 
import pandas as pd 
from scipy import stats
from matplotlib import pyplot as plt
from statsmodels.tsa.seasonal import STL
from scipy.optimize import curve_fit

"""
Script for the statistical analysis and significance tests. 
"""


def calc_kendall_tau(array):
    """
    Simple function to calculate the Kendall tau of an ordered array.
    """
    tau, p_value = stats.kendalltau(np.arange(array.size), array)
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
    
    for i in range(ns):

        popt1 = calc_kendall_tau(nts[i])
        stat[i] = popt1

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