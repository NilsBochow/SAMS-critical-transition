import numpy as np 
import pandas as pd 
from matplotlib import pyplot as plt
import sys
from scipy import stats
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import PathPatch

from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
import os 
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Since Basemap is deprecated, you probably have to set the PROJ_LIB variable manually before loading basemap, replace $HOME with the directory where you have the proj folder
os.environ['PROJ_LIB'] = r'$HOME/anaconda3/pkgs/proj4-5.2.0-he6710b0_1/share/proj/'
from mpl_toolkits.basemap import Basemap, maskoceans

"""
Script for creating a map like Fig. 6. 

The data in the folder files_map are taken from ERA5.
Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2023): ERA5 monthly averaged data on single levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS), DOI: 10.24381/cds.f17050d7 
"""
directory = dir_path = os.path.dirname(os.path.realpath(__file__)) + "/files_map/"
lat_lon_bigger = np.load(directory + "lat_lon_SA_bigger.npz")
lat_bigger = lat_lon_bigger['lat']
lon_bigger = lat_lon_bigger['lon']
la_bigger = len(lat_bigger)
lo_bigger = len(lon_bigger)


lat_lon = np.load(directory + "lat_lon_SA.npz")
lat = lat_lon['lat']
lon = lat_lon['lon']

tlen = 492

lon_min= -80
lon_max = -40
lat_min = -20
lat_max = 10
lat_hatch = lat[np.where(lat == lat_max)[0][0]  : np.where(lat == lat_min)[0][0]+1]
lon_hatch = lon[np.where(lon == lon_min)[0][0] : np.where(lon == lon_max)[0][0] + 1]



def draw_screen_poly( lats, lons, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( list(xy), facecolor='none', edgecolor = "white", lw =3 )
    plt.gca().add_patch(poly)

def calc_kendall_tau(array1, array2):
    tau, p_value = stats.kendalltau(array1, array2)
    return tau, p_value

def kendall_tau_array(array1): 
    """
    Calculate the Kendall Tau rank correlation coefficient for each element along the first dimension of a 3D array.
    Assumes a window length of 20 years (240 months) of the input array.
    """

    array1 = array1[240::, :, :]
    tau = np.zeros((array1.shape[1], array1.shape[2]))
    p_value = np.zeros((array1.shape[1], array1.shape[2]))
    for i in range(array1.shape[1]):
        for j in range(array1.shape[2]):
            tau[i, j], p_value[i, j] = calc_kendall_tau(np.arange(array1.shape[0]),array1[:,i,j])

    return tau

        
def plot_map(gradient_popt, array_string, levels_max, hatched):
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['font.weight'] = 'bold'
    
    
    size = 14
    trend = gradient_popt[:, :]
   
    fig     = plt.figure(figsize=(7, 6.0))
    ax      = fig.add_subplot(111)


    m = Basemap(epsg = '3395', llcrnrlat=-20,urcrnrlat=10,llcrnrlon=-80,urcrnrlon=-40)
    draw_screen_poly( [-12.5,-12.5,-4,-4], [-72.5, -62.5,-62.5, -72.5], m )
    lons, lats = np.meshgrid(lon, lat)
    lons_bigger, lats_bigger = np.meshgrid(lon_bigger, lat_bigger)
    x,y = m(lons, lats)
    xx, yy = m(lons_bigger, lats_bigger)

    v_wind = (np.load(directory + "v_wind_djf_mean.npy"))
    u_wind = (np.load(directory + "u_wind_djf_mean.npy"))


    v_wind = maskoceans(lons_bigger, lats_bigger, v_wind)
    u_wind = maskoceans(lons_bigger, lats_bigger, u_wind)
    m.quiver(xx[::10, ::10], yy[::10, ::10], u_wind[::10, ::10], v_wind[::10, ::10], scale = 120, zorder=2, lw = 5, color = "dimgrey", alpha = 1)


    trend_masked = maskoceans(lons, lats, trend)
    condition = (trend_masked>0)
    lons_hatched, lats_hatched = np.meshgrid(lon_hatch, lat_hatch)
    x_hatch,y_hatch = m(lons_hatched, lats_hatched)
    if hatched == True: 
        if array_string == "autocorrelation": 
            hatch_area = np.load(directory + "ar1_p_value_SA_whole.npy")
            hatch_area = np.ma.masked_greater(hatch_area, 0.05)
            hatch_area = maskoceans(lons_hatched, lats_hatched, hatch_area)
        else:
            hatch_area = np.load(directory + "var_p_value_SA_whole.npy")
            hatch_area = np.ma.masked_greater(hatch_area, 0.05)
            hatch_area = maskoceans(lons_hatched, lats_hatched, hatch_area)
            


    contour = m.contourf(x, y, trend_masked ,levels = np.linspace(-levels_max,levels_max,21), cmap = 'coolwarm', extend = 'both', latlon= False)

    m.drawcoastlines(linewidth = 3)
    m.drawmapscale(-46,-17, -35, -0, 500, linewidth=3, fontsize=12)
    
    for c in contour.collections:
        c.set_edgecolor("face")


    
    
    if hatched == True:
        m.pcolor(x_hatch, y_hatch, hatch_area, hatch='..', alpha=0.)
    divider = make_axes_locatable(ax)

    cax = fig.add_axes([0.2, 0.05, 0.6, 0.03])

    cbar = fig.colorbar(contour, cax =cax, orientation = "horizontal", ticks=ticker.LinearLocator(5))

    if array_string == "autocorrelation":
        cbar.set_label(r'Kendall $\tau$ '+ array_string + ' of P', size =size)
    else: 
        cbar.set_label(r'Kendall $\tau$ '+ array_string + ' of P', size =size)
    cbar.ax.tick_params(labelsize=size, size=6, width=3)


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.savefig(directory + array_string + '_map_precip.pdf', bbox_inches='tight')

kendall_ar1 = np.load(directory + "kendall_ar1.npy")
plot_map(kendall_ar1, "autocorrelation", 1, hatched =True )


kendall_variance = np.load(directory + "kendall_variance.npy")
plot_map(kendall_variance, "variance", 1, hatched =True )