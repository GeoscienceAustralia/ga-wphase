'''
Some utilities to make nice plots.
'''

import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.colors as mcols
from scipy.interpolate import griddata
from obspy.imaging.beachball import beach

try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    print("Warning: mpl_toolkits not installed.")
    pass

def topography_colormap(cpool_sea,cpool_land):
    '''
    Given two colours lists (2D numpy arrays with RGB code) for ocean and land
    respectevily, it constructs a color map (dictionary format) for the global
    topography. The division between land and sea should be evident in this new
    colormap if the colors have been choosen properly in the lists.
    '''

    cpool = np.vstack((cpool_sea,cpool_land ) )

    cpool /= 255.

    n_col_sea = cpool_sea.shape[0]
    n_col_land = cpool_land.shape[0]

    initval = (n_col_sea + n_col_land)*[3*[0]]
    color_map = {}
    color_map = OrderedDict([('red' ,  initval[:]),
                              ('green',  initval[:]),
                              ('blue'  ,  initval[:]) ])

    nodes_sea = np.linspace(0.,0.5,n_col_sea)
    nodes_land= np.linspace(0.5,1,n_col_land)
    nodes = np.concatenate((nodes_sea,nodes_land))

    for rgb, (color, val) in enumerate(color_map.iteritems()):
        for i,elem in enumerate(val):
            elem[0] = nodes[i]
            elem[1] = cpool[i,rgb]
            elem[2] = cpool[i,rgb]
            color_map[color][i] = elem[:]

    for color in color_map:
        color_map[color][n_col_sea-1][2] = color_map[color][n_col_sea][1]
        del color_map[color][n_col_sea]

    return color_map



def get_topography_cpool():
    cpool_sea = np.array( [[52 ,61  ,200   ],
                          [0  ,0   ,127   ],
                          [0  ,0   ,255   ],
                          [102, 205, 200  ],
                          [180, 255, 212  ]],dtype = 'float')

    cpool_land = np.array([ [0  , 159,  79  ],
                             [210, 184, 135  ],
                             [205, 133, 63   ],
                             [139, 69 , 19   ],
                             [180, 180, 180  ],
                             [240,240 ,240   ]],dtype = 'float')

    return(cpool_sea,cpool_land)



def get_boundaries(lats,lons,xpnd = 0.7):
    lnmn = lons.min()
    lnmx = lons.max()
    ltmn = lats.min()
    ltmx = lats.max()
    lllon = 0.5*(lnmn+lnmx)-xpnd*(lnmx-lnmn)
    lllat = 0.5*(ltmn+ltmx)-xpnd*(ltmx-ltmn)
    urlon = 0.5*(lnmn+lnmx)+xpnd*(lnmx-lnmn)
    urlat = 0.5*(ltmn+ltmx)+xpnd*(ltmx-ltmn)
    return lllon,lllat,urlon,urlat



def get_topography_slice(
    lllon,
    lllat,
    urlon,
    urlat,
    topofile = '/home/roberto/data/topography/ETOPO1_Ice_g_gdal.grd'):

    from netCDF4 import Dataset

    etopodata = Dataset(topofile)
    xran = etopodata.variables['x_range'][:]
    yran = etopodata.variables['y_range'][:]
    dim =  etopodata.variables['dimension'][:]
    nlats = dim[1]
    nlons = dim[0]
    dlons = (xran[1] - xran[0])/(nlons-1.)
    dlats = (yran[1] - yran[0])/(nlats-1.)
    lonminind = np.floor((lllon - xran[0])/dlons)
    latminind = np.floor((lllat - yran[0])/dlats)
    lonmaxind = np.ceil((urlon - xran[0])/dlons)
    latmaxind = np.ceil((urlat - yran[0])/dlats)
    lon_inter = lonmaxind  - lonminind
    lat_inter = latmaxind  - latminind
    lats = np.linspace(yran[0]+latminind*dlats,yran[0]+latmaxind*dlats,
                      lat_inter )
    lons = np.linspace(xran[0]+lonminind*dlons,xran[0]+lonmaxind*dlons,
                       lon_inter )
    elev = np.empty((int(lat_inter + 1), int(lon_inter)))
    for i,latind in enumerate(np.arange(nlats- latminind -2 , nlats - latmaxind - 3,
                               -1)):
        trim1 = int(nlons*latind + lonminind + 1 )
        trim2 = int(trim1 + lon_inter)
        elev[i,:] = etopodata.variables['z'][trim1 : trim2]

    return lats,lons,elev



def plot_field(
    latlons,
    field,
    plot_type='scatter',
    show_lats=True,
    show_lons=True,
    ax=None,
    clevs=None,
    topofile = '/home/roberto/data/topography/ETOPO1_Ice_g_gdal.grd',
    **kwargs):
    '''
    Plot the locations considered in the grid search.

    :param latlons: Tuples contain the lat/lon pairs.
    :type latlons: Two column numpy array.
    :param field: The field to plot. Must have the dimensions commensurate to *latlons*.
    :type field: Two dimensional numpy array.
    :param str plot_type: Either "scatter" (for points) or "contour" (for contours).
    :param bool show_lats: Should latitudes be shown on the y-axis.
    :param bool show_lons: Should longitudes be shown on the x-axis.
    :param ax: Axes for figure to add the plot to.
    :type ax: None or :py:class:`matplotlib.Axes`.
    :param list clevs: The levels for the contour plot. Ignored if *plot_type* is "scatter".
    :param str topfile: Path to file to load 'etopo' topography from.
    '''

    from mpl_toolkits.basemap import Basemap

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    cpool_sea,cpool_land = get_topography_cpool()
    cdict = topography_colormap(cpool_sea,cpool_land)
    relief = mcols.LinearSegmentedColormap('my_colormap', cdict, 256)
    # loading topographic data

    lats, lons =  latlons.T
    boundaries = get_boundaries(lats,lons,xpnd = 0.7)
    lllon,lllat,urlon,urlat = boundaries
    #Obtaning nice boundaries for the map
    m = Basemap(*boundaries,projection='merc',resolution='l',ax=ax)
    m.drawcoastlines(linewidth=1.5)
    m.drawcountries()
    if urlon - lllon < 5.:
        dlon = 1.0
    else:
        dlon = 2.0
    if urlat - lllat < 5.:
        dlat = 1.0
    else:
        dlat = 2.0
    if show_lats:
        lats_labs = [1,0,0,1]
    else:
        lats_labs = [0,0,0,0]

    if show_lons:
        lons_labs = [1,0,0,1]
    else:
        lons_labs = [0,0,0,0]
    parallels = np.arange(-90.,90,dlat)
    m.drawparallels(parallels,labels=lats_labs)
    meridians = np.arange(-180,180,dlon)
    m.drawmeridians(meridians,labels=lons_labs)
    #Ploting elevation
    if topofile:
        nxx = int((m.xmax-m.xmin)/500.)+1
        nyy = int((m.ymax-m.ymin)/500.)+1
        tlats,tlons,telev = get_topography_slice(lllon,lllat,urlon,urlat)
        topodat = m.transform_scalar(telev,tlons,tlats,nxx,nyy)
        im = m.imshow(topodat/1000.,relief,vmax=6.5,vmin=-6.5)
    else:
        m.fillcontinents(color='#FF9900',lake_color='#FF9900')
        m.drawmapboundary(fill_color='aqua')
    #Plotting field
    if plot_type == 'scatter':
        field_plot = m.scatter(lons, lats, latlon=True, **kwargs)
    elif plot_type == 'contour':
        grid_lats = np.sort(lats)
        grid_lons = np.sort(lons)
        grid_lats, grid_lons = np.meshgrid(grid_lats, grid_lons)
        grid_inter= griddata((lats,lons),field,(grid_lats,grid_lons),method = 'linear')
        #~ xpts, ypts = m(grid_lons,grid_lats)
        #~ m.contour(xpts,ypts,grid_slip, clevs,
                    #~ linewidths=1.5,colors='k',animated=True
        #~ print clevs
        if not np.any(clevs):
            clevs = np.linspace(field.min(),field.max(),8)
        field_plot = m.contour(grid_lons, grid_lats, grid_inter,clevs,latlon=True, **kwargs)
    output_dic = {'fig' : fig, 'ax': ax, 'field' : field_plot,
                  'limits' : boundaries, 'map' : m }
    if topofile:
        output_dic['topography'] = im
    return output_dic

def stacov(
        (elat, elon),
        lats,
        lons,
        trlist = None,
        mt = None,
        filename = None,
        fig = None,
        ax = None):
    """Plot a station coverage map."""

    if not fig:
        fig = plt.figure()

    if not ax:
        ax = fig.add_subplot()

    cov_map = Basemap(projection='ortho', lat_0=elat, lon_0=elon, resolution='c', area_thresh=10000., ax=ax)

    cov_map.fillcontinents(color='#FF9900', lake_color='#FF9900')
    cov_map.drawmapboundary(fill_color='aqua')
    cov_map.drawmeridians(np.arange(0,360,30), linewidth=0.4, color='grey', dashes=[5,5])
    cov_map.drawparallels(np.arange(-90,90,30), linewidth=0.4, color='grey', dashes=[5,5])

    x, y = cov_map(lons, lats)

    cov_map.scatter(x, y, s=60, c='blue', marker='o', edgecolors='none', zorder=10)

    x1, y1 = cov_map(elon,elat)
    size = 8e5
    if mt is not None:
        ax2 = plt.gca()
        b = beach(mt, xy=(x1, y1), width=size, linewidth=1, facecolor='r')
        b.set_zorder(10)
        ax2.add_collection(b)
    else:
        plt.plot(x1,y1,'y*',markersize=40)

    if trlist is not None:
        stnames = [trid.split(".")[1] for trid in trlist]
        for i in range(len(stnames)):
            plt.text(x[i], y[i], stnames[i], va="top", family="monospace", weight="bold")
    if filename:
        plt.savefig(filename, dpi=100, transparent=True, bbox_inches='tight')
