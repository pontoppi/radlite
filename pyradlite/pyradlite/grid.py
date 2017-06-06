import sys
import glob

import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as ip

import astropy.io.fits as pf
import find_continuum as fc

def read_grid(path):
    parameters = pf.getdata(path+'/run_table.fits')
    nsegments = len(parameters)
    keys = parameters.names
    
    mdata = []
    count = 0
    increment = nsegments/2
    point = nsegments/10
    
    sys.stdout.write("Reading the RADLite grid\n")
    for parameter in parameters:    
        resolutions = glob.glob(path+'/'+parameter[0]+'/model_*.fits')
        #There could be models without a resolution tag
        resolutions = resolutions + glob.glob(path+'/'+parameter[0]+'/model.fits')
        words = parameter[0].split('_')
        ID = words[1]

        for resolution in resolutions:
            words = resolution.split('_')
            try:
                width = float(words[-1][0:-5])
            except:
                width = None
                
            model = pf.getdata(resolution)
        
            mwave = model['wave'].flatten()
            mflux = model['spec'].flatten()
            mcontinuum = fc.find_continuum(mwave,mflux,nan=True)
            mlines = mflux-mcontinuum
            mdict = {'path':resolution,'wave':mwave,'flux':mflux,'continuum':mcontinuum,'lines':mlines,
                     'resolution':width, 'ID':ID}
            for key in keys:
                mdict[key.lower()] = parameter[key]
            mdata.append(mdict)
        percent = 100*(count/(nsegments-1))
        sys.stdout.write("Progress: %d%%   \r" % percent)
        sys.stdout.flush()
        count += 1
    return mdata

def chi2(dwave, dflux, dstddev, mwave, mflux,wbound=None):
    
    if wbound is None:
        mmin = np.min(mwave)
        mmax = np.max(mwave)
        wbound = (mmin,mmax)
    
    gsubs = np.where((dwave>wbound[0]) & (dwave<wbound[1]))
        
    model = ip.interp1d(mwave,mflux)
    chi2 = np.sum((dflux[gsubs]-model(dwave[gsubs]))**2./dstddev[gsubs]**2.)
    ndegrees = len(gsubs[0])
    return chi2, ndegrees


    
