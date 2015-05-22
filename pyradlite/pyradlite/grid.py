import sys
import glob

import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as ip

import astropy.io.fits as pf
import readcol as rc
import find_continuum as fc

def read_grid(path):
    parameters = pf.getdata(path+'/run_table.fits')
    nsegments = len(parameters)
    
    mdata = []
    count = 0
    increment = nsegments/20
    point = nsegments/100
    
    sys.stdout.write("Reading the RADLite grid\n")
    for parameter in parameters:    
        resolutions = glob.glob(path+'/'+parameter[0]+'/model_*.fits')
        words = parameter[0].split('_')
        ID = words[1]

        for resolution in resolutions:
            words = resolution.split('_')
            width = float(words[-1][0:-5])
            
            model = pf.getdata(resolution)
        
            mwave = model['wave'].flatten()
            mflux = model['spec'].flatten()
            mcontinuum = fc.find_continuum(mwave,mflux,nan=True)
            mlines = mflux-mcontinuum
            mdata.append({'path':resolution,'wave':mwave,'flux':mflux,'continuum':mcontinuum,'lines':mlines,
                          'parameter':np.array(parameter),'resolution':width, 'ID':ID})

        if (count % (5*point) == 0):
            sys.stdout.write("\r[" + "=" * (count / increment) +  " " * ((nsegments - count)/ increment) + "]" +  str(count / point) + "%")
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


    
