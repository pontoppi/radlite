import os
import sys
import glob
import warnings
from collections import OrderedDict

import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as ip

import astropy.io.fits as pf
from . import find_continuum as fc

def grid_to_mapping(grid):
    
    fields = grid['grid_pars'].dtype.names[1:]
    pars_asdict = OrderedDict()
    pars_asdict['wave'] = grid['spectra'][0]['wave']
    shape = pars_asdict['wave'].shape
    for field in fields:
        pars_asdict[field.lower()] = np.unique(grid['grid_pars'][field])
        shape += (pars_asdict[field.lower()].shape[0],)
    
    map_values = np.zeros((pars_asdict['wave'].shape[0],len(grid['spectra'])))
    for ii,model in enumerate(grid['spectra']):
        map_values[:,ii] = model['flux']
        
    map_values = map_values.reshape(shape)
    import pdb;pdb.set_trace()
    
    return mapping

def read_grid(path, gridroot=None):
    '''
    Method to read a radlite grid, including spectra and metadata. 
    
    Returns: Dictionary
        Containing a list of individual radlite model dictionaries and the input grid parameters.
    
    '''
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
            
            try:  
                model   = pf.getdata(resolution,1)
                moldata = pf.getdata(resolution,2)
        
                mwave = model['wave'].flatten()
                mflux = model['spec'].flatten()
                mcontinuum = fc.find_continuum(mwave,mflux,nan=True)
                mlines = mflux-mcontinuum
                mdict = {'path':resolution,'wave':mwave,'flux':mflux,'continuum':mcontinuum,'lines':mlines,
                         'resolution':width, 'ID':ID, 'moldata':moldata}
                for key in keys:
                    mdict[key.lower()] = parameter[key]
                mdata.append(mdict)
            except:
                # If the model does not exist, ignore it. 
                warnings.warn('Failed to read model '+ resolution)
                
            
            
        percent = 100*(count/(nsegments-1))
        sys.stdout.write("Progress: %d%%   \r" % percent)
        sys.stdout.flush()
        count += 1
    return {'spectra':mdata,'grid_pars':parameters}

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


    
