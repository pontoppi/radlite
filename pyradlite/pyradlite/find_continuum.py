import astropy.io.fits as pf
import numpy as np
import scipy.signal as signal
import scipy.interpolate as ip
import os

def find_continuum(wavelength,flux,window_length=41,nan=True,iter_max=2,median_filter=True,
                   linetype='emission'):
    """
    Method that calculates a continuum either for an emission- or an absorption-line spectrum. It requires
    scipy >= v0.14.0.
    
    Parameters
    ----------
    wavelength : one-dimensional numpy float array.
    flux : one-dimensional numpy float array.
    window_length : integer scalar, must be odd.
        The size of the smoothing window. Features similar to or wider than this is considered continuum.
    nan : boolean
        If true, any NaNs/infinite values in the input flux array are removed. 
    iter_max : integer scalar
        Number of iterations.
    median_filter : boolean
        If true, a median filter is applied to make the algorithm a little outlier-resistant. 

    Returns
    -------
    continuum : one-dimensional float array.
        The continuum on the same wavelength sampling as the input. 
    
    Example
    -------
    import matplotlib.pylab as plt
    from find_continuum import *
    
    path = '/Users/pontoppi/WORK/OBSPROGS/H2O_IRS/REDUCTIONS/'
    file = 'RNO90_SH_final.fits'

    fullpath = os.path.join(path,file)

    spectrum = pf.getdata(fullpath) 
    wavelength = spectrum['wave'].flatten()
    flux = spectrum['spec'].flatten()

    continuum = find_continuum(wavelength,flux,nan=True)
    
    plt.plot(wavelength,flux)
    plt.plot(wavelength,continuum)
    plt.show()
    
    """
    
    if nan:
        gsubs = np.where(np.isfinite(flux))
        flux = flux[gsubs]
        wave = wavelength[gsubs]
    
    if median_filter:
        flux = signal.medfilt(flux,3)
        
    restricted_flux = np.concatenate((flux[::-1],flux,flux[::-1]))
    wave = np.concatenate((wave-np.max(wave)+np.min(wave),wave,wave-np.min(wave)+np.max(wave)))
    original_flux = np.copy(restricted_flux)
    
    for i in range(iter_max):    
        smoothed_flux = signal.savgol_filter(restricted_flux,window_length=window_length,
                                             polyorder=2,mode='nearest')
                                              
        if linetype is 'absorption':
            gsubs = np.where(original_flux>smoothed_flux)
        elif linetype is 'emission':
            gsubs = np.where(original_flux<smoothed_flux)
        else:
            raise Exception(ValueError)
        
        restricted_flux = original_flux[gsubs]
        restricted_wave = wave[gsubs]
        
        interpolate_flux = ip.interp1d(restricted_wave,restricted_flux,bounds_error=False)
        restricted_flux = interpolate_flux(wave)
            
    continuum = signal.savgol_filter(restricted_flux,window_length=window_length,polyorder=2,
                                               mode='nearest')    

    interpolate_continuum = ip.interp1d(wave,continuum,bounds_error=False)

    return interpolate_continuum(wavelength)
    
    