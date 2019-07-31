##FILE:
##PURPOSE:


##Below Section: IMPORT necessary functions
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import json
import os.path
import os.listdir
import re
from classUtils import func_timer
import astropy.constants as const

#Set helpful constants
dotest = True
if dotest:
    au0 = 1.49597870E13
    c0 = 2.99792458E10
    rsun0     = 6.9599E10
    msun0     = 1.9891E33
    lsun0     = 3.8525E33
    pc0     = 3.085678E18
    h0     = 6.6262000E-27
    kB0     = 1.3807E-16 #1.380658E-16 #From read_molecule_lambda.pro; #1.3807E-16 from natconst.pro
    mp0     = 1.6726231E-24
    G0     = 6.67259E-8
    tsun0     = 5.78E3
    mu0     = 2.3
    amu0    = 1.6605402E-24
else:
    amu0 = const.u.cgs.value #Atomic mass
    au0 = const.au.cgs.value
    c0 = const.c.cgs.value
    G0 = const.G.cgs.value
    mp0 = const.m_p.cgs.value #Mass of proton
    msun0 = const.M_sun.cgs.value
    rsun0 = const.R_sun.cgs.value
    h0 = const.h.cgs.value



##
class RadliteSpectrum():
    def __init__(infilename):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
    ##STORE AND FETCH METHODS
    def get_attr(self, attrname):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        pass
    #


    def change_attr(self, attrname, attrval):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        pass
    #


    def _set_attr(self, attrname, attrval):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        pass
    #



    ##GENERATIVE METHODS
    def gen_spec(self):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        #Print params, if verbose (dist, ssampling, obsres)
        #Read in all RADLite output from given list of runs
        self._read_radliteoutput()
        self._process_radliteoutput()
        #Check equal number of lines in line files and molfiles
        #Remove any duplicate lines
        #Calculate line continuum
        self._calc_continuum()
        #Separate line and continuum
        #Avoid extended line wings?
        #Some wavelength grid stuff for requested velocity sampling
        #Add continuum back
        #Convolve to requested resolving power
        #Resample to requested output grid
        #Noise, maybe
        #Make fits, maybe
        pass
    #



    ##OUTPUT DISPLAY METHODS
    def plot_spec(self):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        pass
    #



    ##READ METHODS
    def _read_core_radliteoutput(self, filelist):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: READ IN + PROCESS data from given RADLite subrun
        specfilename = filelist[0]
        molfilename = filelist[1]
        with open(specfilename, 'r') as openfile: #For spectrum
            specdata = openfile.readlines()
        with open(molfilename, 'r') as openfile: #For mol. data
            moldata = openfile.readlines()

        #FOR SPECTRUM FILE
        #Extract all spectral traits from spectrum file
        numlines = int(specdata[4])
        velinfo = specdata[6].split()
        vwidth = float(velinfo[0])
        vlsr = float(velinfo[1])
        incl = float(velinfo[2])

        #Iterate through molecular lines in subrun
        fluxes_list = [None]*numlines #To hold all fluxes
        vels_list = [None]*numlines #To hold all velocities
        iloc = 8 #Starting index within line files
        for ai in range(0, len(numlines)):
            numpoints = int(specdata[iloc+4]) #Num of freq.
            #Extract section of data for current molecular line
            sechere = specdata[iloc:(iloc+6+numpoints)])
            #Convert into 2D array of floats
            sechere = np.array([lhere.split()
                                for lhere in sechere]).astype(float)
            vels_list[ai] = sechere[:,0] #Velocities
            fluxes_list[ai] = sechere[:,1] #Fluxes
            #Increment place in overall file
            iloc = iloc + 6 + numpoints

        #FOR MOLECULAR DATA FILE
        #Extract molecular traits for this subrun
        molname = moldata[1] #Name of molecule
        molweight = float(moldata[3]) #Weight of molecule
        numlevels = int(moldata[5]) #Number of energy levels

        #Extract transitions
        iloc = 5 + 1 + numlevels + 3 #Starting index of transitions
        sechere = moldata[iloc:len(moldata)] #Section containing all trans.
        sechere = np.array([lhere.split()
                                    for lhere in sechere]) #Split out spaces
        gup_vals = sechere[:,1].astype(float) #Upper degeneracies
        glow_vals = sechere[:,2].astype(float) #Lower degeneracies
        A_vals = sechere[:,3].astype(float) #Einstein A coefficient
        Eup_vals = sechere[:,4].astype(float) #Upper energies
        freq_vals = sechere[:,5].astype(float) #Frequencies
        lvib_vals = sechere[:,6] #Vibrational levels
        lrot_vals = sechere[:,7] #Rotational levels


        ##Below Section: RETURN spectrum and molecular data + EXIT
        return {"flux":fluxes_list, "vel":vels_list, "vwidth":vwidth,
                    "molname":molname, "molweight":molweight,
                    "numlevels":numlevels, "numlines":numlines,
                    "lvib":lvib_vals, "lrot":lrot_vals,
                    "vlsr":vlsr, "incl":incl, "gup":gup_vals, "glow":glow_vals,
                    "A":A_vals, "Eup":Eup_vals, "freq":freq_vals}
    #


    def _read_radliteoutput(self):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: READ IN + PROCESS data from each given RADLite run
        runpath_list = self.get_attr("run_paths")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Starting process of reading RADLite output...")
            print("Extracting all Radlite output from "+str(runpath_list)+"...")

        #Iterate through given RADLite runs
        dict_list = [None]*len(runpath_list) #List to hold data from all subruns
        for ai in range(0, len(runpath_list)):
            pathhere = runpath_list[ai] #Path to current run
            #Assemble subruns (e.g., per core) in this run directory
            suballfiles = [namehere for namehere in os.listdir(pathhere)]
            regex_line = re.compile("linespectrum_moldata_*.dat") #Spec. file
            regex_mol = re.compile("moldata_*.dat") #Mol. data file
            #For all line spectrum files
            subspecfiles = [namehere for namehere in suballfiles
                                if re.match(regex_line, namehere)]
            #For all molecular data files
            submolfiles = [namehere for namehere in suballfiles
                                if re.match(regex_mol, namehere)]
            #Throw an error if not equal number of spectrum and mol. data files
            if len(subspecfiles) != len(submolfiles):
                raise ValueError("Whoa! "+pathhere+" does not contain an equal "
                        +"number of spectrum ('linespectrum_moldata_*.dat') "
                        +"and molecular data ('moldata_*.dat') files. Make "
                        +"you haven't tampered with any of the RADLite output "
                        +"files from running RadliteModel().run_radlite().")

            #Prepare pool of cores to read mol. lines and data for each subrun
            numsubs = len(subspecfiles)
            with mp.Pool(self.get_attr("numcores")) as ppool:
                subdicts = ppool.map(self._read_core_radlite,
                                        [[subspecfiles[bi], submolefiles[bi]]
                                            for bi in range(0, numsubs)])

            #Tack extracted subrun data to overall lists of data
            if self.get_attr("verbose"): #Verbal output, if so desired
                print(pathhere+" has "+str(len(subdicts))+" molecular lines.")
            dict_list[ai] = subdicts


        ##Below Section: GROUP all mol. line data; KEEP only unique lines
        #Flatten list of lists into a 1D list
        dict_list = [inlhere for outlhere in dict_list for inlhere in outlhere]
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done processing all RADLite runs!")
            print("There are a total of "+str(len(dict_list))+" molecular "
                    +"lines across all given RADLite runs.")
            print("Removing any duplicate line occurrences...")
        #Keep only unique lines
        uniqnames_list = [(dhere["molname"]+dhere["lvib"]+dhere["lrot"]
                            +str(dhere["gup"])+str(dhere["glow"])
                            +str(dhere["Eup"])) for dhere in dict_list]
        uniqinds = np.unique(uniqnames_list, return_index=True) #Uniq. line inds
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(len(uniqinds))+" molecular lines are being kept.")
        dict_list = dict_list[uniqinds] #Keep only unique molecular lines


        ##Below Section: STORE RADLite output + EXIT
        self._set_attr(attrname="_radliteoutput", attrval=dict_list)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done with process of reading RADLite output!")
        return
    #



    ##PROCESS METHODS
    def _process_radliteoutput(self):
        """
        DOCSTRING
        WARNING: This function is not intended for direct use by user.
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section:
        dict_list = self.get_attr("_radliteoutput")

    #


start at bottom; minimal commands; min. interpolation done

New proposal:
> Extrapolate each line to max box width
        # - line 'box' width: max( 3*leftspan of vel or 3*obsres ); done per line, then take max across all lines
        # - NOTE: defined new_vel to be 2*max_boxwidth*np.arange(0, num vels)/(num vels - 1) - max_boxwidth -> basically centers at 0 from -width to +width
> Extrapolate continuum for each line from min to max box width
> Subtract continuum from each line
> Set up x_all and y_all
        # - Determine min and max wavelengths from c/center_freq[where != 0]; then broaden by +/- max_vel*1E9*max_mu(or min_mu)/c
        # - Define wavelength grid with const. vel. res.: from min_mu to max_mu, spacing done as: x_all[i+1] = x_all[i] * (1 + res_el/2.9979...E5); threw an error if point count exceeded 1E8
> Interpolate line_only and cont. over x_all that fits them; add into y_all and some cont_all array of sorts
> Fill x_all and y_all
        # - x_mu = new_vel * 1E9 / center_freq + c/center_freq
        # - y_Jy = lines_int here; gsubs = where x_all >= min(x_mu) and x_all <= max(x_mu) - so extract segment of x_all that fits this given line x_mu
        # - y_all[gsubs] = y_all[gsubs] + <interpolate y_jy and x_mu over x_all[gsubs]>; so basically adding this line into overall final spectra, at overall final spectrum's resolution
> Add continuum back in for spec
> Convolve spec and line separately
> Convolve
> Resample
-
        # - line 'box' width: max( 3*leftspan of vel or 3*obsres ); done per line, then take max across all lines
        # - find highest vel-res (line with most vel. points)
        # - interpolate all other lines to have same vel-res (replace their old-vel)
        # - number of vels: taken to be 2*max_boxwidth/res_el + 1, where res_el = vel[1] - vel[0] for any vel_arr, since now same for each line
        # - defined an empty array of size (numlines, number of vels)
        # - defined new_vel to be 2*max_boxwidth*np.arange(0, num vels)/(num vels - 1) - max_boxwidth -> basically centers at 0 from -width to +width
        # - defined specx, specy, each of size numlines*num vels
        # - PER LINE
        # - interpolate cont. from leftmost to rightmost vel, flux of each line separately; interpolate over vel_arr
        # - subtract continuum from line info
        # - calculated line_int as interpolation of line, -maxvel to maxvel across new_vel
        # - calculated freqs??? as (1+vel_arr*1E9/c)*wavenum_central(?)
        # - recorded and stored line_flux = sorted freqs, line,/dist^2
        # - calculated line widths:
        #   - 2.35482*sqrt(TOTAL(line[gsubs]*vel_arr[gsubs]^2)/TOTAL(line[gsubs]) where gsubs = where(abs(line) > max(abs(line)*0.2)
        # - saved l2cs[i] as max(line)/mean(cont)
        # - DONE WITH PER LINE
        # - Determine min and max wavelengths from c/center_freq[where != 0]; then broaden by +/- max_vel*1E9*max_mu(or min_mu)/c
        # - Define wavelength grid with const. vel. res.: from min_mu to max_mu, spacing done as: x_all[i+1] = x_all[i] * (1 + res_el/2.9979...E5); threw an error if point count exceeded 1E8
        # - Define wavelength grid on requested output velocity sampling: x_out, from min_mu to max_mu, with x_out[i+1] = x_out[i] * (1 + sampling / 2.9979...E5); same error if too many points
        # - FOR EACH LINE
        # - x_mu = new_vel * 1E9 / center_freq + c/center_freq
        # - y_Jy = lines_int here; gsubs = where x_all >= min(x_mu) and x_all <= max(x_mu) - so extract segment of x_all that fits this given line x_mu
        # - y_all[gsubs] = y_all[gsubs] + <interpolate y_jy and x_mu over x_all[gsubs]>; so basically adding this line into overall final spectra, at overall final spectrum's resolution
        # - DONE WITH PER LINE
        # - CONVOLUTION
        # - Calculate ngauss as ceil(3*obsres/res_el)
        # - Calculate obsres_sampling as obsres/res_el
        # - Calculate gaus as exp(-(arange(0, ngauss)-(ngauss-1)/2.0)^2/obsres_sampling^2 * 2/log(2)
        # - Sort then interpolate cont. (value at vel=0 for each line), c/center_freqs, across x_all
        # - Set l_only = y_all * 1E23/dist^2
        # - Added cont. as y_all = (y_all+c_all)*1E23/dist^2
        # - Convolved as convol(y_all, gauss, total(gauss),/edge_truncate)
        # - Convolved l_only in same way, using l_only instead of y_all
        # - RESAMPLING
        # - Interpolated y_out, l_only, c_all against x_all to be across x_out; multiplied 1E23/dist^2 into continuum since not done yet
        # - Possibly able to add noise?
        # - Record line, flux, spectrum
        # - Write to fits file, if so desired (but also separate function)




        # - Interpolate (linear) cont. from leftmost - rightmost vel, flux of each line
        # -
        #####Below Section: SORT + EXTRACT all spec. fluxes, cont., and mol. info
        np.argsort(all
        allcont_arr = np.concatenate(



        #Throw error if inclination is not the same across all data



        allspec_arrs = [] #To hold each extracted molecular line
        allEu_vals = [] #To hold extracted upper energies
        allElow_vals = [] #To hold extracted lower energies
        #allnumline_vals = [] #To hold extracted numbers of molecular lines
        allnumfreq_vals = [] #To hold extracted number of frequency points
        alldeltav_vals = [] #To hold extracted delta-velocities
        allvlsr_vals = [] #To hold extracted vlsr values
        allincl_vals = [] #To hold extracted source inclination values

            #Iterate through









        ##Below Section: Extract flux, velocity, frequency, and energy info




    #



    ##CALCULATION METHODS
    #



    ##WRITE METHODS
    def write_fits(self):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        pass
    #
#













#
