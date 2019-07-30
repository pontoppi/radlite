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
        self._read_radlite()
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
    def _read_core_radlite(self, filelist):
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
        gup_vals = sechere[:,1] #Upper degeneracies
        glow_vals = sechere[:,2] #Lower degeneracies
        A_vals = sechere[:,3] #Einstein A coefficient
        Eup_vals = sechere[:,4] #Upper energies
        wavenum_vals = sechere[:,5] #Wavenumbers


        ##Below Section: RETURN spectrum and molecular data + EXIT
        return {"flux":fluxes_list, "vel":vels_list, "vwidth":vwidth,
                    "molname":molname, "molweight":molweight,
                    "numlevels":numlevels, "numlines":numlines,
                    "vlsr":vlsr, "incl":incl, "gup":gup_vals, "glow":glow_vals,
                    "A":A_vals, "Eup":Eup_vals, "wavenum":wavenum_vals}
    #


    def _read_radlite(self):
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
            with mp.Pool(numsubs) as ppool:
                subdicts = ppool.map(self._read_core_radlite,
                                        [[subspecfiles[bi], submolefiles[bi]]
                                            for bi in range(0, numsubs)])

            #Tack extracted subrun data to overall lists of data
            dict_list[ai] = subdicts


        ##Below Section: CONCATENATE all mol. line data; KEEP only unique lines
        #Flatten list of lists into a 1D list
        dict_list = [inlhere for outlhere in dict_list for inlhere in outlhere]

!!!
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




        pass
        filename = 'linespectrum_moldata_'+str(thread)+'.dat'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        dum = file.readline()
        dum = file.readline()
        moldata_base = file.readline()
        moldata_file = file.readline()
        nlines = int(file.readline())
        nfreq = int(file.readline())
        deltav,vzero,incl = file.readline().split()
        lines = []
        for i in np.arange(nlines):
            file.readline()
            upper,lower = file.readline().strip().split()
            upper = int(upper)
            lower = int(lower)
            freq = float(file.readline())
            freq = float(freq)
            vzero = float(file.readline())
            nfreq = int(file.readline())

            vels = np.zeros(nfreq)
            fluxes = np.zeros(nfreq)

            file.readline()
            for j in np.arange(nfreq):
                vel,flux = file.readline().strip().split()
                vels[j] = float(vel)
                fluxes[j] = float(flux)

            lines.append({'velocity':vels,'flux':fluxes,'freq':freq,'upper':upper,'lower':lower})

        return lines
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
