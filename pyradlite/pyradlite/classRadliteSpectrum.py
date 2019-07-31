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
try:
    import astropy.io.fits as fitter
except ImportError:
    import pyfits as fitter

#Set helpful constants
dotest = True
if dotest:
    au0 = 1.49597870E13
    c0 = 2.99792458E10
    cinmu0 = c0*1E4
    cinkm0 = c0/1.0E5
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
    cinmu0 = c0*1.0E4 #mu/s
    cinkm0 = c0/1.0E5 #km/s
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
        ##Below Section: REVIEW user parameters for spectra
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Starting the gen_spec method!")
            print("Reviewing user-defined parameters for spectra...")
            print("Chosen dist [pc]: "+str(self.get_attr("dist")))
            print("Chosen obsres [km/s]: "+str(self.get_attr("obsres")))
            print("Chosen vsampling [km/s]: "+str(self.get_attr("vsampling")))


        ##Below Section: READ IN + PROCESS all RADLite output from given runs
        self._read_radliteoutput()
        self._process_radliteoutput()


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done with gen_spec!\n")
        return
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
        numlines = int(specdata[4]) #Number of mol. lines
        #maxnumpoints = int(specdata[5]) #Max number of data points per line
        incl = float(specdata[6].split()[2]) #Velocity info
        #vwidth = float(velinfo[0])
        #vlsr = float(velinfo[1])
        #incl = float(velinfo[2])

        #Iterate through molecular lines in subrun
        linedict_list = [{}]*numlines #To hold all fluxes
        iloc = 8 #Starting index within line files
        for ai in range(0, len(numlines)):
            freqcen = float(specdata[iloc+2]) #Central frequency
            numpoints = int(specdata[iloc+4]) #Num of freq.
            #Extract section of data for current molecular line
            sechere = specdata[iloc:(iloc+6+numpoints)])
            #Convert into 2D array of floats
            sechere = np.array([lhere.split()
                                for lhere in sechere]).astype(float)
            #Record information for this mol. line
            linedict_list[ai]["vel"] = sechere[:,0] #Velocities
            linedict_list[ai]["flux"] = sechere[:,1] #Fluxes
            linedict_list[ai]["freqcen"] = freqcen #Central frequency
            #Increment place in overall file
            iloc = iloc + 6 + numpoints

        #FOR MOLECULAR DATA FILE
        moldict = {} #Dictionary to hold molecular data
        #Extract molecular traits for this subrun
        moldict["molname"] = moldata[1] #Name of molecule
        moldict["molweight"] = float(moldata[3]) #Weight of molecule
        moldict["numlevels"] = int(moldata[5]) #Number of energy levels

        #Extract transitions
        iloc = 5 + 1 + numlevels + 3 #Starting index of transitions
        sechere = moldata[iloc:len(moldata)] #Section containing all trans.
        sechere = np.array([lhere.split()
                                    for lhere in sechere]) #Split out spaces
        moldict["gup"] = sechere[:,1].astype(float) #Upper degeneracies
        moldict["glow"] = sechere[:,2].astype(float) #Lower degeneracies
        moldict["A"] = sechere[:,3].astype(float) #Einstein A coefficient
        moldict["Eup"] = sechere[:,4].astype(float) #Upper energies
        moldict["freq"] = sechere[:,5].astype(float) #Frequencies
        moldict["lvib"] = sechere[:,6] #Vibrational levels
        moldict["lrot"] = sechere[:,7] #Rotational levels


        ##Below Section: RETURN spectrum and molecular data + EXIT
        #return {"flux":fluxes_list, "vel":vels_list, "vwidth":vwidth,
        #            "maxnumpoints":maxnumpoints,
        #            "molname":molname, "molweight":molweight,
        #            "numlevels":numlevels, "numlines":numlines,
        #            "lvib":lvib_vals, "lrot":lrot_vals,
        #            "vlsr":vlsr, "incl":incl, "gup":gup_vals, "glow":glow_vals,
        #            "A":A_vals, "Eup":Eup_vals, "freq":freq_vals}
        return {"line":linedict_list, "mol":moldict}
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
        linedict_list = [None]*len(runpath_list) #List to hold subrun line data
        moldict_list = [None]*len(runpath_list) #List to hold subrun mol. data
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
            linedict_list[ai] = [shere["line"] for shere in subdicts] #Line data
            moldict_list[ai] = [shere["mol"] for shere in subdicts] #Mol. data

        #Flatten list of lists into a 1D list
        linedict_list = [inlhere for outlhere in linedict_list
                        for inlhere in outlhere]
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done processing all RADLite runs!")
            print("There are a total of "+str(len(linedict_list))+" molecular "
                    +"lines across all given RADLite runs.")


        ##Below Section: KEEP only unique lines
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Removing any duplicate line occurrences...")
        #Keep only unique lines
        uniqnames_list = [(dhere["molname"]+dhere["lvib"]+dhere["lrot"]
                            +str(dhere["gup"])+str(dhere["glow"])
                            +str(dhere["Eup"])) for dhere in linedict_list]
        uniqinds = np.unique(uniqnames_list, return_index=True) #Uniq. line inds
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(len(uniqinds))+" molecular lines are being kept.")
        linedict_list = linedict_list[uniqinds] #Keep only unique mol. lines


        ##Below Section: CHECK for any incompatible values
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Checking for signs of incompatible sources...")
        #Check that all sources have same inclination
        if len(np.unique([dhere["incl"] for dhere in moldict_list])) > 1:
            raise ValueError("Whoa! The RADLite runs you specified don't all "
                            +"have the same source inclination!")


        ##Below Section: STORE RADLite output + EXIT
        findicts = {"line":linedict_list, "mol":moldict_list} #Line and mol data
        self._set_attr(attrname="_radliteoutput", attrval=findicts)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done with process of reading RADLite output!\n")
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
        ##Below Section: ACCESS read-in RADLite output
        linedict_list = self.get_attr("_radliteoutput")["line"]
        moldict_list = self.get_attr("_radliteoutput")["mol"]
        continterpkind = self.get_attr("cont_interp")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Starting to process read-in RADLite output...")


        ##Below Section: DETERMINE largest encompassing mol. line 'box'
        #NOTE: 'box' refers to x-axis span (e.g., velocity) that will be used...
        #...to encompass each molecular line
        #NOTE: 'vel' = velocity
        #Extract max span of all mol. lines (each line is from -span to +span)
        maxvspan = np.max([np.abs(dhere["vel"][0]) for dhere in linedict_list])
        #Choose box width that is 3*maxvspan OR 3*(specified obs. resolution)
        boxwidth = np.max([(3*maxvspan), (3*self.get_attr("obsres"))])
        #Also extract finest velocity resolution used for mol. lines
        indmax = np.argmax([dhere["maxnumpoints"]
                            for dhere in linedict_list]) #Index of maximum
        vres = linedict_list[indmax]["vel"][1] - linedict_list[indmax]["vel"][0]


        ##Below Section: SUBTRACT out continuum for each mol. line
        emonly_list = [None]*len(linedict_list) #Init. array for em. only
        contcen_list = [None]*len(linedict_list) #Init. array for center cont.
        #Iterate through mol. lines
        for ai in range(0, len(linedict_list)):
            dhere = linedict_list[ai] #Current mol. line
            #Linearly interpolate continuum from line edges; calculate at vel=0
            tempys = [dhere["flux"][0], dhere["flux"][-1]] #Fluxes to interp.
            tempxs = [dhere["vel"][0], dhere["vel"][-1]] #Vels. to interp.
            interpfunc = interper(x=tempxs, y=tempys, kind=continterpkind)
            #Record central continuum and isolated emission for this line
            emonly_list[ai] = dhere["flux"] - interpfunc(dhere["vel"]) #Em. only
            contcen_list[ai] = interpfunc(0) #Continuum at central vel. (vel=0)


        ##Below Section: PREPARE arrays to hold full spec. and output spec.
        #Determine min. and max. wavelengths (wavelen.) [mu] for full spectrum
        freqcen_list = np.array([dhere["freqcen"]
                                    for dhere in linedict_list]) #Center freqs.
        freqcen_list = freqcen_list[freqcen_list != 0] #Cut out any freq=0
        #For min. mu
        muminraw = np.min(cinmu0/1.0/freqcen_list) #Unbroadened min. mu value
        mumin = muminraw - (boxwidth*1.0E9*muminraw/cinmu0) #Broad. by box width
        #For max. mu
        mumaxraw = np.max(cinmu0/1.0/freqcen_list) #Unbroadened max. mu value
        mumax = mumaxraw + (boxwidth*1.0E9*mumaxraw/cinmu0) #Broad. by box width

        #Prepare wavelen. [mu] array for full spectrum with uniform vel. spacing
        growthfrac = (1 + (vres/1.0/cinkm0)) #Fraction by which mu will grow
        maxlen = np.floor(np.log(mumax/1.0/mumin)
                            / np.log(growthfrac)) + 1 #mu array size
        fullmu_arr = np.array([(mumin*(growthfrac**ai))
                                for ai in range(0, maxlen)])
        #NOTE: Used compound interest equation to calculate maxlen, fullmu_arr
        fullem_arr = np.zeros(maxlen) #Initialized; will hold full em-spec.

        #Prepare wavelen. [mu] array for output spectrum at desired resolution
        userres = self.get_attr("vsampling") #User-specified sampling res.
        userfrac = (1 + (userres/1.0/cinkm0)) #Fraction by which mu will grow
        userlen = np.floor(np.log(mumax/1.0/mumin)
                            / np.log(userfrac)) + 1 #mu array size
        outmu_arr = np.array([(mumin*(userfrac**ai))
                                for ai in range(0, userlen)])


        ##Below Section: INTERPOLATE + ADD each mol. emission into full em-spec.
        #Iterate through mol. lines
        for ai in range(0, len(linedict_list)):
            dhere = linedict_list[ai] #Current mol. line
            #Extend old vel. range by boxwidth and convert to [mu]
            muolds = ((np.concatenate([[-1*boxwidth], dhere["vel"], [boxwidth]])
                    *1.0E9/dhere["freqcen"]) + (cinmu0/1.0/dhere["freqcen"]))
            #Extend edges of old em. flux range to reach full boxwidth
            emolds = np.concatenate([[emonly_list[0]], emonly_list,
                                    [emonly_list[-1]]]])

            #Find indices in full spec. array where this em. will be added
            newinds = np.where(((fullmu_arr <= muolds[1])
                                & (fullmu_arr >= muolds[0])))[0]
            #Interpolate  mol. em. across full spec. x-axis (span of box-width)
            munews = fullmu_arr[newinds]
            interpfunc = interper(x=muolds, y=emolds) #Interpolated function
            emnews = interpfunc(munews, kind="cubic") #Interpolated y-values

            #Add this mol. em. to the overall full em-spectrum
            fullem_arr[newinds] = fullem_arr[newinds] + emnews


        ##Below Section: INTERPOLATE continuum for entire full spec.
        #Sort central continuum values (at vel=0) by ascending central mu
        tempindsort = np.argsort([(cinmu0/1.0/dhere["freqcen"])
                                for dhere in linedict_list])
        tempmusort = np.array([(cinmu0/1.0/dhere["freqcen"])
                                for dhere in linedict_list])[tempindsort]
        tempcontcensort = np.array(contcen_list)[tempmusort] #Matching cen. cont
        #Interpolate full array of continuum (corresponds to full em. array)
        interpfunc = interper(x=tempmusort, y=tempcontcensort) #Interp. func.
        fullcont_arr = interpfunc(fullmu_arr) #Interpolated cont. values


        ##Below Section: COMBINE interp. line and continuum + UNIT CONVERSIONS
        dist = self.get_attr("dist")
        fullem_arr = fullem_arr * 1.0E23/(dist**2) #Em. -> [Jy] at [dist-pc]
        fullcont_arr = fullcont_arr * 1.0E23/(dist**2) #Cont-> [Jy] at [dist-pc]
        fully_arr = fullcont_arr + fullem_arr #Em. and cont. [Jy] at [dist-pc]


        ##Below Section: CONVOLVE emission-spec. and line-spec.
        #Define a gaussian kernel for convolution
        numgauss = np.ceil(3.0*obsres/vres) #Number of kernel points
        toppart = -1*2*((np.arange(0, numgauss) - ((numgauss-1)/2.0))**2)
        botpart = ((obsres/1.0/vres)**2)*np.log(2.0)
        kerngauss = np.exp(toppart/botpart)

        #Convolve emission-spec. and line-spec.
        #Edge behavior: outer vals reflected at input edge to fill missing vals
        resem_arr = ndimage.convolve(fullem_arr, kerngauss, mode="reflect")
        resy_arr = ndimage.convolve(fully_arr, kerngauss, mode="reflect")


        ##Below Section: RESAMPLE at desired resolution of user
        #For emission-spectrum (em. only)
        interpfunc = interper(x=fullmu_arr, y=fullem_arr) #Interp. func.
        outem_arr = interpfunc(outmu_arr) #Interpolated output emission-spec.
        #For line-spectrum (em.+cont.)
        interpfunc = interper(x=fullmu_arr, y=fully_arr) #Interp. func.
        outy_arr = interpfunc(outmu_arr) #Interpolated output line-spec.
        #For continuum (cont. only)
        interpfunc = interper(x=fullmu_arr, y=fullcont_arr) #Interp. func.
        outcont_arr = interpfunc(outmu_arr) #Interpolated output continuum


        ##Below Section: STORE spectra and molecule information + EXIT
        self._set_attr(attrname="flux_spec", attrval=outy_arr) #Line-spec.
        self._set_attr(attrname="flux_em", attrval=outem_arr) #Em-spec.
        self._set_attr(attrname="flux_cont", attrval=outcont_arr) #Continuum
        self._set_attr(attrname="wavelength", attrval=outmu_arr) #Wavelen. [mu]
        self._set_attr(attrname="molinfo", attrval=moldict_list) #All mol. info
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done processing all RADLite output!\n")
        return
    #



    ##CALCULATION METHODS
    #



    ##WRITE METHODS
    def write_fits(self, fitsname):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: LOAD spectra and molecular data
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Starting the write_fits method!")
            print("Saving spectra and molecular info into .fits file...")
        #Try to load spectra data
        try:
            fluxspec_arr = self.get_attr("flux_spec") #Line-spec.
            fluxem_arr = self.get_attr("flux_em") #Em-spec.
            fluxcont_arr = self.get_attr("flux_cont") #Continuum
            wavelen_arr = self._set_attr("wavelength") #Wavelen. [mu]
            moldict_list = self._set_attr("molinfo") #All mol. info
        except AttributeError: #Throw an error if hasn't been processed yet
            raise AttributeError("Whoa! Looks like you haven't processed "
                        +"any RADLite output data yet. You can do so by "
                        +"running the gen_spec() method for this class.")


        ##Below Section: #ASSEMBLE data for .fits file
        #Create primary header for .fits file
        hdr = fitter.Header()
        hdr["WAVEUNIT"] = "micron"
        hdr["FLUXUNIT"] = "Jy"
        hdr["INCL[deg]"] = moldict_list[0]["incl"]
        hdr["DIST[pc]"] = self.get_attr("dist")
        hdu_hdr = fitter.PrimaryHDU(header=hdr)

        #Create and fill table container with flux data
        c1 = fitter.Column(name="Wavelength", array=wavelen_arr, format='D')
        c2 = fitter.Column(name="Emission", array=fluxem_arr, format='D')
        c3 = fitter.Column(name="Spectrum", array=fluxspec_arr, format='D')
        c4 = fitter.Column(name="Continuum", array=fluxcont_arr, format='D')
        hdu_flux = fitter.BinTableHDU.from_columns([c1, c2, c3, c4])

        #Create and fill table container with molecular data
        namelist = [] !!!HERE
        arrlist = []
        formatlist = []
        collist = [fitter.Column(name=namelist[ai], array=arrlist[ai],
                                    format=formatlist[ai])
                    for ai in range(0, len(moldict_list))]
        hdu_mol = fitter.BinTableHDU.from_columns(collist)


        ##Below Section: WRITE assembled data to .fits file + EXIT
        hduall = fitter.HDUList([hdu_hdr, hdu_flux, hdu_mol])
        hduall.writeto(fitsname)
    #
#













#
