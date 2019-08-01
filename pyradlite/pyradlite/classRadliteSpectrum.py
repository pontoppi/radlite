##FILE:
##PURPOSE:


##Below Section: IMPORT necessary functions
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import json
import os.path
from scipy.interpolate import interp1d as interper
from scipy import ndimage as ndimager
from os import listdir as listdirer
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
    @func_timer
    def __init__(self, infilename):
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
        ##Below Section: READ IN + STORE input file
        #Read in input spectrum data
        with open(os.path.join(infilename)) as openfile:
            inputdict = json.load(openfile)
        #Store in secret dictionary (stripping out comments)
        self._attrdict = {}
        for key in inputdict:
            self._attrdict[key] = inputdict[key]["value"]


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Welcome!  You have successfully initialized an instance "
                    +"of RadliteSpectrum(). You can use this instance to "
                    +"process and plot RADLite output. Start by running the "
                    +"gen_spec() method to process RADLite spectra.\n")
        return
    #



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
        ##Below Section: DETERMINE requested attribute
        #Try accessing it
        try:
            return self._attrdict[attrname]
        except KeyError: #If attribute not yet recorded...
            pass

        #Otherwise, raise an error
        raise AttributeError("'"+attrname+"' doesn't seem to be a valid "
                            +"attribute.  Valid attributes are:\n"
                            +str(np.sort([key for key in self._attrdict]))+".\n"
                            +"Run the method gen_spec() (if you haven't "
                            +"yet) to automatically populate more "
                            +"attributes.\n")
    #


    def change_attr(self, attrname, attrval):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: CHANGE value of given attribute
        #NOTE: Currently this acts as a user-friendly wrapper of _set_attr()...
        #...and could be changed if underlying data structures were changed.
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("You called the change_attr() method for "+attrname+"!")
            print("Note that this will override any previous values.")
        #Set new value for given attribute
        self._set_attr(attrname=attrname, attrval=attrval)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("New value has successfully been set!\n")


        ##Below Section: EXIT
        return
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
        ##Below Section: RECORD given attribute under given name + EXIT
        self._attrdict[attrname] = attrval
        return
    #



    ##GENERATIVE METHODS
    @func_timer
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
            print("")


        ##Below Section: READ IN + PROCESS all RADLite output from given runs
        self._read_radliteoutput()
        self._process_spectrum()


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done running the gen_spec() method!\n")
            print("You can use the get_attr() method to access the fluxes "
                    +"and wavelengths directly, as follows:")
            print("get_attr('spectrum') for the full spectrum")
            print("get_attr('emission') for the emission only")
            print("get_attr('continuum') for the continuum")
            print("get_attr('wavelength') for the wavelengths")
            print("You can also plot them using the plot_spec() method, "
                    +"using the same keywords just written above.\n")
        return
    #



    ##OUTPUT DISPLAY METHODS
    @func_timer
    def plot_spec(self, attrname):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        attrval = self.get_attr(attrname)
        wavelen_arr = self.get_attr("wavelength")
        plt.plot(wavelen_arr, attrval)
        plt.show()
    #



    ##READ METHODS
    def _read_core_radliteoutput(self, filedict):
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
        specfilename = filedict["spec"]
        molfilename = filedict["mol"]
        with open(specfilename, 'r') as openfile: #For spectrum
            specdata = openfile.readlines()
        with open(molfilename, 'r') as openfile: #For mol. data
            moldata = openfile.readlines()


        #FOR MOLECULAR DATA FILE
        #Extract molecular traits for this subrun
        molname = moldata[1] #Name of molecule
        molweight = float(moldata[3]) #Weight of molecule
        numlevels = int(moldata[5]) #Number of energy levels
        numtrans = int(moldata[7+numlevels+1]) #Number of transitions

        #Extract transitions
        iloc = 7 + numlevels + 3 #Starting index of transitions
        sechere = moldata[iloc:(iloc+numtrans)] #Section containing all trans.
        sechere = [lhere.split() for lhere in sechere] #Split out spaces
        #For upper and lower degeneracies
        gup_arr = np.array([lh[0] for lh in sechere]).astype(float)
        glow_arr = np.array([lh[2] for lh in sechere]).astype(float)
        #For einstein coefficients, central wavenumbers, and upper energies
        A_arr = np.array([lh[3] for lh in sechere]).astype(float)
        wavenum_arr = np.array([lh[4] for lh in sechere]).astype(float)
        Eup_arr = np.array([lh[5] for lh in sechere]).astype(float)
        #For vibrational and rotational levels
        lvib_arr = np.array([lh[6] for lh in sechere]).astype(float)
        lrot_arr = np.array([lh[7] for lh in sechere]).astype(float)


        #FOR SPECTRUM FILE
        #Extract all spectral traits from spectrum file
        numlines = int(specdata[4]) #Number of mol. lines
        incl = float(specdata[6].split()[2]) #Velocity info
        maxnumpoints = int(specdata[5]) #Max. number of data points per line
        #vwidth = float(velinfo[0])
        #vlsr = float(velinfo[1])
        #incl = float(velinfo[2])

        #Iterate through molecular lines in subrun
        linedict_list = [{} for ai in range(0, numlines)] #To hold all fluxes
        moldict_list = [{} for ai in range(0, numlines)] #To hold all mol. info
        iloc = 6 #Starting index within line files
        for ai in range(0, numlines):
            freqcen = float(specdata[iloc+3]) #Central frequency
            numpoints = int(specdata[iloc+5]) #Num of freq.

            #Extract section of data for current molecular line
            sechere = specdata[(iloc+7):(iloc+7+numpoints)]
            #Convert into 2D array of floats
            sechere = np.array([lhere.split()
                                for lhere in sechere]).astype(float)

            #Record data for this mol. line
            linedict_list[ai]["vel"] = sechere[:,0] #Velocities
            linedict_list[ai]["flux"] = sechere[:,1] #Fluxes
            linedict_list[ai]["freqcen"] = freqcen #Central frequency
            linedict_list[ai]["maxnumpoints"] = maxnumpoints #Max. # data points
            #Increment place in overall file
            iloc = iloc + 7 + numpoints - 1

            #Record information about this mol. line
            #For basic and molecular info
            moldict_list[ai]["molname"] = molname
            moldict_list[ai]["molweight"] = molweight
            moldict_list[ai]["incl"] = incl #Source inclination
            #For upper and lower degeneracies
            moldict_list[ai]["gup"] = gup_arr[ai]
            moldict_list[ai]["glow"] = glow_arr[ai]
            #For einstein coefficient, central wavenumber, and upper energy
            moldict_list[ai]["A"] = A_arr[ai]
            moldict_list[ai]["wavenumcen"] = wavenum_arr[ai]
            moldict_list[ai]["Eup"] = Eup_arr[ai]
            #For vibrational and rotational levels
            moldict_list[ai]["lvib"] = lvib_arr[ai]
            moldict_list[ai]["lrot"] = lrot_arr[ai]


        ##Below Section: RETURN mol. line data + EXIT
        return {"line":linedict_list, "mol":moldict_list}
    #


    @func_timer
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
            print("Extracting all Radlite output from the following list of "
                    +"files: "+str(runpath_list)+"...")

        #Iterate through given RADLite runs
        linedict_list = [None]*len(runpath_list) #List to hold subrun line data
        moldict_list = [None]*len(runpath_list) #List to hold subrun mol. data
        for ai in range(0, len(runpath_list)):
            pathhere = runpath_list[ai] #Path to current run
            #Assemble subrun files (e.g., per core) in this run directory
            suballfiles = listdirer(pathhere)
            subspecfiles = np.sort([os.path.join(pathhere, namehere)
                            for namehere in suballfiles
                            if (namehere.startswith("linespectrum_moldata_")
                                and namehere.endswith(".dat"))])
            submolfiles = np.sort([os.path.join(pathhere, namehere)
                            for namehere in suballfiles
                            if (namehere.startswith("moldata_")
                                and namehere.endswith(".dat"))])

            #Throw an error if not equal number of spectrum and mol. data files
            if len(subspecfiles) != len(submolfiles):
                raise ValueError("Whoa! "+pathhere+" does not contain an equal "
                        +"number of spectrum ('linespectrum_moldata_*.dat') "
                        +"and molecular data ('moldata_*.dat') files. Make "
                        +"you haven't tampered with any of the RADLite output "
                        +"files from running RadliteModel().run_radlite().")

            #Prepare pool of cores to read mol. lines and data for each subrun
            numsubs = len(subspecfiles)
            numcores = self.get_attr("numcores")
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Starting a pool of "+str(numcores)+" cores "
                        +"to extract RADLite output...")
            subfileshere = [{"spec":subspecfiles[bi], "mol":submolfiles[bi]}
                            for bi in range(0, numsubs)] #Files for this subrun
            with mp.Pool(numcores) as ppool:
                subdictshere = ppool.map(self._read_core_radliteoutput,
                                                subfileshere) #Run cores
            linedict_list[ai] = [shere["line"] for shere in subdictshere]
            moldict_list[ai] = [shere["mol"] for shere in subdictshere]
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("The cores have finished extracting RADLite output!\n")


        ##Below Section: Flatten 3D (1=run, 2=subrun, 3=line) lists to 1D lists
        #For mol. line data
        linedict_list = [inlhere for outlhere in linedict_list
                        for inlhere in outlhere] #Flatten to 2D
        linedict_list = [inlhere for outlhere in linedict_list
                        for inlhere in outlhere] #Flatten to 1D
        #For molecular info
        moldict_list = [inlhere for outlhere in moldict_list
                        for inlhere in outlhere] #Flatten to 2D
        moldict_list = [inlhere for outlhere in moldict_list
                        for inlhere in outlhere] #Flatten to 1D
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done processing all RADLite runs!")
            print("There are a total of "+str(len(linedict_list))+" molecular "
                    +"lines across all given RADLite runs.\n")


        ##Below Section: KEEP only unique lines
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Removing any duplicate line occurrences...")
        #Keep only unique lines
        uniqnames_list = [(dhere["molname"]
                            +str(dhere["lvib"])+str(dhere["lrot"])
                            +str(dhere["gup"])+str(dhere["glow"])
                            +str(dhere["Eup"])) for dhere in moldict_list]
        uniqinds = np.unique(uniqnames_list, return_index=True)[1] #Uniq. inds
        linedict_list = [linedict_list[uind]
                            for uind in uniqinds] #Keep only unique mol. lines
        moldict_list = [moldict_list[uind]
                            for uind in uniqinds] #Keep only unique mol. info
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(len(uniqinds))+" molecular lines are being kept.")
            print("Now there are "+str(len(linedict_list))+" lines.\n")


        ##Below Section: CHECK for any incompatible values
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Doing light checks for signs of incompatible sources...")
        #Check that all sources have same inclination
        if len(np.unique([dhere["incl"] for dhere in moldict_list])) > 1:
            raise ValueError("Whoa! The RADLite runs you specified don't all "
                            +"have the same source inclination!")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Light checks passed!\n")


        ##Below Section: STORE RADLite output + EXIT
        self._set_attr(attrname="_radlitelineoutput", attrval=linedict_list)
        self._set_attr(attrname="molinfo", attrval=moldict_list)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done with process of reading RADLite output!\n")
        return
    #



    ##PROCESS METHODS
    @func_timer
    def _process_spectrum(self):
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
        linedict_list = self.get_attr("_radlitelineoutput")
        interpkind = self.get_attr("interpolation")
        obsres = self.get_attr("obsres")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Starting to process read-in RADLite molecular line data...")
            print("NOTE: All interpolation will be done using user-"
                    +"specified interpolation scheme ('"+interpkind+"').")


        ##Below Section: DETERMINE largest encompassing mol. line 'box'
        #NOTE: 'box' refers to x-axis span (e.g., velocity) that will be used...
        #...to encompass each molecular line
        #NOTE: 'vel' = velocity
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Broadening velocity span for each line...")
        #Extract max span of all mol. lines (each line is from -span to +span)
        maxvspan = np.max([np.abs(dhere["vel"][0]) for dhere in linedict_list])
        #Choose box width that is 3*maxvspan OR 3*(specified obs. resolution)
        boxwidth = np.max([(3*maxvspan), (3*self.get_attr("obsres"))])
        #Also extract finest velocity resolution used for mol. lines
        indmax = np.argmax([dhere["maxnumpoints"]
                            for dhere in linedict_list]) #Index of maximum
        vres = linedict_list[indmax]["vel"][1] - linedict_list[indmax]["vel"][0]


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
        maxlen = int(np.floor(np.log(mumax/1.0/mumin)
                            / np.log(growthfrac)) + 1) #mu array size
        fullmu_arr = np.array([(mumin*(growthfrac**ai))
                                for ai in range(0, maxlen)])
        #NOTE: Used compound interest equation to calculate maxlen, fullmu_arr
        fullem_arr = np.zeros(maxlen) #Initialized; will hold full em-spec.

        #Prepare wavelen. [mu] array for output spectrum at desired resolution
        userres = self.get_attr("vsampling") #User-specified sampling res.
        userfrac = (1 + (userres/1.0/cinkm0)) #Fraction by which mu will grow
        userlen = int(np.floor(np.log(mumax/1.0/mumin)
                            / np.log(userfrac)) + 1) #mu array size
        outmu_arr = np.array([(mumin*(userfrac**ai))
                                for ai in range(0, userlen)])


        ##Below Section: SUBTRACT out continuum for each mol. line
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Interpolating and subtracting out continuum...")
        emonly_list = [None]*len(linedict_list) #Init. array for em. only
        contcen_list = [None]*len(linedict_list) #Init. array for center cont.
        #Iterate through mol. lines
        for ai in range(0, len(linedict_list)):
            dhere = linedict_list[ai] #Current mol. line
            #Linearly interpolate continuum from line edges; calculate at vel=0
            tempys = [dhere["flux"][0], dhere["flux"][-1]] #Fluxes to interp.
            tempxs = [dhere["vel"][0], dhere["vel"][-1]] #Vels. to interp.
            interpfunc = interper(x=tempxs, y=tempys, kind=interpkind)
            #Record central continuum and isolated emission for this line
            emonly_list[ai] = dhere["flux"] - interpfunc(dhere["vel"]) #Em. only
            contcen_list[ai] = interpfunc(0) #Continuum at central vel. (vel=0)


        ##Below Section: INTERPOLATE + ADD each mol. emission into full em-spec.
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Combining molecular lines to form full spectrum...")
        #Iterate through mol. lines
        for ai in range(0, len(linedict_list)):
            dhere = linedict_list[ai] #Current mol. line
            #Extend old vel. range by boxwidth and convert to [mu]
            muolds = ((np.concatenate([[-1*boxwidth], dhere["vel"], [boxwidth]])
                    *1.0E9/dhere["freqcen"]) + (cinmu0/1.0/dhere["freqcen"]))
            #Extend edges of old em. flux range to reach full boxwidth
            emolds = np.concatenate([[emonly_list[ai][0]], emonly_list[ai],
                                    [emonly_list[ai][-1]]])

            #Find indices in full spec. array where this em. will be added
            newinds = np.where(((fullmu_arr <= muolds[-1])
                                & (fullmu_arr >= muolds[0])))[0]
            #Interpolate  mol. em. across full spec. x-axis (span of box-width)
            munews = fullmu_arr[newinds]
            interpfunc = interper(x=muolds, y=emolds,
                                    kind=interpkind) #Interp func
            emnews = interpfunc(munews) #Interpolated y-values

            #Add this mol. em. to the overall full em-spectrum
            fullem_arr[newinds] = fullem_arr[newinds] + emnews


        ##Below Section: INTERPOLATE continuum for entire full spec.
        #Sort central continuum values (at vel=0) by ascending central mu
        tempindsort = np.argsort([(cinmu0/1.0/dhere["freqcen"])
                                for dhere in linedict_list])
        tempmusort = np.array([(cinmu0/1.0/dhere["freqcen"])
                                for dhere in linedict_list])[tempindsort]
        tempcontcensort = np.array(contcen_list)[tempindsort] #Corres. cen. cont
        #Linearly stretch continuum to reach the edges of the spectrum
        tempmusort = np.concatenate([[mumin], tempmusort, [mumax]])
        tempcontcensort = np.concatenate([[tempcontcensort[0]], tempcontcensort,
                                            [tempcontcensort[-1]]])
        #Interpolate full array of continuum (corresponds to full em. array)
        interpfunc = interper(x=tempmusort, y=tempcontcensort,
                                kind=interpkind) #Interp. func.
        fullcont_arr = interpfunc(fullmu_arr) #Interpolated cont. values


        ##Below Section: COMBINE interp. line and continuum + UNIT CONVERSIONS
        dist = self.get_attr("dist")
        fullem_arr = fullem_arr * 1.0E23/(dist**2) #Em. -> [Jy] at [dist-pc]
        fullcont_arr = fullcont_arr * 1.0E23/(dist**2) #Cont-> [Jy] at [dist-pc]
        fully_arr = fullcont_arr + fullem_arr #Em. and cont. [Jy] at [dist-pc]

        ##Below Section: CONVOLVE emission-spec. and line-spec.
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Convolving emission and (emission+continuum) to user-"
                    +"specified observation resolution...")
        #Define a gaussian kernel for convolution
        numgauss = np.ceil(3.0*obsres/vres) #Number of kernel points
        toppart = -1*2*((np.arange(0, numgauss) - ((numgauss-1)/2.0))**2)
        botpart = ((obsres/1.0/vres)**2)*np.log(2.0)
        kerngauss = np.exp(toppart/botpart)

        #Convolve emission-spec. and line-spec.
        #Edge behavior: outer vals reflected at input edge to fill missing vals
        resem_arr = (ndimager.convolve(fullem_arr, kerngauss, mode="reflect")
                        /1.0/np.sum(kerngauss)) #Convolved, norm. by kernel sum
        resy_arr = (ndimager.convolve(fully_arr, kerngauss, mode="reflect")
                        /1.0/np.sum(kerngauss)) #Convolved, norm. by kernel sum


        ##Below Section: RESAMPLE at desired resolution of user
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Resampling to user-specified spectrum resolution...")
        #For emission-spectrum (em. only)
        interpfunc = interper(x=fullmu_arr, y=resem_arr) #Interp. func.
        outem_arr = interpfunc(outmu_arr) #Interpolated output emission-spec.
        #For line-spectrum (em.+cont.)
        interpfunc = interper(x=fullmu_arr, y=resy_arr) #Interp. func.
        outy_arr = interpfunc(outmu_arr) #Interpolated output line-spec.
        #For continuum (cont. only)
        interpfunc = interper(x=fullmu_arr, y=fullcont_arr) #Interp. func.
        outcont_arr = interpfunc(outmu_arr) #Interpolated output continuum


        ##Below Section: STORE spectra and molecule information + EXIT
        self._set_attr(attrname="spectrum", attrval=outy_arr) #Line-spec.
        self._set_attr(attrname="emission", attrval=outem_arr) #Em-spec.
        self._set_attr(attrname="continuum", attrval=outcont_arr) #Continuum
        self._set_attr(attrname="wavelength", attrval=outmu_arr) #Wavelen. [mu]
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done processing RADLite spectrum!\n")
        return
    #



    ##WRITE METHODS
    @func_timer
    def write_fits(self, fitsname, overwrite):
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
            print("Saving spectra and molecular info in a .fits file...")
        #Try to load spectra data
        try:
            fluxspec_arr = self.get_attr("spectrum") #Line-spec.
            fluxem_arr = self.get_attr("emission") #Em-spec.
            fluxcont_arr = self.get_attr("continuum") #Continuum
            wavelen_arr = self.get_attr("wavelength") #Wavelen. [mu]
            moldict_list = self.get_attr("molinfo") #All mol. info
        except AttributeError: #Throw an error if hasn't been processed yet
            raise AttributeError("Whoa! Looks like you haven't processed "
                        +"any RADLite output data yet. You can do so by "
                        +"running the gen_spec() method for this class.")


        ##Below Section: #ASSEMBLE data for .fits file
        #Create primary header for .fits file
        hdr = fitter.Header()
        hdr["WAVEUNIT"] = "micron"
        hdr["FLUXUNIT"] = "Jy"
        hdr["INCL_deg"] = moldict_list[0]["incl"]
        hdr["DIST_pc"] = self.get_attr("dist")
        hdu_hdr = fitter.PrimaryHDU(header=hdr)

        #Create and fill table container with flux data
        c1 = fitter.Column(name="Wavelength", array=wavelen_arr)#, format='D')
        c2 = fitter.Column(name="Emission", array=fluxem_arr)#, format='D')
        c3 = fitter.Column(name="Spectrum", array=fluxspec_arr)#, format='D')
        c4 = fitter.Column(name="Continuum", array=fluxcont_arr)#, format='D')
        hdu_flux = fitter.BinTableHDU.from_columns([c1, c2, c3, c4])

        #Create and fill table container with molecular data
        namelist = [key for key in moldict_list[0]] #Names of mol. data
        arrlist = [np.array([dhere[key] for dhere in moldict_list[0]])
                    for key in namelist] #Arrays corresponding to data names
        collist = [fitter.Column(name=namelist[ai], array=arrlist[ai])#,
                                    #format=formatlist[ai])
                    for ai in range(0, len(moldict_list))] #Data columns
        hdu_mol = fitter.BinTableHDU.from_columns(collist) #Data container


        ##Below Section: WRITE assembled data to .fits file + EXIT
        hduall = fitter.HDUList([hdu_hdr, hdu_flux, hdu_mol])
        hduall.writeto(fitsname, overwrite=overwrite)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Data has been successfully saved to "+fitsname+"!")
        return
    #
#












#
