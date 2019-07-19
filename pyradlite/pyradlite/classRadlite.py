##FILE:
##PURPOSE:


##Below Section: IMPORT necessary functions
import subprocess
import multiprocessing as mp
import numpy as np
import json
import os.path
from datetime import datetime as dater
from classUtils import func_timer
import astropy.constants as const
import radmc

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



##
class Radlite():
    @func_timer
    def __init__(self, infilepath="./"):
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
        ##Below Section: READ IN + STORE input files
        #Read in input RADLite data
        with open(os.path.join(infilepath, "input_radlite.json")) as openfile:
            self._valdict = json.load(openfile)
        #Read in HITRAN data and then extract desired molecule
        with open(os.path.join(infilepath, "data_hitran.json")) as openfile:
            hitrandict = json.load(openfile)
        #Store molecular data for specified molecule
        moldict = hitrandict[self._valdict["molname"]] #Extract molecular data
        for key in moldict: #Store molecule-specific data
            self._valdict[key] = moldict[key]
        self._valdict["molmass"] = moldict["molweight"]*amu0 #Tack on mol. mass


        ##Below Section: PRINT a welcome message, if so desired
        if self.get_value("verbose"):
            print("--------------------------------------------------")
            print("Welcome to RADLite Version 1.2, wrapped in Python.")
            print("")
            print("RADLite Version 1.2 was written by:")
            print("Klaus Pontoppidan (pontoppi@stsci.edu)")
            print("Kees Dullemond")
            print("Alex Lockwood")
            print("Rowin Meijerink")
            print("--------------------------------------------------")
            print("")
            print("")


        ##Below Section: CHECK inputs for any user error; otherwise record them
        #Make sure that desired image cube output is valid
        validimage = ["spec", "circ"]
        if self._valdict["image"] not in validimage:
            raise ValueError("Sorry, the image you have chosen ("
                    +str(inpdict["image"])+") is not a valid image.  The value "
                    +"of image must be one of the following: "
                    +str(validimage)+".")
        ##Below Section: CHOOSE RADLite exec. based on desired output
        #Prepare executable for desired image cube output
        if self._valdict["image"] == "spec": #If desired cube output is spectrum
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Will prepare a spectrum-formatted image cube...")
                print("")
            self._valdict["executable"] = self._valdict["exe_path"]+"RADlite"
        elif self._valdict["image"] == "circ": #If desired output is circular
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Will prepare a circular-formatted image cube...")
                print("")

        #Make sure that desired LTE format is valid
        if not isinstance(self._valdict["lte"], bool):
            raise ValueError("Sorry, the lte you have chosen ("
                    +str(self._valdict["lte"])+") is not a valid lte.  The "
                    +"value of lte must be a boolean (True or False).")


        ###TEMP. BLOCK START
        ##Below Section: TEMPORARY - USE PYRADMC TO READ IN RADMC MODEL
        modradmc = radmc.radmc_model("./")
        self._valdict["radius"] = modradmc.radius
        self._valdict["theta"] = modradmc.theta
        self._valdict["dusttemperature"] = modradmc.dusttemperature[:,:,0].T
        #NOTE: !!! Assuming one species read in; hence the [:,:,0]
        ###TEMP. BLOCK END



        ##Below Section: EXIT
        return
    #



    ##STORE AND FETCH METHODS
    #@func_timer
    def get_value(self, valname):
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
        ##Below Section: DETERMINE requested value
        #Try accessing it
        try:
            return self._valdict[valname]
        except KeyError: #If value not yet recorded...
            pass

        #If that doesn't work, try calculating it
        try:
            eval("self._calc_"+valname+"()") #Try calculating it
            return self._valdict[valname]
        except AttributeError: #If value not calculable...
            pass

        #If that doesn't work, try reading it in
        try:
            eval("self._read_"+valname+"()") #Try reading it in
            return self._valdict[valname]
        except AttributeError: #If value not readable
            raise NameError("'"+valname+"' doesn't seem to be a valid value.  "
                            +"Valid values are:\n"
                            +str([key for key in self._valdict])+".\n"
                            +"Run the method run_radlite() (if you haven't "
                            +"yet) to automatically populate more values.\n"
                            +"Alternatively, you can pass the name of a "
                            +"supported physics component (e.g., 'velocity') "
                            +"to the get_value() method to populate that "
                            +"component on its own.")
    #


    def _get_value_perprocess(self, valname, pind):
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
        return self._valdict[valname][self.get_value("_splitinds")[pind]]
    #


    def _set_value(self, valname, val):
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
        ##Below Section: RECORD given value under given name + EXIT
        self._valdict[valname] = val
        return
    #


    @func_timer
    def run_radlite(self):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: SET UP lists of initial and input files
        initfilelist = ["problem_params.pro", "line_params.ini"]
        cpfilelist = ["radius.inp",
                        "theta.inp", "frequency.inp", "density.inp",
                        "dustdens.inp", "dusttemp_final.dat",
                        "dustopac.inp", "dustopac_*.inp", "abundance.inp",
                        "temperature.inp", "velocity.inp"]


        ##Below Section: WRITE RADLite input files
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Writing radlite input files...")
        self._write_abundanceinp() #Abundance
        self._write_turbulenceinp() #Turbulence
        self._write_velocityinp() #Velocity


        ##Below Section: SET UP result directory
        rundir = os.path.join(self.get_value("run_dir"), self.get_value("run_name"))
        if self.get_value("dodate"): #If date should be appended to dir. name
            currtime = dater.now()
            timestamp = "_{0}_{1d}h{2d}m{3d}s".format(currtime.date,
                            currtime.hour, currtime.minute, currtime.second)
            rundir += timestamp

        #Make the desired directory, if nonexistent
            try:
                comm = subprocess.call(["mkdir", rundir]) #Create directory
            except (comm != 0): #Will override files within, otherwise
                pass
        if self.get_value("verbose"): #Verbal output, if so desired
            print("All input files and final data will be saved to the "
                    +"following directory: "+rundir)

        #Copy initial files into final directory
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Copying over initial data files into "+rundir+"...")
        for iname in initfilelist:
            comm = subprocess.call(["cp",
                            os.path.join((self.get_value("inp_path"),
                                            iname)),
                            os.path.join((rundir, "/"))])


        ##Below Section: PROCESS data from the LTE or NLTE data file
        #Read in either LTE or NLTE data
        if self.get_value("lte"): #If LTE treatment desired
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Extracting LTE molecular data...")
                print("")
            self._read_hitran(numcores)
        else: #Else, if non-LTE treatment desired
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Extracting NLTE molecular data...")
                print("")
            molfile = "molfile.dat"
            self._read_lambda(numcores)


        ##Below Section: RUN RADLITE on single/multiple cores in subfolders
        numcores = self.get_value("ncores")
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Running RADLite on "+str(numcores)+" core(s)...")
        #Prepare pool of cores
        plist = []
        for ai in range(0, numcores):
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Prepping "+str(ai)+"th core...")

            #Make a directory for this run
            cpudir = "./workingdir_cpu"+str(ai)+"/" #core-specific dir.
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Generating directory: "+cpudir)
            try:
                comm = subprocess.call(["mkdir", cpudir]) #Create subdirectory
            except (comm != 0): #If directory already exists, replace it
                comm = subprocess.call(["rm", "-r", cpudir]) #Erase prev. dir.
                comm = subprocess.call(["mkdir", cpudir]) #New empty dir.
            #Copy initial files into this new core subdirectory
            for iname in initfilelist:
                comm = subprocess.call(["cp", rundir+"/"+iname, cpudir])

            #Call core routine
            phere = mp.Process(target=self._run_core,
                                    args=(pind, cpudir, rundir,))
            plist.append(phere)
            #Start process
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Starting "+str(ai)+"th core in "+cpudir+"...")
            phere.start()

        #Close pool of cores
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Done running core(s)!")
        for ai in range(0, numcores):
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Closing "+str(ai)+"th core...")
            plist[ai].join()


        ##Below Section: FINISH and EXIT
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Done running run_radlite()!")
        return
    #


    @func_timer
    def _run_core(self, pind, cpudir, rundir):
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
        ##Below Section: GENERATE radlite input files
        #NOTE: SOME OF THESE ARE GENERAL ENOUGH TO NOT BE DONE PER core !!!!!!!!!!
        self._write_radliteinp(cpudir) #Direct radlite input file
        self._write_levelpopinp(cpudir) #Level population input file
        self._write_moldatadat(filepathandname=cpudir, pind=pind,
            outfilename=(cpudir+"moldata_"+str(pind)+".dat")) #Mol. datafiles


        ##Below Section: CHECK that all files are in order
        #BELOW FROM IDL !!! - read_vel, 'velocity.inp', vel
        #Check passband width
        #velmax = np.max(np.abs(vel.vphi))
        #if self.get_value("verbose"): #Verbal output, if so desired
        #    print("Max. velocity = {0:.2f}km/s".format(velmax/1E5))


        ##Below Section: RUN RADLITE
        logfile = open(cpudir+"RADLITE_core"+str(pind)+".log", 'w')
        comm = subprocess.call([self.executable], #Call RADLite
                                cwd=cpudir+"/", #Call within core subdir.
                                stdout=logfile) #Send output to log file


        ##Below Section: EXTRACT output files
        comm = subprocess.call(["mv",
                            cpudir+"/moldata_"+str(pind)+".dat",
                            rundir+"/"]) #Molecular data used by this core
        comm = subprocess.call(["mv",
                            cpudir+"/linespectrum_moldata_"+str(pind)+".dat",
                            rundir+"/"]) #Line spectrum output
        if self.image == 2:
            comm = subprocess.call(["mv",
                            cpudir+"/lineposvelcirc_moldata_"+str(pind)+".dat",
                            rundir+"/"]) #Circular 3D image cube output


        ##Below Section: DELETE subdir. + EXIT
        comm = subprocess.call(["rm", "-r", cpudir]) #Erase core dir.
        return
    #



    ##READ METHODS
    @func_timer
    def _read_hitran(self):
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
        ##Below Section: PROCESS all lines in file
        filepathandname = os.path.join(
                                    self.get_value("hit_path"),
                                    self.get_value("hitran_file"))
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Extracting molecular lines from file "+filepathandname+"...")
        #Read in all filelines
        with open(filepathandname, 'r') as openfile:
            alllines = openfile.readlines()
        if self.get_value("verbose"): #Verbal output, if so desired
            print("There are "+str(len(alllines))+" molecular lines in total.")

        #Set HITRAN entries and number of characters per HITRAN file entry
        hitranchars = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15, 15, 15,
                            6, 12, 1, 7, 7] #Characters per HITRAN entry
        hitrannames = ["mol", "iso", "wavenum", "ins", "A", "air", "self",
                            "Elow", "temp", "press", "vup", "vlow",
                            "qup", "qlow", "err", "ref", "flag", "gup", "glow"]
        hitrandtypes = [int]*2 +[float]*8 +[str]*4 +[int]*2 +[str] +[float]*2

        #Record *all* molecular line into a convenient dictionary
        #NOTE: Any undesired lines will be removed later
        hitrandict = {} #Initialize dictionary to hold line values
        ii = 0 #Index to keep track of place within entries
        for ai in range(0, len(hitrannames)):
            namehere = hitrannames[ai] #Name of current entry
            lenhere = hitranchars[ai] #Character length of current entry
            datahere = np.array([linehere[ii:ii+hitranchars[ai]]
                                    for linehere in alllines], #All entries
                                    dtype=hitrandtypes[ai]) #Entry datatype
            hitrandict[namehere] = datahere #Record the extracted entry data
            ii = ii + hitranchars[ai] #Update place within fileline characters
        #Also calculate and record upper energy levels
        hitrandict["Eup"] = hitrandict["Elow"] + hitrandict["wavenum"] #E_up

        #Throw error if any signs of incomplete data
        if 0 in hitrandict["gup"]:
            raise ValueError("Something's wrong!  Incomplete molecular "
                        +"HITRAN data at wavenum="
                        +str(hitrandict["wavenum"][
                                                hitrandict["gup"]==0])+".")


        ##Below Section: REMOVE any lines outside of desired range
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Removing molecular lines outside of specified criteria...")
        #Extract indices of lines that fall within criteria
        keepinds = np.where( ~(
                    #If incorrect isotopologue
                    (hitrandict["isonum"] != self.isotop)
                    #If wavenumber outside of desired wavenumber range
                    | (hitrandict["wavenum"] < self.wavenumrange[0])
                    | (hitrandict["wavenum"] > self.wavenumrange[1])
                    #If intensity, level, or upper energy beyond given cutoffs
                    | (hitrandict["ins"] < self.cutoff)
                    | (hitrandict["Eu"] > self.Eupmax)
                    | (float(hitrandict["vup"]) > self.vupmax)))[0]
        numlines = np.sum(keepinds) #Number of lines that fall within criteria

        #Keep only those lines that fall within criteria; delete all other lines
        for key in hitrandict:
            hitrandict[key] = hitrandict[key][keepinds]
        self.hitrandict = hitrandict #Record final set of molecular lines
        self.numlines = numlines #Record final count of lines
        if self.get_value("verbose"): #Verbal output, if so desired
            print("There are "+str(numlines)+" molecular lines "
                    +"that fall within specified criteria.")


        ##Below Section: SPLIT data across given number of cores
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Dividing up the lines for "
                        +self.numcores+" cores...")
        #Determine indices for splitting up the data
        numpersplit = self.numlines // self.numcores #No remainder
        splitinds = [[(ai*numpersplit),((ai+1)*numpersplit)]
                        for ai in range(0, self.numcores)] #Divide indices
        splitinds[-1][1] = self.numlines #Tack leftovers onto last core
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Here are the chosen line intervals per core:")
            print([("core "+str(ehere[0])+": Interval "+str(ehere[1]))
                                    for ehere in enumerate(splitinds)])
        self._set_value(valname="_splitinds", val=splitinds) #Record split lines


        ##Below Section: EXIT
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Done extracting molecular data!")
        return
    #


    @func_timer
    def _read_starinfo(self):
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
        ##Below Section: PROCESS all lines in file
        filepathandname = os.path.join(
                                    self.get_value("inp_path"), "starinfo.inp")
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Extracting molecular lines from file "+filepathandname+"...")
        #Read in all filelines
        with open(filepathandname, 'r') as openfile:
            alllines = openfile.readlines()


        ##Below Section: STORE read-in data + EXIT
        rstar = float(alllines[1]) #Stellar radius
        mstar = float(alllines[2]) #Stellar mass
        teff = float(alllines[3]) #Stellar eff. temperature
        starinfodict = {"mstar":mstar, "rstar":rstar, "teff":teff}
        #Store together and individually
        self._set_value(valname="starinfo", val=starinfodict)
        self._set_value(valname="mstar", val=mstar)
        self._set_value(valname="rstar", val=rstar)
        self._set_value(valname="teff", val=teff)
        return
    #



    ##CALCULATION METHODS
    ###NOTE: MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING METHODS
    @func_timer
    def _calc_abundance(self):
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
        ##Below Section: EXTRACT structural information
        temparr = self.get_value("dusttemperature")
        rlen = len(self.get_value("radius"))
        tlen = len(self.get_value("theta"))
        gamval = self.get_value("gamma")
        abundrange = [self.get_value("min_abun"), self.get_value("max_abun")]
        temp_fr = self.get_value("temp_fr")
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Calculating abundance...")


        ##Below Section: CALCULATE abundance based on specified abundance mode
        if temp_fr is False: #For constant abundance
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Setting a constant abundance...")
            abundarr = np.ones(shape=(tlen, rlen))*abundrange[1]
            self._set_value(valname="abund_collpartner", val=0.0)
        else: #For abundance that changes below freeze-out temperature
            if self.get_value("verbose"): #Verbal output, if so desired
                print("Setting abundance to min. below {0:f}K..."
                                .format(temp_fr))
            abundarr = np.ones(shape=(tlen, rlen))*abundrange[1]
            abundarr[temparr < temp_fr] = abundrange[0]
            self._set_value(valname="collpartner", val=0.0)


        ##NOTE: IN MAKE_ABUNDANCE.PRO, THERE WAS PART HERE ABOUT MOL_DESTRUCT

        ##Below Section: RECORD calculated abundance + EXIT
        self._set_value(valname="abundance", val=abundarr)
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Done calculating abundance!")
        return
    #


    ###NOTE: _CALC_TURBULENCE WILL BE MOVED TO BASEMODEL CLASS; MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING METHODS
    @func_timer
    def _calc_turbulence(self):
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
        ##Below Section: EXTRACT gas information
        temparr = self.get_value("dusttemperature")
        gamval = self.get_value("gamma")
        muval = self.get_value("mu")
        weightval = self.get_value("molweight")
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Calculating turbulence...")
            print("Using first dust component temperature to determine "
                    +"turbulent velocities...")


        ##Below Section: CALCULATE turbulence (alpha-model) using sound speed
        csarr = np.sqrt(gamval*kB0*temparr/1.0/(muval*mp0)) #Sound speed
        turbarr = self.get_value("alpha")*csarr #Turbulence
        #Add thermal broadening in quadrature
        thermbroadarr = np.sqrt(2.0*kB0*temparr/(weightval*mp0)) #Therm. broad.
        turbarr = np.sqrt((turbarr**2) + (thermbroadarr**2)) #Updated turbulence


        ##Below Section: RECORD calculated turbulence + EXIT
        self._set_value(valname="turbulence", val=turbarr)
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Done calculating turbulence!")
        return
    #


    ###NOTE: _CALC_VELOCITY WILL BE MOVED TO BASEMODEL CLASS; MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING _CALC_VELOCITY METHODS
    @func_timer
    def _calc_velocity(self):
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
        ##Below Section: EXTRACT stellar information
        rarr = self.get_value("radius") #Radii
        rlen = len(rarr) #Number of radius points
        tlen = len(self.get_value("theta"))//2 #Number of theta points
        rexparr = np.resize(rarr, (tlen, rlen)) #Radii, expended to 2D
        mstar = self.get_value("starinfo")["mstar"] #Stellar mass
        rstar = self.get_value("starinfo")["rstar"] #Stellar radius
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Calculating velocity field...")
            print("Used starinfo.inp file for mstar and rstar...")


        ##Below Section: CALCULATE Keplerian velocity
        vdict = {} #Dictionary to hold different velocity dimensions
        vdict["r"] = np.zeros(shape=(tlen, rlen)) #Radial velocity
        vdict["th"] = np.zeros(shape=(tlen, rlen)) #Theta velocity
        vdict["phi"] = np.sqrt(G0*mstar/1.0/rexparr) #Phi velocity


        ##Below Section: RECORD velocity and EXIT
        self._set_value("velocity", vdict) #Record velocity
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Done calculating velocity.")
        return
    #



    ##WRITE METHODS
    @func_timer
    def _write_abundanceinp(self):
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
        ##Below Section: BUILD string containing abundance information
        #Extract abundance
        abundarr = self.get_value("abundance") #Abundance data
        collpartner = self.get_value("collpartner")
        rlen = len(self.get_value("radius")) #Length of radius array
        tlen = len(self.get_value("theta"))//2 #Half-length of theta array
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\n".format(rlen, tlen)
        #Fill in string with abundance information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8e}\t{1:.8f}\n".format(abundarr[ti, ri],
                                                collpartner)

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_value("inp_path"), "abundance.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    @func_timer
    def _write_gasdensityinp(self, cpudir):
        pass
    #


    @func_timer
    def _write_levelpopinp(self, cpudir):
        pass
    #


    @func_timer
    def _write_linespectruminp(self, numlines, molfilename, outfilename):
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
        ##Below Section: BUILD string to form the spectrum input file
        writestr = ""
        writestr += "{0:<8d}Format number\n".format(1)
        writestr += "{0:<8d}Spectrum output style\n".format(1)
        writestr += ("-"*63)+"\n"
        writestr += "{0:<8d}Format number\n".format(2)
        writestr += ("-"*63)+"\n"
        writestr += "{0:<8.2f}{1:<50s}\n".format(
                            self.passband, "Width of line passband [km/s]")
        writestr += "{0:<8.2f}{1:<50s}\n".format(
                            self.vsampling, "Velocity sampling [km/s]")
        writestr += ("-"*65)+"\n"
        writestr += "{0:<8d}Format number\n".format(2)
        writestr += "{0:<8d}{1:<50s}\n".format(
                            self.image,
                            "Command (0=spectrum, 2=image[3-D P/V cube])")
        writestr += "{0:<8.1f}{1:<50s}\n".format(
                            self.dist, "Distance in [pc]")
        writestr += "{0:<8.1f}{1:<50s}\n".format(
                            self.incl, "Inclination [deg]")
        writestr += "{0:<8.1f}{1:<60s}\n".format(
                            self.vlsr, "Radial velocity, rel. to local "
                                        +"standard of rest [km/s]")
        writestr += "{0:<8d}{1:<50s}\n".format(
                            numlines, "Nr of lines to make spectrum/image")
        writestr += "{0:<8d}Starting line to make spectrum/image\n".format(1)
        if self.image == 2: #For 3D image cube
            npix = int(np.ceil(self.imwidth_au/1.0/self.ssampling))
            writestr += "{0:<8d}{1:<50s}\n".format(npix, "Nr of x pixels")
            writestr += "{0:<8d}{1:<50s}\n".format(npix, "Nr of y pixels")
            writestr += "{0:<8d}image size in cm?\n".format(1)
            writestr += "{0:<8.1e}{1:<50s}\n".format(
                            self.imwidth_cm, "size x direction")
            writestr += "{0:<8.1e}{1:<50s}\n".format(
                            self.imwidth_cm, "size y direction")
            writestr += "{0:<8d}Phi offset?\n".format(0)
            writestr += "{0:<8d}x offset?\n".format(0)
            writestr += "{0:<8d}y offset?\n".format(0)
            writestr += "{0:<8d}add star?\n".format(1)


        ##Below Section: WRITE the results to file + EXIT function
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    @func_timer
    def _write_moldatadat(self, filepathandname, pind, outfilename):
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
        #Below Section: EXTRACT lists of info for line transitions
        molname = self._get_value_perprocess("molname", pind=pind)
        molweight = self._get_value_perprocess("molweight", pind=pind)
        Euplist = self._get_value_perprocess("Eup", pind=pind)
        Elowlist = self._get_value_perprocess("Elow", pind=pind)
        wavenumlist = self._get_value_perprocess("wavenum", pind=pind)
        Alist = self._get_value_perprocess("A", pind=pind)
        vuplist = self._get_value_perprocess("vup", pind=pind)
        vlowlist = self._get_value_perprocess("vlow", pind=pind)
        quplist = self._get_value_perprocess("qup", pind=pind)
        qlowlist = self._get_value_perprocess("qlow", pind=pind)


        ##Below Section: COMBINE levels + REMOVE duplicates to get unique levels
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Writing "+outfilename+"...")
            print("Counting up unique levels for "+outfilename+"...")
        #Combine energies, transitions, and degeneracies
        Ealllist = np.concatenate((Elowlist, Euplist))
        valllist = np.concatenate((vlowlist, vuplist))
        qalllist = np.concatenate((qlowlist, quplist))
        galllist = np.concatenate((glowlist, guplist))
        #Sort combined lists by energies
        sortinds = np.argsort(Ealllist)
        Ealllist = Eallist[sortinds]
        valllist = vallist[sortinds]
        qalllist = qallist[sortinds]
        galllist = gallist[sortinds]

        #Extract indices for unique levels only
        uniqEinds = np.where( ~( #Indices of unique Elow, Eu value combinations
                    (np.abs(Ealllist[0:-1]-Ealllist[1:len(Ealllist)])
                            /1.0/(Ealllist[1:len(Ealllist)]+0.1)) < 0.0001))[0]
        uniqvinds = np.where( ~( #Indices of unique vlow, vup combinations
                        vuplist[0:-1] == vuplist[1:len(Ealllist)]))[0]
        uniqginds = np.where( ~( #Indices of unique glow, gup combinations
                        glowlist[0:-1] == glowlist[1:len(Ealllist)]))[0]

        #Combine and apply indices
        if self.get_value("verbose"): #Verbal output, if so desired
            print("Removing any duplicate levels for "+outfilename
                    +"...  To start, there are "+str(len(Ealllist))+" levels.")
        uniqinds = np.unique(np.concatenate((uniqEinds, uniqvinds, uniqginds)))
        Ealllist = Ealllist[uniqinds]
        valllist = valllist[uniqinds]
        qalllist = qalllist[uniqinds]
        galllist = galllist[uniqinds]
        if self.get_value("verbose"): #Verbal output, if so desired
            print(str(np.sum(~uniqinds))+" levels have been removed for "
                        +outfilename+", leaving "+str(len(Eallist))
                        +" unique levels.")


        ##Below Section: BUILD string to form the molecular data file
        writestr = ""
        writestr += "!MOLECULE\n{0:s}\n".format(molname)
        writestr += "!MOLECULAR WEIGHT\n{0:4.1f}\n".format(molweight)
        writestr += "!NUMBER OF ENERGY LEVELS\n{0:6d}\n".format(len(Ealllist))
        #Tack on unique levels, energies, and degeneracies
        writestr += "!LEVEL + ENERGIES(cm^-1) + WEIGHT + v + Q\n"
        for ai in range(0, len(Ealllist)):
            writestr += "{0:5d}{1:12.4f}{2:7.1f}{3:>15s}{4:>15s}\n".format(
                            (ai+1), Ealllist[ai], galllist[ai], valllist[ai],
                            qalllist[ai])
        #Tack on transitions
        writestr += "!NUMBER OF RADIATIVE TRANSITIONS\n"
        writestr += "{0:6d\n}".format(len(wavenumlist))
        writestr += ("!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm^-1) + "
                        +"E_u(cm^-1) + v_l + Q_p + Q_pp\n")
        for ai in range(0, len(wavenumlist)):
            levu = np.where((np.abs(Ealllist - Euplist[ai])
                                /1.0/Euplist[ai]) < 1E-4)[0]
            if Elowlist[ai] != 0: #If not down to 0-level
                levl = np.where((np.abs(Ealllist - Elowlist[ai])
                                /1.0/Elowlist[ai]) < 1E-4)[0]
            else:
                levl = np.where(Ealllist == 0)[0]
            writestr += ("{0:5d}{1:5d}{2:5d}{3:12.3e}{4:16.7f}".format(
                                (ai+1), levu[0]+1, levl[0]+1, Alist[ai],
                                wavenumlist[ai])
                        +"{5:12.5f}{6:>15}{7:>15}{8:>15}{9:>15}\n".format(
                                Euplist[ai], vuplist[ai], vlowlist[ai],
                                quplist[ai], qlowlist[ai]))


        ##Below Section: WRITE the results to file + EXIT function
        with open(outputfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    @func_timer
    def _write_radliteinp(self, cpudir):
        pass
    #


    @func_timer
    def _write_turbulenceinp(self):
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
        ##Below Section: BUILD string containing turbulence information
        #Extract turbulence
        turbarr = self.get_value("turbulence") #Turbulence data
        rlen = len(self.get_value("radius")) #Length of radius array
        tlen = len(self.get_value("theta"))//2 #Length of theta array
        #Set up string
        writestr = "1\n" #Initialize string
        writestr += "{0:d}\t{1:d}\n".format(rlen, tlen)
        #Fill in string with turbulence information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8f}\n".format(
                                        (turbarr[ti, ri]/1.0E5)) #[cm/s]->[km/s]

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_value("inp_path"), "turbulence.inp")
        with  open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    @func_timer
    def _write_velocityinp(self):
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
        ##Below Section: BUILD string containing velocity information
        #Extract velocity
        velarr = self.get_value("velocity") #Velocity data
        rlen = len(self.get_value("radius")) #Length of radius array
        tlen = len(self.get_value("theta"))//2 #Length of theta array
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\n".format(rlen, tlen)
        #Fill in string with velocity information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8f}\t{1:.8f}\t{2:.8f}\n".format(
                                velarr["r"][ti, ri], velarr["th"][ti, ri],
                                velarr["phi"][ti, ri])

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_value("inp_path"), "velocity.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #






#









#
