##FILE: Radlite
##PURPOSE: Contains the classes RadliteModel() and RadliteSpectrum().


##Below Section: IMPORT necessary packages
import astropy.constants as const
from classUtils import func_timer
from datetime import datetime as dater
import json
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os.path
from os import listdir as listdirer
from scipy.interpolate import interp1d as interper
from scipy import ndimage as ndimager
import subprocess
import radmc #!!!!!!!!!!!!!!!!!!
try:
    import astropy.io.fits as fitter
except ImportError:
    import pyfits as fitter

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
    #BELOW FROM READ_MOLECULE_LAMBDA
    #c0     = 2.99792458E10
    #h0     = 6.62620755E-27
    #kB0     = 1.380658E-16
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
class RadliteModel():
    @func_timer
    def __init__(self, infilename, hitranfilename):
        """
        Method: __init__
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Initializes an instance of the RadliteModel() class.
            > Checks that input file (infilename) has minimum required
              ...parameters to run the run_radlite() method.
            > Processes all molecular lines that fit within the user's defined
              ...criteria.
            > Splits work for running RADLite for the processed lines across the
              ...given number of cores.
        Inputs: 2 required
            > infilename (required)
                - Type: string
                - Example: "/User/path/to/file/input_radlite.json"
                - Description: A .json file containing ALL input parameters
                  ...for this class instance.
            > hitranfilename (required)
                - Type: string
                - Example: "/User/path/to/file/data_hitran.json"
                - Description: A .json file containing ALL molecular information
                  ...and names of the associated HITRAN file.
        Outputs: N/A
        Notes:
            > Creates an underlying dictionary to hold all input parameters...
              ...and future attributes.
        """
        ##Below Section: READ IN + STORE input files
        #Read in input RADLite data
        with open(os.path.join(infilename)) as openfile:
            inputdict = json.load(openfile)
        #Store in secret dictionary (stripping out comments)
        self._attrdict = {}
        for key in inputdict:
            self._attrdict[key] = inputdict[key]["value"]
        self._attrdict["units"] = {} #Inner dictionary to hold units
        #NOTE: Unit dictionary for calculated quantities only! NOT inputs

        #Read in HITRAN data and then extract desired molecule
        with open(hitranfilename) as openfile:
            hitrandict = json.load(openfile)
        #Store molecular data for specified molecule
        moldict = hitrandict[self._attrdict["molname"]] #Extract mol. data
        for key in moldict: #Store molecule-specific data
            self._attrdict[key] = moldict[key]
        self._attrdict["molmass"] = moldict["molweight"]*amu0 #Mol. mass


        ##Below Section: PRINT a welcome message, if so desired
        if self.get_attr("verbose"):
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
            print("")
            print("-"*10+"\n"+"Now preparing all RADLite input files...")
            print("")


        ##Below Section: CHECK inputs for user error
        self._check_inputs()


        ###TEMP. BLOCK START
        ##Below Section: TEMPORARY - USE PYRADMC TO READ IN RADMC MODEL
        modradmc = radmc.radmc_model("./")
        self._set_attr(attrname="radius", attrval=modradmc.radius,
                        attrunit="cm")
        self._set_attr(attrname="theta", attrval=modradmc.theta, attrunit="rad")
        self._set_attr(attrname="dustdensity",
                            attrval=modradmc.read_dustdensity()[:,:,0].T,
                            attrunit=r"cm$^{-3}$")
        self._set_attr(attrname="dusttemperature",
                            attrval=modradmc.dusttemperature[:,:,0].T,
                            attrunit="K")
        self._set_attr(attrname="gastodust",
                        attrval=float(modradmc.pars['gastodust'].split(";")[0]))
        #NOTE: !!! Assuming one species read in; hence the [:,:,0] for duststuff
        ###TEMP. BLOCK END


        ##Below Section: WRITE RADLite input files
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Writing RADLite input files...")
            print("")
        self._write_abundanceinp() #Abundance
        self._write_densityinp() #Gas density
        self._write_gastemperatureinp() #Gas temperature
        self._write_radliteinp() #Radlite input
        self._write_turbulenceinp() #Turbulence
        self._write_velocityinp() #Velocity


        ##Below Section: PREPARE all molecular line and population data
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("-"*10+"\n"+"Assembling molecular line data...")
            print("")
        self._prep_mol_forall() #Assemble all molecular line data


        ##Below Section: PREPARE molecular and population data per core
        numcores = self.get_attr("numcores")
        #Initialize private structure to hold level population info per core
        #self._set_attr(attrname="_levelpopdicts", attrval={})

        #Prepare pool of cores
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("-"*10+"\n"+"Preparing molecular line data per core...")
            print("")
        with mp.Pool(numcores) as ppool:
            lpopdicts = ppool.map(self._prep_mol_forcore,
                                    range(0, numcores))
        self._set_attr(attrname="_levelpopdicts", attrval=lpopdicts)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done preparing molecular line data per core!")
            print("")


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("-"*10+"\n"+"RADLite preparation is complete!")
            print("You can now run RADLite using the run_radlite() method.")
            print("-"*20)
            print("")
        return
    #



    ##STORE AND FETCH METHODS
    def get_attr(self, attrname):
        """
        Method: get_attr
        Purpose:
            > Fetches the value of the given attribute.
        Inputs: 1 required
            > attrname (required)
                - Type: string
                - Example: "mstar"
                - Description: Name of the given attribute (such as "mstar" for
                  ...the stellar mass).
        Outputs: 1
            > <Attribute value>
                - Type: Varies
                - Example: 1.0E33
                - Description: The value of the given attribute.
        Notes:
            > If an invalid name is given, then the code will return an error
              ...listing all available attributes.
        """
        ##Below Section: DETERMINE requested attribute
        #Try accessing it
        try:
            return self._attrdict[attrname]
        except KeyError: #If attribute not yet recorded...
            pass

        #If that doesn't work, try calculating it
        try:
            eval("self._calc_"+attrname+"()") #Try calculating it
            return self._attrdict[attrname]
        except AttributeError:
            pass
        except KeyError: #If attribute not calculable...
            pass

        #If that doesn't work, try reading it in
        try:
            eval("self._read_"+attrname+"()") #Try reading it in
            return self._attrdict[attrname]
        except AttributeError:
            pass
        except KeyError: #If attribute not readable
            pass

        #Otherwise, raise an error
        raise AttributeError("'"+attrname+"' doesn't seem to be a valid "
                            +"attribute.  Valid attributes are "
                            +"(in alphabetical order):\n"
                            +str(np.sort([key for key in self._attrdict]))+".\n"
                            +"Run the method run_radlite() (if you haven't "
                            +"yet) to automatically populate more "
                            +"attributes.\n"
                            +"Alternatively, you can pass the name of a "
                            +"supported physics component (e.g., 'velocity_phi'"
                            +") to the get_attr() method to populate that "
                            +"component on its own.")
    #


    def get_unit(self, attrname):
        """
        Method: get_unit
        Purpose:
            > Fetches the automatic unit of the given attribute.
        Inputs: 1 required
            > attrname (required)
                - Type: string
                - Example: "mstar"
                - Description: Name of the given attribute (such as "mstar"
                  ...for the stellar mass).
        Outputs: 1
            > <Attribute value>
                - Type: string
                - Example: "g"
                - Description: The automatic unit of the given attribute.
        Notes:
            > If the given attribute is unitless, then the code will return
              ...an empty string.
            > If an invalid name is given, then the code will return an error
              ...listing all available attributes.
        """
        ##Below Section: RETURN unit of attribute under given name + EXIT
        try:
            return self._attrdict["units"][attrname]
        except KeyError:
            raise KeyError("No unit returned because this attribute is either "
                            +"invalid or has not been populated yet.  If the "
                            +"latter is the case, then you can run "
                            +"the method run_radlite() (if you haven't "
                            +"yet) to automatically populate more "
                            +"attributes.\n"
                            +"Alternatively, you can pass the name of a "
                            +"supported physics component (e.g., 'velocity'"
                            +") to the get_attr() method to populate that "
                            +"component, and therefore its unit, on its own.")
    #


    def _get_core_attr(self, attrname, dictname, pind):
        """
        Method: _get_core_attr
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Fetches the portion of the given attribute that has been
              ...partitioned for core number #pind.
        Inputs: 3 required
            > attrname (required)
                - Type: string
                - Example: "Eup"
                - Description: Name of a desired attribute (such as "Eup" for
                  ...upper energy levels).
            > dictname (required)
                - Type: string
                - Example: "_hitrandict"
                - Description: Name of the underlying dictionary containing
                  ...the desired attribute (such as "_hitrandict" for the
                  ...underlying HITRAN database).
            > pind (required)
                - Type: int
                - Example: 2
                - Description: Index (starting from 0) of the core.
        Outputs: 1
            > <Attribute value>
                - Type: Varies
                - Example: <Array of upper energies>
                - Description: Value of the given attribute that has been
                  ...partitioned for this core.
        Notes:
            > Only certain attributes can/should be accessed with this method.
        """
        boundshere = self.get_attr("_splitinds")[pind]
        return self.get_attr(dictname)[attrname][boundshere[0]:boundshere[1]]
    #


    def _set_attr(self, attrname, attrval, attrunit=""):
        """
        Method: _set_attr
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Sets the given value for the given attribute in the underlying
              ...database.
        Inputs: 2 required, 1 optional
            > attrname (required)
                - Type: string
                - Example: "mstar"
                - Description: Name of a desired attribute (such as "mstar" for
                  ...the stellar mass).
            > attrval (required)
                - Type: Varies
                - Example: 1.0E33
                - Description: Desired value for the given attribute.
            > attrunit (optional)
                - Type: string
                - Example: 2
                - Description: Index (starting from 0) of the core.
        Outputs: N/A
        Notes:
            > Only certain attributes can/should be accessed with this method.
        """
        ##Below Section: RECORD given attribute under given name + EXIT
        self._attrdict[attrname] = attrval
        self._attrdict["units"][attrname] = attrunit #Assign unit
        return
    #



    ##CHECK METHODS
    def _check_inputs(self):
        ##Below Section: CHECK user inputs; make sure they are valid
        #Make sure that desired image cube output is valid
        validimage = [0, 2]
        if self.get_attr("image") not in validimage:
            raise ValueError("Sorry, the image you have chosen ("
                    +str(inpdict["image"])+") is not a valid image.  The "
                    +"value of image must be one of the following: "
                    +str(validimage)+".")

        #Truncate number of cores, if too many requested
        maxcores = mp.cpu_count()
        if self.get_attr("numcores") > (maxcores - 1):
            newnumcores = max([(maxcores - 1), 1])
            print("Whoa, looks like you requested too many cores (you "
                    +"requested "+str(self.get_attr("numcores"))+")!")
            print("Looks like you only have "+str(maxcores)+" available in "
                    +"total, so we're reducing the number of cores down to "
                    +str(newnumcores)+".")
            self._set_attr(attrnum="numcores", attrval=newnumcores)

        #Make sure that desired LTE format is valid
        if not isinstance(self.get_attr("lte"), bool):
            raise ValueError("Sorry, the lte you have chosen ("
                    +str(self.get_attr("lte"))+") is not a valid lte.  "
                    +"The value of lte must be a boolean (True or False).")


        ##Below Section: CHOOSE RADLite exec. based on desired output
        #Prepare executable for desired image cube output
        if self.get_attr("image") == 0: #If desired output is spectrum
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Will prepare a spectrum-formatted image cube...")
                print("")
            self._set_attr(attrname="executable",
                            attrval=self.get_attr("exe_path")+"RADlite")
        elif self.get_attr("image") == 2: #If desired output is circ.
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Will prepare a circular-formatted image cube...")
                print("")
    #



    ##RUN METHODS
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
        #initfilelist = ["problem_params.pro", "line_params.ini"]
        cpfilelist = ["radius.inp",
                        "theta.inp", "frequency.inp", "density.inp",
                        "dustdens.inp", "dusttemp_final.dat",
                        "dustopac.inp", "dustopac_*.inp", "abundance.inp",
                        "temperature.inp", "velocity.inp"]


        ##Below Section: SET UP result directory
        rundir = self.get_attr("run_dir")
        if self.get_attr("dodate"): #If date should be appended to dir. name
            currtime = dater.now()
            timestamp = "_{0}_{1}h{2}m{3}s".format(str(currtime.date()),
                            str(currtime.hour), str(currtime.minute),
                            str(currtime.second))
            rundir += timestamp

        #Make the desired directory, if nonexistent
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Attempting to create run directory "+rundir+"...")
        comm = subprocess.call(["mkdir", rundir]) #Create directory
        if (comm != 0): #Will override files within, otherwise
            if self.get_attr("verbose"): #Verbal output, if so desired
                print(rundir+" already exists. Will overwrite files within.")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("All input files and final data will be saved to the "
                    +"following folder in the current directory: "+rundir)

        #Copy initial files into final directory
        #if self.get_attr("verbose"): #Verbal output, if so desired
        #    print("Copying over initial data files into "+rundir+"...")
        #for iname in initfilelist:
        #    comm = subprocess.call(["cp",
        #                    os.path.join(self.get_attr("inp_path"), iname),
        #                    rundir+"/"])
        #if self.get_attr("verbose"): #Verbal output, if so desired
        #    print("Done copying over initial data files!\n")


        ##Below Section: PROCESS data for LTE/NLTE and mol. line treatment
        #Read in either LTE or NLTE data
        if self.get_attr("lte"): #If LTE treatment desired
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Extracting LTE molecular data...")
                print("")
            self._read_hitran()
        else: #Else, if non-LTE treatment desired
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Extracting NLTE molecular data...")
                print("")
            molfile = "molfile.dat"
            self._read_lambda()


        ##Below Section: RUN RADLITE on single/multiple cores in subfolders
        numcores = self.get_attr("numcores")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Running RADLite on "+str(numcores)+" core(s)...")

        #Prepare working directory for cores
        workingdir = os.path.join('./', rundir, "workingdir")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Generating working directory for cores called: "+workingdir)
        comm = subprocess.call(["mkdir", workingdir]) #Create subdirectory
        if (comm != 0): #If directory already exists, replace it
            comm = subprocess.call(["rm", "-r", workingdir]) #Erase prev. dir.
            comm = subprocess.call(["mkdir", workingdir]) #New empty dir.

        #Prepare final directory for core output
        outputdir = os.path.join('./', rundir, "outputdir")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("All core outputs will be stored in: "+outputdir)
        comm = subprocess.call(["mkdir", outputdir]) #Create subdirectory
        if (comm != 0): #If directory already exists, replace it
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Replacing previous "+outputdir+"...")
            comm = subprocess.call(["rm", "-r", outputdir]) #Erase prev. dir.
            comm = subprocess.call(["mkdir", outputdir]) #New empty dir.

        #Prepare pool of cores
        plist = []
        for ai in range(0, numcores):
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Prepping "+str(ai)+"th core...")

            #Make a core-specific directory for this run
            cpudir = os.path.join(workingdir, ("workingdir_cpu"+str(ai)+"/"))
            comm = subprocess.call(["mkdir", cpudir]) #Create subdirectory
            #Copy initial files into this new core subdirectory
            #for iname in initfilelist:
            #    comm = subprocess.call(["cp", rundir+"/"+iname, cpudir])

            #Call core routine
            phere = mp.Process(target=self._run_core,
                                    args=(ai, cpudir, outputdir,))
            plist.append(phere)
            #Start this core
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Starting "+str(ai)+"th core in "+cpudir+"...")
            phere.start()

        #Close pool of cores
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done running core(s)!")
        for ai in range(0, numcores):
            plist[ai].join()
            if self.get_attr("verbose"): #Verbal output, if so desired
                print(str(ai)+"th core finished and closed!")

        #Delete working directory for cores
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Deleting working directory used for cores...")
        comm = subprocess.call(["rm", "-r", workingdir]) #Erase working dir.


        ##Below Section: FINISH and EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done running run_radlite()!")
        return
    #


    @func_timer
    def _run_core(self, pind, cpudir, outputdir):
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
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(pind)+"th core has started working...")
        ##Below Section: COPY OVER radlite physics/structure input files
        copyfiles = ["abundance.inp", "density.inp", "dustdens.inp", "dustopac.inp", "dustopac_1.inp", "dusttemp.info", "dusttemp_final.dat", "frequency.inp", "line.inp", "radius.inp", "radlite.inp", "starspectrum.inp", "scatsource.dat", "theta.inp", "temperature.inp", "turbulence.inp", "velocity.inp"]
        for filehere in copyfiles:
            comm = subprocess.call(["cp",
                            os.path.join(self.get_attr("inp_path"), filehere),
                            cpudir+"/"]) #Input file needed for RADLite


        ##Below Section: GENERATE radlite molecular line input files
        self._write_core_moldatadat(cpudir=cpudir, pind=pind) #Mol. data file
        self._write_core_levelpopinp(cpudir=cpudir, pind=pind) #Level pop. file
        self._write_core_linespectruminp(cpudir=cpudir, pind=pind) #Spec. file


        ##Below Section: RUN RADLITE
        with open(cpudir+"RADLITE_core.log", 'w') as openlogfile:
            comm = subprocess.call([self.get_attr("executable")], #Call RADLite
                                cwd=cpudir+"/", #Call within core subdir.
                                stdout=openlogfile) #Send output to log file


        ##Below Section: EXTRACT output files
        #Move log file produced by this core
        comm = subprocess.call(["mv",
                            cpudir+"/RADLITE_core.log",
                            outputdir+"/RADLITE_core_"+str(pind)+".log"])
        #Move molecular data used by this core
        comm = subprocess.call(["mv",
                            cpudir+"/moldata.dat",
                            outputdir+"/moldata_"+str(pind)+".dat"])
        #Move line spectrum output produced by this core
        comm = subprocess.call(["mv",
                            cpudir+"/linespectrum_moldata.dat",
                            outputdir+"/linespectrum_moldata_"+str(pind)+".dat"])
        #Move level population output produced by this core
        comm = subprocess.call(["mv",
                            cpudir+"/levelpop_moldata.dat",
                            outputdir+"/levelpop_moldata_"+str(pind)+".dat"])
        if self.get_attr("image") == 2:
            #Move circular 3D image cube output produced by this core
            comm = subprocess.call(["mv",
                            cpudir+"/lineposvelcirc_moldata.dat",
                            outputdir+"/lineposvelcirc_moldata_"+str(pind)+".dat"])


        ##Below Section: EXIT
        #comm = subprocess.call(["rm", "-r", cpudir]) #Erase core dir.
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(pind)+"th core has finished working!")
        return
    #



    ##OUTPUT DISPLAY METHODS
    def plot_attr(self, yattrname, xattrname=None, fig=None, figsize=(10,10),
        s=30, linewidth=3, linestyle="-", marker="o", color="black",
        xlog=False, ylog=False, xscaler=1.0, yscaler=1.0, alpha=1.0,
        xlim=None, ylim=None,
        xunit=None, yunit=None, cbarunit=None,
        xlabel=None, ylabel=None, cbarlabel=None,
        cbarrotation=270, cbarlabelpad=25,
        axisfontsize=16, titlefontsize=18, legfontsize=16,
        tickfontsize=14, title="",
        dolegend=False, leglabel="", legloc="best",
        dopart=False, dosave=False, savename="testing.png"):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: INITIALIZE empty plot, if no existing plot given
        if fig is None:
            fig = plt.figure(figsize=figsize)


        ##Below Section: FETCH x and y-axis values
        #Fetch y-axis values
        yvals = self.get_attr(yattrname)
        #Generate numerical x-axis if none given
        if xattrname is None:
            xattrname = ""
            xvals = np.arange(0, len(yvals))
            xunit = ""
        else:
            xvals = self.get_attr(xattrname)


        ##Below Section: PLOT as either 2D gradient or 1D line+scatter
        ndim = len(np.asarray(yvals).shape) #Number of dimensions for plot
        xshape = np.asarray(xvals).shape #Shape of x-array
        if ndim == 2: #If 2D quantity (assumed axes are radius vs. theta)
            #Set desired 2D quantity to new variable z
            zattrname = yattrname
            zvals = self.get_attr(zattrname) #2D quantity
            zunit = self.get_unit(zattrname) #Unit for 2D quantity
            zscaler = yscaler
            #Extract radius and theta for x and y-axes
            yattrname = "theta"
            yvals = self.get_attr(yattrname) #y-axis values
            yunit = self.get_unit(yattrname) #Unit for y-axis
            xattrname = "radius"
            xvals = self.get_attr(xattrname) #x-axis values
            xunit = self.get_unit(xattrname) #Unit for x-axis
            #Plot gradient
            grad = plt.imshow(zvals*zscaler, cmap=color)
        elif (ndim == 1) and (len(xshape) == 1): #If 1D quantities
            #Plot line plot
            plt.plot(xvals*xscaler, yvals*yscaler, color=color,
                    linewidth=linewidth, linestyle=linestyle,
                    alpha=alpha, label=leglabel)
            #Plot scatter plot
            plt.scatter(xvals*xscaler, yvals*yscaler, color=color,
                    marker=marker, s=s, alpha=alpha, label=leglabel)
        else: #If neither 1D nor 2D
            raise ValueError("Oh no!  Make sure you've passed in y and/or x "
                            +"attributes to plot_attr() that have either 1D or "
                            +"2D values.")


        ##Below Section: SCALE plot axes, if so desired
        #Log scale, if so desired
        if xlog: #For x-axis
            plt.xscale("log")
        if ylog: #For y-axis
            plt.yscale("log")

        #Axis limits, if so desired
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)


        ##Below Section: GENERATE colorbar, if this is a 2D plot
        if ndim == 2:
            if cbarlabel is None:
                #Determine the unit, if not given
                if cbarunit is None:
                    cbarunit = self.get_unit(zattrname) #Automatic unit
                cbarlabel = zattrname.capitalize()
                if cbarunit != "": #Tack on unit, if exists
                    cbarlabel = cbarlabel + " ["+cbarunit+"]"
            cbar = plt.colorbar(grad)
            cbar.set_label(label=cbarlabel, fontsize=titlefontsize,
                            rotation=cbarrotation, labelpad=cbarlabelpad)
            cbar.ax.tick_params(labelsize=tickfontsize)


        ##Below Section: LABEL plot axes, if so desired
        #x-axis labels (with units), if so desired
        if xlabel is None:
            #Set the x-axis label
            xlabel = xattrname.capitalize()
        #Determine the unit, if not given
        if xunit is None:
            xunit = self.get_unit(xattrname) #Automatic unit
        if xunit != "": #Tack on unit, if exists
            xlabel = xlabel +" ["+xunit+"]"
        plt.xlabel(xlabel, fontsize=axisfontsize)

        #y-axis labels (with units), if so desired
        if ylabel is None:
            #Set the y-axis label
            ylabel = yattrname.capitalize()
        #Determine the unit, if not given
        if yunit is None: #For y-axis
            yunit = self.get_unit(yattrname) #Automatic unit
        if yunit != "": #Tack on unit, if exists
            ylabel = ylabel +" ["+yunit+"]"
        plt.ylabel(ylabel, fontsize=axisfontsize)


        ##Below Section: SET title + legend + tick label size
        #For legend, if so desired
        if dolegend:
            plt.legend(loc=legloc, frameon=False, fontsize=legfontsize)

        #For title
        plt.title(title, fontsize=titlefontsize)

        #Set font size of tick labels
        plt.xticks(fontsize=tickfontsize)
        plt.yticks(fontsize=tickfontsize)


        ##Below Section: FINALIZE + EXIT
        #Stop here, if these commands are part of an external plot
        if dopart:
            return

        #Otherwise, save the figure, if so desired
        if dosave:
            plt.savefig(savename)
            plt.close()
            return

        #Otherwise, display the figure
        plt.show()
        return
    #



    ##PREPARATION METHODS
    def _prep_mol_forcore(self, pind):
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
        dictname = "_hitrandict"
        Euparr = self._get_core_attr("Eup", dictname=dictname, pind=pind)
        Elowarr = self._get_core_attr("Elow", dictname=dictname, pind=pind)
        guparr = self._get_core_attr("gup", dictname=dictname, pind=pind)
        glowarr = self._get_core_attr("glow", dictname=dictname, pind=pind)
        wavenumarr = self._get_core_attr("wavenum",
                                                dictname=dictname, pind=pind)
        vuparr = self._get_core_attr("vup", dictname=dictname, pind=pind)
        vlowarr = self._get_core_attr("vlow", dictname=dictname, pind=pind)
        quparr = self._get_core_attr("qup", dictname=dictname, pind=pind)
        qlowarr = self._get_core_attr("qlow", dictname=dictname, pind=pind)
        #
        tlen = len(self.get_attr("theta"))//2
        tempgasarr = self.get_attr("gastemperature")[0:tlen,:]
        psum = self.get_attr("psum")
        psumtemp = self.get_attr("psum_temp")


        ##Below Section: COMBINE levels + REMOVE duplicates to get unique levels
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Counting up unique levels for "+str(pind)+"th core...")
        #Combine energies, transitions, and degeneracies
        Eallarr = np.concatenate((Elowarr, Euparr))
        vallarr = np.concatenate((vlowarr, vuparr))
        qallarr = np.concatenate((qlowarr, quparr))
        gallarr = np.concatenate((glowarr, guparr))
        #Sort combined lists by energies
        sortinds = np.argsort(Eallarr)
        Eallarr = Eallarr[sortinds]
        vallarr = vallarr[sortinds]
        qallarr = qallarr[sortinds]
        gallarr = gallarr[sortinds]

        #Extract indices for unique levels only
        keepinds = [True]*len(Eallarr)
        for ai in range(1, len(Eallarr)):
            if ( #For non-unique E, v, and g value combinations
                        ((np.abs(Eallarr[ai-1]-Eallarr[ai])
                                /1.0/(Eallarr[ai]+0.1)) < 0.0001)
                        and (vallarr[ai-1] == vallarr[ai])
                        and (gallarr[ai-1] == gallarr[ai])):
                keepinds[ai] = False
        uniqinds = np.unique(np.where(keepinds)[0])

        #Combine and apply indices
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Removing any duplicate levels...  "
                    +"To start, there are "+str(len(Eallarr))+" levels.")
        Euniqarr = Eallarr[uniqinds]
        vuniqarr = vallarr[uniqinds]
        quniqarr = qallarr[uniqinds]
        guniqarr = gallarr[uniqinds]
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(len(uniqinds))+" unique indices found, "
                        +"meaning "+str(len(Euniqarr))+" unique levels.")


        ##Below Section: CALCULATE partition sum and level populations
        numlevels = len(Euniqarr)
        numtrans = len(wavenumarr)
        Euniqarr_K = Euniqarr*h0*c0/1.0/kB0 #[K]; Upper energy levels
        #Interpolate partition sum
        psumfunc = interper(psumtemp, psum, kind=self.get_attr("interpolation"))
        psuminterped = psumfunc(tempgasarr)

        #Calculate population levels
        npoparr = np.array([(guniqarr[ai]
                                    *np.exp(-1.0*Euniqarr_K[ai]/tempgasarr))
                            for ai in range(0, numlevels)]) /1.0/psuminterped
        #Trim any ridiculously-low values
        npoparr[npoparr < 1E-99] = 0.0
        #print(np.sum(npoparr))
        #print(npoparr.shape)
        #print(np.sum(tempgasarr))
        #print(np.sum(Euniqarr_K))
        #print("!!!")


        ##Below Section: STORE + RETURN unique level results + EXIT
        lpopdict = {"npop":npoparr, "Euniq":Euniqarr, "guniq":guniqarr,
                            "quniq":quniqarr, "vuniq":vuniqarr,
                            "numlevels":numlevels, "numtrans":numtrans}
        return lpopdict
    #


    @func_timer
    def _prep_mol_forall(self):
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
        ##Below Section: CALL functions to read molecular line data
        self._read_hitran()
        self._read_psum()


        ##Below Section: SPLIT data across given number of cores
        numlines = self.get_attr("numlines") #Number of lines
        numcores = self.get_attr("numcores") #Number of cores
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Dividing up the lines for "
                        +str(numcores)+" cores...")
        #Determine indices for splitting up the data
        tempsplits = np.array([numlines // numcores]*numcores) #Split per core
        remainder = numlines % numcores #The remainder
        tempsplits[0:remainder] += 1 #Spread remainder over cores
        tempcumusplits = np.concatenate((np.array([0]),
                                np.cumsum(tempsplits))) #Cumulate over splits
        splitinds = [[tempcumusplits[ai], tempcumusplits[ai+1]]
                                for ai in range(0, numcores)] #Split per core
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Here are the chosen line intervals per core:")
            print([("Core "+str(ehere[0])+": Interval "+str(ehere[1]))
                                    for ehere in enumerate(splitinds)])
        self._set_attr(attrname="_splitinds", attrval=splitinds) #Split lines


        ##Below Section: EXIT
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
        ##Below Section: EXTRACT necessary parameters
        isonum = self.get_attr("isonum")
        wavenumrange = [1.0E4/self.get_attr("max_mu"),
                                            1.0E4/self.get_attr("min_mu")]
        insmin = self.get_attr("min_ins")
        Eupmax = self.get_attr("max_Eup")
        vupmax = self.get_attr("max_vup")


        ##Below Section: PROCESS all lines in HITRAN file
        filepathandname = os.path.join(
                                    self.get_attr("hit_path"),
                                    self.get_attr("hitran_file"))
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Extracting molecular lines from file "+filepathandname+"...")
        #Read in all filelines
        with open(filepathandname, 'r') as openfile:
            alllines = openfile.readlines()
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("There are "+str(len(alllines))+" molecular lines in total.")

        #Set HITRAN entries and number of characters per HITRAN file entry
        hitranchars = [3, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15, 15, 15,
                            6, 12, 1, 7, 7] #Characters per HITRAN entry
        hitrannames = ["isonum", "wavenum", "ins", "A", "air", "self",
                            "Elow", "temp", "press", "vup", "vlow",
                            "qup", "qlow", "err", "ref", "flag", "gup", "glow"]
        hitrandtypes = [int]*1 +[float]*8 +[str]*7 +[float]*2

        #Record *all* molecular line into a convenient dictionary
        #NOTE: Any undesired lines will be removed later
        hitrandict = {} #Initialize dictionary to hold line attributes
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
        #And also convert v levels into ints (removing any whitespaces)
        vupnum = np.array(
                            [int(hitrandict["vup"][ai].replace(" ", ""))
                                for ai in range(0, len(hitrandict["vup"]))])

        #Throw error if any signs of incomplete data
        if 0 in hitrandict["gup"]:
            raise ValueError("Something's wrong!  Incomplete molecular "
                        +"HITRAN data at wavenum="
                        +str(hitrandict["wavenum"][
                                                hitrandict["gup"]==0])+".")


        ##Below Section: REMOVE any lines outside of desired range
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Removing molecular lines outside of specified criteria...")
        #Extract indices of lines that fall within criteria
        keepinds = np.where( ~(
                    #If incorrect isotopologue
                    (hitrandict["isonum"] != isonum)
                    #If wavenumber outside of desired wavenumber range
                    | (hitrandict["wavenum"] < wavenumrange[0])
                    | (hitrandict["wavenum"] > wavenumrange[1])
                    #If intensity, level, or upper energy beyond given cutoffs
                    | (hitrandict["ins"] < insmin)
                    | (hitrandict["Eup"] > Eupmax)
                    | (vupnum > vupmax)))[0]
        numlines = len(keepinds) #Number of lines that fall within criteria

        #Keep only those lines that fall within criteria; delete all other lines
        for key in hitrandict:
            hitrandict[key] = hitrandict[key][keepinds]
        self._set_attr(attrname="_hitrandict", attrval=hitrandict) #Final set
        self._set_attr(attrname="numlines", attrval=numlines) #Final count
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("There are "+str(numlines)+" molecular lines "
                    +"that fall within specified criteria.")


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done extracting molecular data!\n")
        return
    #


    def _read_mstar(self):
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
        ##Below Section: CALL general stellar info reader
        self._read_starinfo()
        return
    #


    #@func_timer
    def _read_psum(self):
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
        ##Below Section: PROCESS all filelines in partition sum file
        #Read in all filelines
        with open(self.get_attr("psumfile"), 'r') as openfile:
            alllines = openfile.readlines()

        #Record information for all chemical species from filelines
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Extracting partition function info from "
                                        +self.get_attr("psumfile")+"...")
        specieslist = alllines[1].split() #List of all chemical species' names
        numspecies = len(specieslist) #Number of chemical species
        molmassarr = np.asarray(alllines[2].split()
                                        ).astype(float) #List of all mol. mass
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("There are a total of "+str(numspecies)+" chemical species "
                                        +"in "+self.get_attr("psumfile")+".")

        #Extract info for requested molecule (chemical species)
        try: #Try to determine index of requested molecule
            iind = specieslist.index(self.get_attr("molname")) #Mol's index
        except ValueError: #If requested molecule not present in list
            raise ValueError("Oh no! No partition sum table found in "
                            +self.get_attr("psumfile")+" for the requested "
                            +"molecule "+self.get_attr("molname")+".")
        psumtemparr = np.array([float(alllines[ai].split()[0])
                                    for ai in range(3, len(alllines))]) #Temp.
        psumarr = np.array([float(alllines[ai].split()[iind+1])
                                    for ai in range(3, len(alllines))]) #Part.


        ##Below Section: STORE partition info + EXIT
        self._set_attr(attrname="psum", attrval=psumarr)
        self._set_attr(attrname="psum_temp", attrval=psumtemparr,
                        attrunit="K")
        self._set_attr(attrname="psum_molmass", attrval=molmassarr[iind],
                        attrunit="g")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Partition information for "+self.get_attr("molname")+" "
                                        +"has been successfully extracted!\n")
        return
    #


    def _read_rstar(self):
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
        ##Below Section: CALL general stellar info reader
        self._read_starinfo()
        return
    #


    #@func_timer
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
                                    self.get_attr("inp_path"), "starinfo.inp")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Extracting stellar info from file "+filepathandname+"...")
        #Read in all filelines
        with open(filepathandname, 'r') as openfile:
            alllines = openfile.readlines()


        ##Below Section: STORE read-in data + EXIT
        rstar = float(alllines[1]) #Stellar radius
        mstar = float(alllines[2]) #Stellar mass
        teff = float(alllines[3]) #Stellar eff. temperature
        starinfodict = {"mstar":mstar, "rstar":rstar, "teff":teff}
        #Store together and individually
        self._set_attr(attrname="starinfo", attrval=starinfodict)
        self._set_attr(attrname="mstar", attrval=mstar, attrunit="g")
        self._set_attr(attrname="rstar", attrval=rstar, attrunit="cm")
        self._set_attr(attrname="teff", attrval=teff, attrunit="K")
        return
    #


    def _read_teff(self):
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
        ##Below Section: CALL general stellar info reader
        self._read_starinfo()
        return
    #



    ##CALCULATION METHODS
    ###NOTE: MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING METHODS
    #@func_timer
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
        temparr = self.get_attr("dusttemperature")
        rlen = len(self.get_attr("radius"))
        tlen = len(self.get_attr("theta"))
        gamval = self.get_attr("gamma")
        abundrange = [self.get_attr("min_abun"), self.get_attr("max_abun")]
        temp_fr = self.get_attr("temp_fr")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Calculating abundance...")


        ##Below Section: CALCULATE abundance based on specified abundance mode
        if temp_fr is False: #For constant abundance
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Setting a constant abundance...")
            abundarr = np.ones(shape=(tlen, rlen))*abundrange[1]
        else: #For abundance that changes below freeze-out temperature
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Setting abundance to min. below {0:f}K..."
                                .format(temp_fr))
            abundarr = np.ones(shape=(tlen, rlen))*abundrange[1]
            abundarr[temparr < temp_fr] = abundrange[0]
        self._set_attr(attrname="collpartner", attrval=0.0)


        ##NOTE: IN MAKE_ABUNDANCE.PRO, THERE WAS PART HERE ABOUT MOL_DESTRUCT

        ##Below Section: RECORD calculated abundance + EXIT
        self._set_attr(attrname="abundance", attrval=abundarr)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating abundance!\n")
        return
    #


    ###NOTE: _CALC_GASDENSITY WILL BE MOVED TO BASEMODEL CLASS; MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING METHODS
    #@func_timer
    def _calc_gasdensity(self):
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
        ##Below Section: CALCULATE gas density based on dust density
        gastodustratio = self.get_attr("gastodust") #Gas to dust ratio
        gasdensarr = self.get_attr("dustdensity").copy()*gastodustratio
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Calculating gas density from dust density...")


        ##Below Section: RECORD calculated gas temperature + EXIT
        self._set_attr(attrname="gasdensity", attrval=gasdensarr,
                        attrunit=r"cm$^{-3}$")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating gas density!\n")
        return
    #


    ###NOTE: _CALC_GASTEMPERATURE WILL BE MOVED TO BASEMODEL CLASS; MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING METHODS
    #@func_timer
    def _calc_gastemperature(self):
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
        gastemparr = self.get_attr("dusttemperature").copy()
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Calculating gas temperature...")
            print("Assuming gas temperature = dust temperature...")


        ##Below Section: RECORD calculated gas temperature + EXIT
        self._set_attr(attrname="gastemperature", attrval=gastemparr,
                        attrunit="K")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating gas temperature!\n")
        return
    #


    ###NOTE: _CALC_TURBULENCE WILL BE MOVED TO BASEMODEL CLASS; MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING METHODS
    #@func_timer
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
        temparr = self.get_attr("dusttemperature")
        gamval = self.get_attr("gamma")
        muval = self.get_attr("mu")
        weightval = self.get_attr("molweight")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Calculating turbulence...")
            print("Using first dust component temperature to determine "
                    +"turbulent velocities...")


        ##Below Section: CALCULATE turbulence (alpha-model) using sound speed
        csarr = np.sqrt(gamval*kB0*temparr/1.0/(muval*mp0)) #Sound speed
        turbarr = self.get_attr("alpha")*csarr #Turbulence
        #Add thermal broadening in quadrature
        thermbroadarr = np.sqrt(2.0*kB0*temparr/(weightval*mp0)) #Therm. broad.
        turbarr = np.sqrt((turbarr**2) + (thermbroadarr**2)) #Updated turbulence


        ##Below Section: RECORD calculated turbulence + EXIT
        self._set_attr(attrname="turbulence", attrval=turbarr,
                        attrunit=r"cm s$^{-1}$")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating turbulence!\n")
        return
    #


    ###NOTE: _CALC_VELOCITY WILL BE MOVED TO BASEMODEL CLASS; MORE COMPLEX STRUCTURES WILL HAVE OVERRIDING _CALC_VELOCITY METHODS
    #@func_timer
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
        rarr = self.get_attr("radius") #Radii
        rlen = len(rarr) #Number of radius points
        tlen = len(self.get_attr("theta"))//2 #Number of theta points
        rexparr = np.resize(rarr, (tlen, rlen)) #Radii, expended to 2D
        mstar = self.get_attr("starinfo")["mstar"] #Stellar mass
        rstar = self.get_attr("starinfo")["rstar"] #Stellar radius
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Calculating velocity field components...")
            print("Used starinfo.inp file for mstar and rstar...")


        ##Below Section: CALCULATE Keplerian velocity
        vel_radial = np.zeros(shape=(tlen, rlen)) #Radial velocity
        vel_theta = np.zeros(shape=(tlen, rlen)) #Theta velocity
        vel_phi = np.sqrt(G0*mstar/1.0/rexparr) #Phi velocity


        ##Below Section: RECORD velocity components and EXIT
        self._set_attr(attrname="velocity_theta", attrval=vel_theta,
                            attrunit=r"cm s$^{-1}$")
        self._set_attr(attrname="velocity_radial", attrval=vel_radial,
                            attrunit=r"cm s$^{-1}$")
        self._set_attr(attrname="velocity_phi", attrval=vel_phi,
                            attrunit=r"cm s$^{-1}$")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating velocity components!\n")
        return
    #


    def _calc_velocity_radial(self):
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
        ##Below Section: CALL general velocity calculator
        self._calc_velocity()
        return
    #


    def _calc_velocity_theta(self):
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
        ##Below Section: CALL general velocity calculator
        self._calc_velocity()
        return
    #


    def _calc_velocity_phi(self):
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
        ##Below Section: CALL general velocity calculator
        self._calc_velocity()
        return
    #



    ##WRITE METHODS
    #@func_timer
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
        abundarr = self.get_attr("abundance") #Abundance data
        collpartner = self.get_attr("collpartner")
        rlen = len(self.get_attr("radius")) #Length of radius array
        tlen = len(self.get_attr("theta"))//2 #Half-length of theta array
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\n".format(rlen, tlen)
        #Fill in string with abundance information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8e}\t{1:.8f}\n".format(abundarr[ti, ri],
                                                collpartner)

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_attr("inp_path"), "abundance.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    @func_timer
    def _write_core_levelpopinp(self, cpudir, pind):
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
        ##Below Section: EXTRACT level population and partition sum information
        #Extract level population information
        lpopdict = self.get_attr("_levelpopdicts")[pind] #Level population info
        numlevels = lpopdict["numlevels"] #Number of energy levels
        npoparr = lpopdict["npop"] #Level populations
        Euniqarr = lpopdict["Euniq"] #Unique energy levels
        guniqarr = lpopdict["guniq"] #Unique degeneracy levels
        tlen = len(self.get_attr("theta"))//2 #Number of theta points
        rlen = len(self.get_attr("radius")) #Number of radius points


        ##Below Section: BUILD string containing level population information
        #FOR LEVELPOP_MOLDATA
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\t{2:d}\t1\n".format(rlen, tlen, numlevels)

        #Convert energies and degeneracies into strings
        Estr = "\t".join((Euniqarr*c0*h0).astype(str))
        gstr = "\t".join(guniqarr.astype(str))
        #Tack energies and degeneracies onto overall string
        writestr += Estr + "\n"
        writestr += gstr + "\n"

        #Fill in string with level populations
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                npopstrhere = "\n".join(npoparr[:, ti, ri].astype(str))
                writestr += npopstrhere + "\n"

        #Write the results to file
        outfilename = os.path.join(cpudir, "levelpop_moldata.dat")
        with  open(outfilename, 'w') as openfile:
            openfile.write(writestr)

        #FOR LEVELPOP INFO
        #Set up string
        writestr = "" #Initialize string
        writestr += "-3\nlevelpop_moldata.dat\n0"

        #Write the results to file
        outfilename = os.path.join(cpudir, "levelpop.info")
        with  open(outfilename, 'w') as openfile:
            openfile.write(writestr)


        ##Below Section: EXIT function
        return
    #


    #@func_timer
    def _write_core_linespectruminp(self, pind, cpudir):
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
                    self.get_attr("passband"), "Width of line passband [km/s]")
        writestr += "{0:<8.2f}{1:<50s}\n".format(
                    self.get_attr("vsampling"), "Velocity sampling [km/s]")
        writestr += ("-"*65)+"\n"
        writestr += "{0:<8d}Format number\n".format(2)
        writestr += "moldata.dat\tMolecular data file\n"
        writestr += "{0:<8d}{1:<50s}\n".format(self.get_attr("image"),
                                "Command (0=spectrum, 2=image[3-D P/V cube])")
        writestr += "{0:<8.1f}{1:<50s}\n".format(
                    1.0, "Distance in [pc]")
        writestr += "{0:<8.1f}{1:<50s}\n".format(
                    self.get_attr("incl"), "Inclination [deg]")
        writestr += "{0:<8.1f}{1:<60s}\n".format(
                    self.get_attr("vlsr"), "Radial velocity, rel. to local "
                                        +"standard of rest [km/s]")
        writestr += "{0:<8d}{1:<50s}\n".format(
                    len(self._get_core_attr("wavenum",
                                        pind=pind, dictname="_hitrandict")),
                    "Nr of lines to make spectrum/image")
        writestr += "{0:<8d}Starting line to make spectrum/image\n".format(1)
        #
        if self.get_attr("image") == 2: #For 3D image cube
            npix = int(np.ceil(self.get_attr("imwidth")
                                    /1.0/self.get_attr("ssampling")))
            imwidth_cm = self.get_attr("imwidth")/1.0/au0
            writestr += "{0:<8d}{1:<50s}\n".format(npix, "Nr of x pixels")
            writestr += "{0:<8d}{1:<50s}\n".format(npix, "Nr of y pixels")
            writestr += "{0:<8d}image size in cm?\n".format(1)
            writestr += "{0:<8.1e}{1:<50s}\n".format(
                            imwidth_cm, "size x direction")
            writestr += "{0:<8.1e}{1:<50s}\n".format(
                            imwidth_cm, "size y direction")
            writestr += "{0:<8d}Phi offset?\n".format(0)
            writestr += "{0:<8d}x offset?\n".format(0)
            writestr += "{0:<8d}y offset?\n".format(0)
            writestr += "{0:<8d}add star?\n".format(1)


        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(cpudir, "linespectrum.inp")
        with  open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    #@func_timer
    def _write_core_moldatadat(self, cpudir, pind):
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
        #Extract molecular line information
        molname = self.get_attr("molname")
        molweight = self.get_attr("molweight")
        dictname = "_hitrandict"
        Euparr = self._get_core_attr("Eup", pind=pind, dictname=dictname)
        Elowarr = self._get_core_attr("Elow", pind=pind, dictname=dictname)
        guparr = self._get_core_attr("gup", pind=pind, dictname=dictname)
        glowarr = self._get_core_attr("glow", pind=pind, dictname=dictname)
        wavenumarr = self._get_core_attr("wavenum", pind=pind,
                                                dictname=dictname)
        Aarr = self._get_core_attr("A", pind=pind, dictname=dictname)
        vuparr = self._get_core_attr("vup", pind=pind, dictname=dictname)
        vlowarr = self._get_core_attr("vlow", pind=pind, dictname=dictname)
        quparr = self._get_core_attr("qup", pind=pind, dictname=dictname)
        qlowarr = self._get_core_attr("qlow", pind=pind, dictname=dictname)
        #Sort the lists by E_low
        sortinds = np.argsort(Elowarr)
        Euparr = Euparr[sortinds]
        Elowarr = Elowarr[sortinds]
        guparr = guparr[sortinds]
        glowarr = glowarr[sortinds]
        wavenumarr = wavenumarr[sortinds]
        Aarr = Aarr[sortinds]
        vuparr = vuparr[sortinds]
        vlowarr = vlowarr[sortinds]
        quparr = quparr[sortinds]
        qlowarr = qlowarr[sortinds]

        #Extract level population information
        lpopdict = self.get_attr("_levelpopdicts")[pind]
        Euniqarr = lpopdict["Euniq"]
        guniqarr = lpopdict["guniq"]
        vuniqarr = lpopdict["vuniq"]
        quniqarr = lpopdict["quniq"]
        numlevels = lpopdict["numlevels"]
        numtrans = lpopdict["numtrans"]


        ##Below Section: BUILD string to form the molecular data file
        writestr = ""
        writestr += "!MOLECULE\n{0:s}\n".format(molname)
        writestr += "!MOLECULAR WEIGHT\n{0:4.1f}\n".format(molweight)
        writestr += "!NUMBER OF ENERGY LEVELS\n{0:6d}\n".format(numlevels)
        #Tack on unique levels, energies, and degeneracies
        writestr += "!LEVEL + ENERGIES(cm^-1) + WEIGHT + v + Q\n"
        for ai in range(0, numlevels):
            writestr += "{0:5d}{1:12.4f}{2:7.1f}{3:>15s}{4:>15s}\n".format(
                            (ai+1), Euniqarr[ai], guniqarr[ai], vuniqarr[ai],
                            quniqarr[ai])
        #Tack on transitions
        writestr += "!NUMBER OF RADIATIVE TRANSITIONS\n"
        writestr += "{0:6d}\n".format(numtrans)
        writestr += ("!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm^-1) + "
                        +"E_u(cm^-1) + v_l + Q_p + Q_pp\n")
        for ai in range(0, numtrans):
            levu = np.where((np.abs(Euniqarr - Euparr[ai])
                                /1.0/Euparr[ai]) < 1E-4)[0]
            if Elowarr[ai] != 0: #If not down to 0-level
                levl = np.where((np.abs(Euniqarr - Elowarr[ai])
                                /1.0/Elowarr[ai]) < 1E-4)[0]
            else:
                levl = np.where(Euniqarr == 0)[0]
            writestr += ("{0:5d}{1:5d}{2:5d}{3:12.3e}{4:16.7f}".format(
                                (ai+1), levu[0]+1, levl[0]+1, Aarr[ai],
                                wavenumarr[ai])
                        +"{0:12.5f}{1:>15}{2:>15}{3:>15}{4:>15}\n".format(
                                Euparr[ai], vuparr[ai], vlowarr[ai],
                                quparr[ai], qlowarr[ai]))


        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(cpudir, "moldata.dat")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    #@func_timer
    def _write_densityinp(self):
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
        ##Below Section: BUILD string containing gas density information
        #Extract gas density
        gasdensarr = self.get_attr("gasdensity") #Gas density data
        rlen = len(self.get_attr("radius")) #Length of radius array
        tlen = len(self.get_attr("theta"))//2 #Half-length of theta array
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\t1\n".format(rlen, tlen)
        #Fill in string with abundance information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8e}\n".format(gasdensarr[ti, ri])

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_attr("inp_path"), "density.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    #@func_timer
    def _write_gastemperatureinp(self):
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
        ##Below Section: BUILD string containing gas temperature information
        #Extract temperature
        gastemparr = self.get_attr("gastemperature") #Gas temperature data
        rlen = len(self.get_attr("radius")) #Length of radius array
        tlen = len(self.get_attr("theta"))//2 #Half-length of theta array
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\t1\n".format(rlen, tlen)
        #Fill in string with abundance information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8f}\n".format(gastemparr[ti, ri])

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_attr("inp_path"), "temperature.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    #@func_timer
    def _write_radliteinp(self):
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
        ##Below Section: BUILD string to form the input file
        writestr = "" #Initialize string
        writestr += "104\t\tInput format version\n"
        writestr += "="*40+"\n"
        writestr += str("-1")+"\t\tMaximum number; OBSOLETE NITER PARAMETER\n"
        writestr += "2\t\tIteration method\t\t(0=LI,1=ALI,2=ALI+Ng,-1=IS)\n"
        writestr += "0\t\tFlux conservation trick\t\t(0=disable)\n"
        writestr += "1.000e-06\t\tConvergence tolerance\n"
        writestr += "0\t\tConvergence crit type\t\t(0=J_nu,1=j_nu)\n"
        writestr += "1\t\tInitial guess type\t\t(0=triv)\n"
        writestr += "="*40+"\n"
        writestr += "41\t\tNr of mu-angle points\n"
        writestr += "16\t\tNr of phi-angle points\n"
        writestr += "2\t\tType of mu gridding\t\t(0=regular,1=dr/r,2=dr/rspecial)\n"
        writestr += "0.1\t\tLargest allowed error\n"
        writestr += "1.0\t\tDmu wrt dR\t\t(1.0=tangent ray method)\n"
        writestr += "3\t\tExtra mus around mu=0\n"
        writestr += "="*40+"\n"
        if self.get_attr("noisrf"):
            writestr += "0\t\tType of outer boundary\t\t(0=empty,2=microwavebg)\n"
        else:
            writestr += "3\t\tType of outer boundary\t\t(0=empty,2=microwavebg)\n"
        writestr += "2\t\tType of inner boundary\t\t(0=solid,1=empty?,2=central star,3=mirror)\n"
        writestr += "0\t\tType of equator boundary\t\t(0=no,1=disk)\n"
        writestr += "="*40+"\n"
        writestr += "57.3\t\tDefault inclination angle\n"
        writestr += str(self.get_attr("cir_np"))+"\t\tNr Phi-points circular CCD\n"
        writestr += str(self.get_attr("b_per_r"))+"\t\tThe number of b_i of the image per r_i of the grid\n"
        writestr += str(self.get_attr("b_extra"))+"\t\tThe number of extra b_i inside the inner radius\n"
        writestr += "="*40+"\n"
        writestr += "1\t\tSave NLTE\n"
        writestr += "0\t\tSave intensity at inu\n"
        writestr += "-1\t\tSave source at inu\n"
        writestr += "0\t\tSave moments\t\t(-1=save,2=append,3=onlylast)\n"
        writestr += "-1\t\tSave dusttemp/etc\t\t(-1=save,2=append,3=onlylast)\n"
        writestr += "0\t\tSave ALI operator\n"
        writestr += "1\t\tSave phys vars of medium\n"
        writestr += "0\t\tSave flux conservation data\n"
        writestr += "="*40+"\n"
        writestr += "-551\t\tType of setup\n"
        writestr += "0\t\tDo dust?\n"
        writestr += "1\t\tDo lines?\n"
        writestr += "1\t\tInclude dust in lines?\n"
        writestr += "1\t\tInclude star pumping?\n"
        writestr += "="*40+"\n"


        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_attr("inp_path"), "radlite.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    #@func_timer
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
        turbarr = self.get_attr("turbulence") #Turbulence data
        rlen = len(self.get_attr("radius")) #Length of radius array
        tlen = len(self.get_attr("theta"))//2 #Length of theta array
        #Set up string
        writestr = "1\n" #Initialize string
        writestr += "{0:d}\t{1:d}\n".format(rlen, tlen)
        #Fill in string with turbulence information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8f}\n".format(
                                        (turbarr[ti, ri]/1.0E5)) #[cm/s]->[km/s]

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_attr("inp_path"), "turbulence.inp")
        with  open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #


    #@func_timer
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
        vel_radial = self.get_attr("velocity_radial") #Radial velocity data
        vel_theta = self.get_attr("velocity_theta") #Theta velocity data
        vel_phi = self.get_attr("velocity_phi") #Phi velocity data
        rlen = len(self.get_attr("radius")) #Length of radius array
        tlen = len(self.get_attr("theta"))//2 #Length of theta array
        #Set up string
        writestr = "" #Initialize string
        writestr += "{0:d}\t{1:d}\n".format(rlen, tlen)
        #Fill in string with velocity information
        for ri in range(0, rlen):
            for ti in range(0, tlen):
                writestr += "{0:.8f}\t{1:.8f}\t{2:.8f}\n".format(
                                vel_radial[ti, ri],
                                vel_theta[ti, ri],
                                vel_phi[ti, ri])

        ##Below Section: WRITE the results to file + EXIT function
        outfilename = os.path.join(self.get_attr("inp_path"), "velocity.inp")
        with open(outfilename, 'w') as openfile:
            openfile.write(writestr)
        return
    #
#





##
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
        self._attrdict["units"] = {} #Inner dictionary to hold units
        #NOTE: Unit dictionary for calculated quantities only! NOT inputs


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


    def get_unit(self, attrname):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: RETURN unit of attribute under given name + EXIT
        try:
            return self._attrdict["units"][attrname]
        except KeyError:
            raise KeyError("No unit returned because this attribute is either "
                            +"invalid or has not been populated yet.  If the "
                            +"latter is the case, then you can run "
                            +"the method gen_spec() (if you haven't "
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


    def _set_attr(self, attrname, attrval, attrunit=""):
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
        self._attrdict["units"][attrname] = attrunit #Assign unit
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
    def plot_spec(self, yattrname, xattrname=None, fig=None, figsize=(10,10),
        linewidth=3, linestyle="-", color="black",
        xlog=False, ylog=False, xscaler=1.0, yscaler=1.0, alpha=1.0,
        xlim=None, ylim=None, xunit=None, yunit=None,
        xlabel=None, ylabel=None,
        axisfontsize=16, titlefontsize=18, legfontsize=16,
        tickfontsize=14, title="",
        dolegend=False, leglabel="", legloc="best",
        dopart=False, dosave=False, savename="testing.png"):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: INITIALIZE empty plot, if no existing plot given
        if fig is None:
            fig = plt.figure(figsize=figsize)


        ##Below Section: FETCH x and y-axis values
        #Fetch y-axis values
        yvals = self.get_attr(yattrname)
        #Fetch x-axis values
        if xattrname is None:
            xattrname = "wavelength" #Sets x-axis to wavelength, if none given
        xvals = self.get_attr(xattrname)


        ##Below Section: PLOT the desired attributes
        plt.plot(xvals*xscaler, yvals*yscaler, color=color,
                    linewidth=linewidth, linestyle=linestyle,
                    alpha=alpha, label=leglabel)


        ##Below Section: SCALE plot axes, if so desired
        #Log scale, if so desired
        if xlog: #For x-axis
            plt.xscale("log")
        if ylog: #For y-axis
            plt.yscale("log")

        #Axis limits, if so desired
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)


        ##Below Section: LABEL plot axes, if so desired
        #x-axis labels (with units), if so desired
        if xlabel is None:
            #Set the x-axis label
            xlabel = xattrname.capitalize()
        #Determine the unit, if not given
        if xunit is None:
            xunit = self.get_unit(xattrname) #Automatic unit
        if xunit != "": #Tack on unit, if exists
            xlabel = xlabel +" ["+xunit+"]"
        plt.xlabel(xlabel, fontsize=axisfontsize)

        #y-axis labels (with units), if so desired
        if ylabel is None:
            #Set the y-axis label
            ylabel = yattrname.capitalize()
        #Determine the unit, if not given
        if yunit is None: #For y-axis
            yunit = self.get_unit(yattrname) #Automatic unit
        if yunit != "": #Tack on unit, if exists
            ylabel = ylabel +" ["+yunit+"]"
        plt.ylabel(ylabel, fontsize=axisfontsize)


        ##Below Section: SET title + legend + tick label size
        #For legend, if so desired
        if dolegend:
            plt.legend(loc=legloc, frameon=False, fontsize=legfontsize)

        #For title
        plt.title(title, fontsize=titlefontsize)

        #Set font size of tick labels
        plt.xticks(fontsize=tickfontsize)
        plt.yticks(fontsize=tickfontsize)


        ##Below Section: FINALIZE + EXIT
        #Stop here, if these commands are part of an external plot
        if dopart:
            return

        #Otherwise, save the figure, if so desired
        if dosave:
            plt.savefig(savename)
            plt.close()
            return

        #Otherwise, display the figure
        plt.show()
        return
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
        self._set_attr(attrname="spectrum", attrval=outy_arr,
                        attrunit="Jy") #Line-spec.
        self._set_attr(attrname="emission", attrval=outem_arr,
                        attrunit="Jy") #Em-spec.
        self._set_attr(attrname="continuum", attrval=outcont_arr,
                        attrunit="Jy") #Continuum
        self._set_attr(attrname="wavelength", attrval=outmu_arr,
                        attrunit=r"$\mu$m") #Wavelen. [mu]
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
