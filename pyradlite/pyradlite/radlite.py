##FILE: Radlite
##PURPOSE: Contains the classes RadliteModel() and RadliteSpectrum().


##Below Section: IMPORT necessary packages
import astropy.constants as const
from exutils import func_timer
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
import radmc
try:
    import astropy.io.fits as fitter
except ImportError:
    import pyfits as fitter

#Set helpful constants
dotest = False
if dotest:
    au0 = 1.49597870E13
    c0 = 2.99792458E10
    rsun0     = 6.9599E10
    msun0     = 1.9891E33
    lsun0     = 3.8525E33
    pc0     = 3.085678E18
    h0     = 6.6262000E-27
    kB0     = 1.3807E-16
    mp0     = 1.6726231E-24
    G0     = 6.67259E-8
    tsun0     = 5.78E3
    mu0     = 2.3
    amu0    = 1.6605402E-24
    cinmu0 = c0*1.0E4 #mu/s
    cinkm0 = c0/1.0E5 #km/s
else:
    amu0 = const.u.cgs.value #Atomic mass
    au0 = const.au.cgs.value
    c0 = const.c.cgs.value
    G0 = const.G.cgs.value
    mp0 = const.m_p.cgs.value #Mass of proton
    msun0 = const.M_sun.cgs.value
    rsun0 = const.R_sun.cgs.value
    h0 = const.h.cgs.value
    kB0 = const.k_B.cgs.value
    cinmu0 = c0*1.0E4 #mu/s
    cinkm0 = c0/1.0E5 #km/s
#

#Prepare the welcome message
welcome_message = (
    "--------------------------------------------------\n"
    +"Welcome to RADLite Version 1.3.0!\n\n"
    +"RADLite Version 1.2 (in IDL) was written by:\n"
    +"Klaus Pontoppidan (pontoppi@stsci.edu)\n"
    +"Kees Dullemond\n"
    +"Alex Lockwood\n"
    +"Rowin Meijerink\n\n"
    +"RADLite Version 1.3.0 (in Python) was written by:\n"
    +"Jamila Pegues (jamila.pegues@cfa.harvard.edu)\n"
    +"Klaus Pontoppidan (pontoppi@stsci.edu)\n"
    +"--------------------------------------------------\n\n")



##
class RadliteModel():
    @func_timer
    def __init__(self, infilename, hitranfilename, radmcfilepath):
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
        with open(infilename, 'r') as openfile:
            inputdict = json.load(openfile)
        #Store in secret dictionary (stripping out comments)
        self._attrdict = {}
        for key in inputdict:
            self._attrdict[key] = inputdict[key]["value"]
        self._attrdict["units"] = {} #Inner dictionary to hold units
        #NOTE: Unit dictionary for calculated quantities only! NOT inputs

        #Read in HITRAN data and then extract desired molecule
        with open(hitranfilename, 'r') as openfile:
            hitrandict = json.load(openfile)
        #Store molecular data for specified molecule
        moldict = hitrandict[self._attrdict["molname"]] #Extract mol. data
        for key in moldict: #Store molecule-specific data
            self._attrdict[key] = moldict[key]
        self._attrdict["molmass"] = moldict["molweight"]*amu0 #Mol. mass


        ##Below Section: PRINT the welcome message, if so desired
        if self.get_attr("verbose"):
            print(welcome_message)
            print("")
            print("-"*10+"\n"+"Now preparing all RADLite input files...")
            print("")


        ##Below Section: CHECK inputs for user error
        self._check_inputs()


        ###TEMP. BLOCK START - UNTIL UPDATED TO RUN RADMC/RADMC-3D
        ##Below Section: TEMPORARY - USE PYRADMC TO READ IN RADMC MODEL
        modradmc = radmc.radmc_model(radmcfilepath)
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
        numcores = self.get_attr("numcores") #Number of cores
        #Prepare+run pool of cores
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


    ##CHECK METHODS
    def _check_inputs(self):
        """
        Method: _check_inputs
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Checks that the user has initialized this class instance with
              ...sensible inputs.  Generally we assume that users will use
              ...sensible inputs, so this method does not actually check all
              ...inputs at the moment.  Can easily be updated/changed as
              ...absolutely needed.
        Inputs: N/A
        Outputs: N/A
        Notes:
            > Will raise an error if any of the necessary inputs are missing.
        """
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
            print("Oh no, looks like you requested too many cores (you "
                    +"requested "+str(self.get_attr("numcores"))+")!")
            print("Looks like you only have "+str(maxcores)+" available in "
                    +"total, so we're reducing the number of cores down to "
                    +str(newnumcores)+".")
            self._set_attr(attrname="numcores", attrval=newnumcores)

        #!!! TEMPORARY UNTIL N-LTE SUPPORTED
        if not self.get_attr("lte"):
            raise ValueError("Sorry, currently we only support LTE "
                                +"calculations.  Please set lte to false "
                                +"in the input file.")


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
            > If the desired attribute has not been stored before, then the code
              ...will try to read in its value.  If that doesn't work, then the
              ...code will try to calculate its value.  And if that doesn't
              ...work, then the code will stop trying (see note below).
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
            > Fetches the automatic unit of the given attribute.  Will not
              ...return any units for inputs from initialization.
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
              ...attribute database.
        Inputs: 2 required, 1 optional
            > attrname (required)
                - Type: string
                - Example: "mstar"
                - Description: Name of the desired attribute (such as "mstar"
                  ...for the stellar mass).
            > attrval (required)
                - Type: Varies
                - Example: 1.0E33
                - Description: Desired value for the given attribute.
            > attrunit (optional; default="")
                - Type: string
                - Example: "g"
                - Description: Unit of the given attribute value.
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: RECORD given attribute under given name + EXIT
        self._attrdict[attrname] = attrval
        self._attrdict["units"][attrname] = attrunit #Assign unit
        return
    #



    ##RUN METHODS
    @func_timer
    def run_radlite(self):
        """
        Method: run_radlite
        Purpose:
            > Runs RADLite on the number of cores requested during
              ...initialization of this class.
            > Organizes RADLite input/working/output files in the run directory
              ...specified during initiialization.
        Inputs: N/A
        Outputs: N/A
        Notes:
            > Working directories will be created and deleted during this
              ...process.
            > If the specified run directory already exists, then old files
              ...in that directory will be overwritten!
            > After running run_radlite(), use the gen_spec() method in the
              ...RadliteSpectrum() class to read and process the RADLite output.
        """
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

            #Call core routine
            phere = mp.Process(target=self._run_core,
                                args=(ai, cpudir, outputdir,))
            plist.append(phere)
            #Start this core
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Starting "+str(ai)+"th core in "+cpudir+"...")
            phere.start()
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done starting core(s)!")

        #Close pool of cores
        for ai in range(0, numcores):
            plist[ai].join()
            if self.get_attr("verbose"): #Verbal output, if so desired
                print(str(ai)+"th core finished and closed!")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("All core(s) have finished!")

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
        Method: _run_core
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > A submethod called by the run_radlite() method.
            > Runs RADLite on core #pind.
            > Organizes RADLite input/working/output files for core #pind.
        Inputs: 3 rqeuired
            > cpudir (required)
                - Type: string
                - Example: "workingdir"
                - Description: Path to the working directory for this core.
            > outputdir (required)
                - Type: string
                - Example: "outputdir"
                - Description: Path to the output directory for the RADLite run.
            > pind (required)
                - Type: int
                - Example: 2
                - Description: Index (starting from 0) of the core.
        Outputs: 1 (written, not returned)
            > "RADLITE_core_<pind>.log", "moldata_<pind>.dat",
              ..."linespectrum_moldata_<pind>.dat", and
              ..."levelpop_moldata_<pind>.dat" are placed in outputdir.
        Notes: N/A
        """
        if self.get_attr("verbose"): #Verbal output, if so desired
            print(str(pind)+"th core has started working...")
        ##Below Section: STOP this core if it has no lines to process
        if self.get_attr("_levelpopdicts")[pind] is None: #If no level pop. info
            if self.get_attr("verbose"): #Verbal output, if so desired
                print(str(pind)+"th core does not have any lines to process.  "
                        +"Exiting...")
            return


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
    def plot_attr(self, yattrname, fig=None, figsize=(10,10), markersize=30, linewidth=3, linestyle="-", linecolor="black", markerstyle="o", markercolor="blue", cmap=plt.cm.bone_r, xlog=False, ylog=False, zlog=False, xscaler=1.0, yscaler=1.0, zscaler=1.0, alpha=1.0, xlim=None, ylim=None, xunit=None, yunit=None, zunit=None, xlabel=None, ylabel=None, cbarlabel=None, cbarrotation=270, cbarlabelpad=25, levels=100, axisfontsize=16, titlefontsize=18, legfontsize=16, tickfontsize=14, title="", dolegend=False, leglabel="", legloc="best", dopart=False, dosave=False, savename="testing.png"):

        """
        Method: plot_attr
        Purpose:
            > Plots the given attribute(s).
        Inputs: 1 required, 24 optional
            > yattrname (required)
                - Type: string
                - Example: "velocity_phi"
                - Description: Name of the attribute (such as "velocity_phi" for
                  ...the phi component of the velocity) to either plot as a
                  ...gradient (if 2D) or plot as the y-axis of a line+scatter
                  ...plot (if 1D).
                - Type: integer OR float; in [0,1]
                - Example: 0.5
                - Description: Measure of transparency of the line+scatter plot.
                  ...1 is fully opaque and 0 is fully transparent.
            > axisfontsize (optional; default=16)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the x-axis and y-axis labels and
                  ...for the colorbar label (if yattrname attribute is 2D).
            > cbarlabel (optional; default=None)
                - Type: string OR None
                - Example: "Z-axis Values"
                - Description: This parameter is used only when the attribute
                  ...given by yattrname is 2D.
                  ...It is the label for the colorbar.  If None, then
                  ...the capitalized name of the 2D attribute will be used for
                  ...the label.
            > cbarlabelpad (optional; default=25)
                - Type: integer OR float
                - Example: -5
                - Description: This parameter is used only when the attribute
                  ...given by yattrname is 2D.
                  ...It is the space between the colorbar label and the
                  ...colorbar ticks.  Negative and positive numbers will move
                  ...the label left and right, respectively.
            > cbarrotation (optional; default=270)
                - Type: string OR None
                - Example: 90
                - Description: This parameter is used only when the attribute
                  ...given by yattrname is 2D.
                  ...It is the rotation of the colorbar label in [degrees].
            > cmap (optional; default=plt.cm.bones_r)
                - Type: <matplotlib.pyplot Colormap instance>
                - Example: plt.cm.afmhot_r
                - Description: This should be the name of a colormap available
                  ...from matplotlib.pyplot and will be the colormap used for
                  ...the gradient.  2D case only.
            > dolegend (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will show a legend on the plot.
                - Notes: The value for leglabel will be shown on the legend.
            > dopart (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will exit the function without saving
                  ...or displaying the plot.
                - Notes: Useful for overplotting different attributes and/or
                  ...different class instances onto the same figure, which can
                  ...be done using multiple plot_attr() calls.  Requires that a
                  ...<matplotlib.pyplot Figure instance> be passed for fig.
            > dosave (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will save the plot as savename.
            > fig (optional; default=None)
                - Type: <matplotlib.pyplot Figure instance> OR None
                - Example: matplotlib.pyplot.figure()
                - Description: A figure upon which to plot the given attribute.
                  ...See matplotlib.pyplot documentation for Figure.
            > figsize (optional; default=(10,10))
                - Type: list-like, with 2 values
                - Example: (10,10)
                - Description: 2 values indicating the figure size along the
                  ...x-axis and y-axis.
            > legfontsize (optional; default=16)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the legend (if dolegend=True).
            > leglabel (optional; default="")
                - Type: string OR None
                - Example: "Y-axis Values"
                - Description: The label for the line+scatter plot.  Will be
                  ...displayed only if dolegend=True.
            > legloc (optional; default="best")
                - Type: string
                - Example: "lower left"
                - Description: The location for the plot legend.  It should
                  ...should be a location supported by matplotlib.pyplot.legend.
                  ...Will be used only if dolegend=True.
            > levels (optional; default=50)
                - Type: int OR list-like
                - Example: [10, 20, 30, 40, 50]
                - Description: Contour levels for the gradient.  2D case only.
            > linecolor (optional; default="black")
                - Type: string
                - Example: "orange"
                - Description: This should be the name of a color available
                  ...from matplotlib.pyplot and will be used for the plotted
                  ...line.
            > linestyle (optional; default="-")
                - Type: string
                - Example: "-"
                - Description: Style of the line in the line+scatter plot.
                - Notes: All values accepted by matplotlib.pyplot.plot() are
                ...accepted here.
            > linewidth (optional; default=3)
                - Type: integer OR float
                - Example: 3
                - Description: Thickness of the line in the line+scatter plot.
                - Notes: Set to 0 if no line is desired.
            > markercolor (optional; default="blue")
                - Type: string
                - Example: "blue"
                - Description: Color of the markers for the line+scatter plot.
                - Notes: All values accepted by matplotlib.pyplot are accepted
                  ...here.
            > markersize (optional; default=30)
                - Type: integer OR float
                - Example: 30
                - Description: Size of the markers for the line+scatter plot.
                - Notes: Set to 0 if no markers are desired.
            > markerstyle (optional; default="o")
                - Type: string
                - Example: "o"
                - Description: Style of the markers for the line+scatter plot.
                - Notes: All values accepted by matplotlib.pyplot.scatter() are
                ...accepted here.
            > savename (optional; default="testing.png")
                - Type: string
                - Example: "plot_folder/snazzy_plot.png"
                - Description: Filepath+filename for saving the figure.  Only
                  ...used if dosave=True.
            > tickfontsize (optional; default=14)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the x-axis and y-axis ticks,
                  ...and for the colorbar (if the yattrname attribute is 2D).
            > title (optional; default="")
                - Type: string
                - Example: "A Snazzy Plot"
                - Description: The title for the overall plot.
            > titlefontsize (optional; default=18)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the figure title.
            > xlabel (optional; default=None)
                - Type: string OR None
                - Example: "X-axis Values"
                - Description: The label for the x-axis.  If None, then
                  ...the capitalized name of the x-axis attribute (xattrname)
                  ...will be used for the label.
            > xlim (optional; default=None)
                - Type: list-like with 2 values OR None
                - Example: [50,200]
                - Description: The x-axis range for the plot.  If None, then
                  ...the x-axis range will not be changed.
            > xlog (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will set the x-axis to log-scale.  If
                  ...False, will leave the x-axis as is.
            > xscaler (optional; default=1.0)
                - Type: integer OR float
                - Example: 0.01
                - Description: A factor by which to multiply the x-axis values.
                - Notes: Useful for changing the unit of the x-axis (e.g., a
                ...value of 0.01 to change from cm to m).
            > xunit (optional; default=None)
                - Type: string OR None
                - Example: "centimeters"
                - Description: If xunit is None, then the code will use the
                  ...default unit (if it exists) for the requested x-axis
                  ...attribute.  If xunit is a string, then the code will use
                  ...the given string as part of the label for the x-axis.
                  ...It will be wrapped in square brackets.
                - Notes: Useful if, for example, the user scaled the x-axis with
                  ...xscaler and thus changed the unit of the x-axis values.
                  ...If the user does not want a unit shown at all, then the
                  ...user should set xunit to an empty string ("").
            > ylabel (optional; default=None)
                - Type: string OR None
                - Example: "Y-axis Values"
                - Description: The label for the y-axis.  If None, then
                  ...the capitalized name of the y-axis attribute (yattrname)
                  ...will be used for the label.
            > ylim (optional; default=None)
                - Type: List-like with 2 values OR None
                - Example: [50,200]
                - Description: The y-axis range for the plot.  If None, then
                  ...the y-axis range will not be changed.
            > ylog (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will set the y-axis to log-scale.  If
                  ...False, will leave the y-axis as is.
            > yscaler (optional; default=1.0)
                - Type: integer OR float
                - Example: 0.01
                - Description: This is a factor by which to multiply the y-axis
                  ...values.
                - Notes: Useful for changing the unit of the y-axis
                  ...(e.g., a value of 0.01 to change from cm to m).
            > yunit (optional; default=None)
                - Type: string OR None
                - Example: "radians"
                - Description: If yunit is None, then the code will use the
                  ...default unit (if it exists) for the requested y-axis
                  ...attribute.  If yunit is a string, then the code will use
                  ...the given string as part of the label for the y-axis.
                  ...It will be wrapped in square brackets.
                  ...If the user does not want a unit shown at all, then the
                  ...user should set yunit to an empty string ("").
                - Notes: Useful if, for example, the user scaled the y-axis with
                  ...yscaler and thus changed the unit of the y-axis values.
            > zscaler (optional; default=1.0)
                - Type: integer OR float
                - Example: 0.01
                - Description: This is a factor by which to multiply the z-axis
                  ...values.  2D case only.
                - Notes: Useful for changing the unit of the z-axis
                  ...(e.g., a value of 0.01 to change from cm to m).
            > zlog (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will set the z-axis to log-scale.  If
                  ...False, will leave the z-axis as is.  2D case only.
            > zunit (optional; default=None)
                - Type: string OR None
                - Example: "centimeters/second"
                - Description: If zunit is None, then the code will use the
                  ...default unit (if it exists) for the requested z-axis
                  ...attribute.  If xunit is a string, then the code will use
                  ...the given string as part of the label for the z-axis.
                  ...It will be wrapped in square brackets.  2D case only.
                - Notes: Useful if, for example, the user scaled the z-axis with
                  ...zscaler and thus changed the unit of the z-axis values.
                  ...If the user does not want a unit shown at all, then the
                  ...user should set zunit to an empty string ("").
        Outputs: Nothing OR displays a plot OR saves a plot
        Notes:
            > If a 2D attribute is given for yattrname, then the code will plot
              ...that attribute as a gradient with theta and radius along the
              ...axes.
            > If a 1D attribute is given for yattrname, then the code will plot
              ...that attribute along the y-axis in a line+scatter plot.
            > If a figure is passed for fig and dopart is True, then the figure
              ...will be populated in the background and will not be erased or
              ...shown.
        """

        ##Below Section: INITIALIZE empty plot, if no existing plot given
        if fig is None:
            fig = plt.figure(figsize=figsize)


        ##Below Section: FETCH y-axis values
        yvals = self.get_attr(yattrname)


        ##Below Section: PLOT as either 2D gradient or 1D line+scatter
        ndim = len(np.asarray(yvals).shape) #Number of dimensions for plot
        if ndim == 2: #If 2D quantity (assumed axes are radius vs. theta)
            #Set desired 2D quantity to new variable z
            zattrname = yattrname
            zvals = self.get_attr(zattrname) #2D quantity
            #Extract radius and theta for x and y-axes
            yattrname = "theta"
            yvals = self.get_attr(yattrname) #y-axis values
            xattrname = "radius"
            xvals = self.get_attr(xattrname) #x-axis values
            #For cases where z-axis covers only half of thetas, trim thetas
            if zvals.shape[0] == len(yvals)/2.0:
                halflen = len(yvals)//2
                yvals = yvals[0:halflen] #Cut thetas in half (since mirrored)
            #Plot gradient
            grad = plt.contourf(xvals*xscaler, yvals*yscaler, zvals*zscaler,
                                cmap=cmap, levels=levels)
        elif ndim == 1: #If 1D quantities
            #Generate numerical x-axis
            xvals = np.arange(0, len(yvals))
            #Plot line plot
            plt.plot(xvals*xscaler, yvals*yscaler, color=linecolor,
                    markerfacecolor=markercolor, markeredgecolor=markercolor,
                    markersize=markersize,
                    linewidth=linewidth, linestyle=linestyle,
                    alpha=alpha, label=leglabel, marker=markerstyle)
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
                if zunit is None:
                    zunit = self.get_unit(zattrname) #Automatic unit
                cbarlabel = zattrname.capitalize()
                if zunit != "": #Tack on unit, if exists
                    cbarlabel = cbarlabel + " ["+zunit+"]"
            cbar = plt.colorbar(grad)
            cbar.set_label(label=cbarlabel, fontsize=titlefontsize,
                            rotation=cbarrotation, labelpad=cbarlabelpad)
            cbar.ax.tick_params(labelsize=tickfontsize)


        ##Below Section: LABEL plot axes, if so desired
        #Extract units, if not given
        if yunit is None:
            yunit = self.get_unit(yattrname) #Unit for y-axis
        if xunit is None:
            if ndim == 1: #For 1D case, no unit needed for x-axis
                xunit = ""
            elif ndim == 2: #For 2D case, take automatic unit
                xunit = self.get_unit(xattrname) #Unit for x-axis
        #x-axis labels (with units), if so desired
        if xlabel is None:
            #Set blank x-axis label for 1D case
            if ndim == 1:
                xlabel = ""
            #Otherwise, set attribute name as label for 2D case
            elif ndim == 2:
                xlabel = xattrname.capitalize()
        #Determine the unit for 2D case, if not given
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
        Method: _prep_mol_forcore
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Determines unique levels (e.g., energy) for core number #pind.
            > Calculates level populations for this core.
        Inputs: 1 required
            > pind (required)
                - Type: int
                - Example: 2
                - Description: Index (starting from 0) of this core.
        Outputs: 1
            > <dictionary>
                - Key-Value Pairs:
                    o "Euniq": Array of unique energy levels
                    o "guniq": Array of corresponding degeneracies
                    o "quniq": Array of corresponding rotational levels
                    o "vuniq": Array of corresponding vibrational levels
                    o "npop": Array of level populations
                    o "numlevels": Number of unique energy levels
                    o "numtrans": Number of unique molecular transitions (lines)
                - Description: Dictionary containing level population
                  ...information for this core.
        Notes: N/A
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


        ##Below Section: CALCULATE level populations
        numlevels = len(Euniqarr)
        if numlevels == 0: #If no molecular lines partitioned to this core
            return None #Nothing for this core to do later
        numtrans = len(wavenumarr)
        Euniqarr_K = Euniqarr*h0*c0/1.0/kB0 #[K]; Upper energy levels
        #Fetch interpolated partition sum
        psuminterped = self.get_attr("_psum_interped")
        #Calculate population levels
        npoparr = np.array([(guniqarr[ai]
                                    *np.exp(-1.0*Euniqarr_K[ai]/tempgasarr))
                            for ai in range(0, numlevels)]) /1.0/psuminterped
        #Trim any ridiculously-low values
        npoparr[npoparr < 1E-99] = 0.0


        ##Below Section: STORE + RETURN unique level results + EXIT
        lpopdict = {"npop":npoparr, "Euniq":Euniqarr, "guniq":guniqarr,
                            "quniq":quniqarr, "vuniq":vuniqarr,
                            "numlevels":numlevels, "numtrans":numtrans}
        return lpopdict
    #


    def _prep_mol_forall(self):
        """
        Method: _prep_mol_forall
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Reads in HITRAN molecular data and partition sum information.
            > Splits molecular lines across the number of cores specified by
              ...the user at initialization.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: Determine molecular line and partition sum information
        #Read molecular line and partition data
        self._read_hitran()
        self._read_psum()
        #Set up partition sum calculation
        psum = self.get_attr("_psum")
        psumtemp = self.get_attr("_psum_temp")
        tlen = len(self.get_attr("theta"))//2
        tempgasarr = self.get_attr("gastemperature")[0:tlen,:]
        #Interpolate and store partition sum
        psumfunc = interper(psumtemp, psum, kind=self.get_attr("interpolation"),bounds_error=False,fill_value="extrapolate")
        self._set_attr(attrname="_psum_interped", attrval=psumfunc(tempgasarr))


        ##Below Section: SPLIT data across given number of cores
        numlines = self.get_attr("numlines") #Number of lines
        numcores = self.get_attr("numcores") #Number of cores
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Dividing up a total of "+str(numlines)+" lines for "
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
        Method: _read_hitran
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Reads the HITRAN file for the molecule specified by the user
              ...during initialization.
            > Extracts molecular lines that adhere to the criteria (e.g.,
              ...minimum and maximum abundances) specified by the user during
              ...initialization.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: EXTRACT necessary parameters
        isonum = self.get_attr("isonum")
        wavenumrange = [1.0E4/self.get_attr("max_mu"),
                                            1.0E4/self.get_attr("min_mu")]
        insmin = self.get_attr("min_ins")
        Eupmax = self.get_attr("max_Eup")
        #vupmax = self.get_attr("max_vup")


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

        #Record *all* molecular lines into a convenient dictionary
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

        #Print a note if any signs of incomplete data
        gupis0inds = np.array([False]*len(hitrandict["gup"])) #Initialization
        if 0 in hitrandict["gup"]:
            gupis0inds = (hitrandict["gup"] == 0)
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Whoa!  Incomplete molecular HITRAN data at "
                        +"wavenum="+str(hitrandict["wavenum"][gupis0inds])
                        +".  Skipping these entries...")


        ##Below Section: REMOVE any lines outside of desired range
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Removing molecular lines outside of specified criteria...")

        #Extract boolean indices of lines that fall within criteria
        keepinds = ~np.array(
                    #If gup is 0
                    (gupis0inds)
                    #If incorrect isotopologue
                    | (hitrandict["isonum"] != isonum)
                    #If wavenumber outside of desired wavenumber range
                    | (hitrandict["wavenum"] < wavenumrange[0])
                    | (hitrandict["wavenum"] > wavenumrange[1])
                    #If intensity, level, or upper energy beyond given cutoffs
                    | (hitrandict["ins"] < insmin)
                    | (hitrandict["Eup"] > Eupmax))
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Looks like {0} out of {1} total lines will be kept..."
                    .format(keepinds.sum(), len(keepinds)))

        #If ortho or para requested, then extract indices of ortho or para lines
        if (self.get_attr("oandp")) and (self.get_attr("whichop") != "both"):
            whichop = self.get_attr("whichop")
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("Ortho (o) vs. para (p) lines also requested!  "
                        +"Now extracting only the "+whichop+" lines...")
            #Split up the rotational strings
            qups = hitrandict["qup"]
            qup2s = np.array([qups[ai].split()[1]
                            for ai in range(0, len(qups))]).astype(float)
            qup3s = np.array([qups[ai].split()[2]
                            for ai in range(0, len(qups))]).astype(float)
            #Take last numerical entry in vibrational string
            vup3s = np.array([hitrandict["vup"][ai]
                            .replace(" ","")[-2:len(hitrandict["vup"][ai])]
                            for ai in range(0, len(hitrandict["vup"]))]
                            ).astype(int)
            #Vib. string *should* have up to 6 characters; throw error if not
            if max([len(hitrandict["vup"][ai].replace(" ", ""))
                        for ai in range(0, len(hitrandict["vup"]))]) > 6:
                raise ValueError("SERIOUS READ-IN ERROR HAS OCCURRED!!!  "
                                "PLEASE CONTACT YOUR CODE PROVIDER!!!")
            #Add certain components together to form o. vs. p. criterion
            opnums = qup2s + qup3s + vup3s #Criterion

            #Extract indices of specified form
            if whichop == "o":
                opinds = np.array((opnums % 2) == 1) #If odd, then ortho
            elif whichop == "p":
                opinds = np.array((opnums % 2) == 0) #If even, then para
            else:
                raise ValueError("Oh no!  Please make sure that the user-"
                            +"specified value of 'whichop' in the input file "
                            +"(processed during initialization) is either "
                            +"'o', 'p', or 'both'.")

            #Boolean-combine the o. vs. p. indices with general-criteria indices
            if self.get_attr("verbose"): #Verbal output, if so desired
                print("There are "+str(opinds.sum())+" "+whichop+" lines "
                        +"across the entire database.  So only these will be "
                        +"considered.")
            keepinds = keepinds * opinds #Updated boolean array

        #Note the total number of lines that fall within criteria
        numlines = keepinds.sum()

        #Keep only those lines that fall within criteria; delete all other lines
        for key in hitrandict:
            hitrandict[key] = hitrandict[key][keepinds]
        self._set_attr(attrname="_hitrandict", attrval=hitrandict) #Final set
        self._set_attr(attrname="numlines", attrval=numlines) #Final count
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("There are "+str(numlines)+" molecular lines "
                    +"that fall within ALL specified criteria.")


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done extracting molecular data!\n")
        return
    #


    def _read_mstar(self):
        """
        Method: _read_mstar
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > PURE WRAPPER: Calls _read_starinfo.
            > Allows an avenue for the stellar mass ("mstar") to be 'read-in'
              ...using the get_attr() framework.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: CALL general stellar info reader
        self._read_starinfo()
        return
    #


    def _read_psum(self):
        """
        Method: _read_psum
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Reads in the partition sum data from the psumfile specified by the
              ...user during initialization.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
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
        self._set_attr(attrname="_psum", attrval=psumarr)
        self._set_attr(attrname="_psum_temp", attrval=psumtemparr,
                        attrunit="K")
        self._set_attr(attrname="_psum_molmass", attrval=molmassarr[iind],
                        attrunit="g")
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Partition information for "+self.get_attr("molname")+" "
                                        +"has been successfully extracted!\n")
        return
    #


    def _read_rstar(self):
        """
        Method: _read_rstar
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > PURE WRAPPER: Calls _read_starinfo.
            > Allows an avenue for the stellar radius ("rstar") to be 'read-in'
              ...using the get_attr() framework.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: CALL general stellar info reader
        self._read_starinfo()
        return
    #


    def _read_starinfo(self):
        """
        Method: _read_starinfo
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Reads in stellar information from the RADMC output file
              ...("starinfo.inp").
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
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
        Method: _read_teff
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > PURE WRAPPER: Calls _read_starinfo.
            > Allows an avenue for the stellar effective temperature ("teff")
              ...to be 'read-in' using the get_attr() framework.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: CALL general stellar info reader
        self._read_starinfo()
        return
    #



    ##CALCULATION METHODS
    def _calc_abundance(self):
        """
        Method: _calc_abundance
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Calculates the abundance of the molecule given the criteria
              ...(e.g., freeze-out temperature) specified by the user during
              ...initialization.
              ...If the freeze-out temperature is None, then the abundance will
              ...be constant across the source (set to "max_abund").  If the
              ...freeze-out temperature is a number, then the abundance will be
              ...set to "max_abund" at temperatures >= freeze-out temperature
              ...and set to "min_abund" at temperatures < freeze-out
              ...temperature.
        Inputs: N/A
        Outputs: N/A
        Notes:
            > This method could most certainly be overwritten in a super-class
              ...(which would inherit this base class).
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


        ##Below Section: RECORD calculated abundance + EXIT
        self._set_attr(attrname="abundance", attrval=abundarr)
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done calculating abundance!\n")
        return
    #


    def _calc_gasdensity(self):
        """
        Method: _calc_gasdensity
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Calculates the gas density of the source, assuming that the gas
              ...density is equal to (dust density * ("gastodust")).
        Inputs: N/A
        Outputs: N/A
        Notes:
            > This method could most certainly be overwritten in a super-class
              ...(which would inherit this base class).
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


    def _calc_gastemperature(self):
        """
        Method: _calc_gastemperature
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Calculates the gas temperature of the source, assuming that it
              ...is equal to the dust temperature.
        Inputs: N/A
        Outputs: N/A
        Notes:
            > This method could most certainly be overwritten in a super-class
              ...(which would inherit this base class).
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


    def _calc_turbulence(self):
        """
        Method: _calc_turbulence
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Calculates the turbulence of the source, assuming that it is equal
              ...to the scaled sound speed and thermal broadening added in
              ...quadrature.
        Inputs: N/A
        Outputs: N/A
        Notes:
            > This method could most certainly be overwritten in a super-class
              ...(which would inherit this base class).
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


    def _calc_velocity(self):
        """
        Method: _calc_abundance
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Calculates the gas velocity of the source, assuming that the gas
              ...is Keplerian.
              ...Three components are calculated: phi ("velocity_phi"),
              ...theta ("velocity_theta"), and radial ("velocity_radial").
        Inputs: N/A
        Outputs: N/A
        Notes:
            > This method could most certainly be overwritten in a super-class
              ...(which would inherit this base class).
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
        Method: _calc_velocity_radial
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > PURE WRAPPER: Calls _calc_velocity.
            > Allows an avenue for the radial velocity component
              ...("velocity_radial") to be 'calculated' using the get_attr()
              ...framework.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: CALL general velocity calculator
        self._calc_velocity()
        return
    #


    def _calc_velocity_theta(self):
        """
        Method: _calc_velocity_theta
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > PURE WRAPPER: Calls _calc_velocity.
            > Allows an avenue for the theta velocity component
              ...("velocity_theta") to be 'calculated' using the get_attr()
              ...framework.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: CALL general velocity calculator
        self._calc_velocity()
        return
    #


    def _calc_velocity_phi(self):
        """
        Method: _calc_velocity_phi
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > PURE WRAPPER: Calls _calc_velocity.
            > Allows an avenue for the phi velocity component
              ...("velocity_phi") to be 'calculated' using the get_attr()
              ...framework.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
        """
        ##Below Section: CALL general velocity calculator
        self._calc_velocity()
        return
    #



    ##WRITE METHODS
    def _write_abundanceinp(self):
        """
        Method: _write_abundanceinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "abundance.inp" input file for RADLite.
        Inputs: N/A
        Outputs: 1 (written, not returned)
            > "abundance.inp" is written to the inp_path directory specified
              ...during initialization.
        Notes: N/A
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
        Method: _write_core_levelpopinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "levelpop_moldata.dat" and "levelpop.info" input
              ...files for core #pind to use for running RADLite.
        Inputs: 2 required
            > cpudir (required)
                - Type: string
                - Example: "workingdir"
                - Description: Path to the working directory for this core.
            > pind (required)
                - Type: int
                - Example: 2
                - Description: Index (starting from 0) of the core.
        Outputs: 2 (written, not returned)
            > "levelpop_moldata.dat" and "levelpop.info" are written to cpudir.
        Notes: N/A
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


        ##Below Section: EXIT method
        return
    #


    def _write_core_linespectruminp(self, pind, cpudir):
        """
        Method: _write_core_linespectruminp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "linespectrum.inp" input file for core #pind to use
              ...for running RADLite.
        Inputs: 2 required
            > cpudir (required)
                - Type: string
                - Example: "workingdir"
                - Description: Path to the working directory for this core.
            > pind (required)
                - Type: int
                - Example: 2
                - Description: Index (starting from 0) of the core.
        Outputs: 1 (written, not returned)
            > "linespectrum.inp" is written to cpudir.
        Notes: N/A
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


    def _write_core_moldatadat(self, cpudir, pind):
        """
        Method: _write_core_moldatadat
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "moldata.data" input files for core #pind to use
              ...for running RADLite.
        Inputs: 2 required
            > cpudir (required)
                - Type: string
                - Example: "workingdir"
                - Description: Path to the working directory for this core.
            > pind (required)
                - Type: int
                - Example: 2
                - Description: Index (starting from 0) of the core.
        Outputs: 1 (written, not returned)
            > "moldata.dat" is written to cpudir.
        Notes: N/A
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


    def _write_densityinp(self):
        """
        Method: _write_densityinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "density.inp" input file for RADLite.
        Inputs: N/A
        Outputs: 1 (written, not returned)
            > "density.inp" is written to the inp_path directory specified
              ...during initialization.
        Notes: N/A
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


    def _write_gastemperatureinp(self):
        """
        Method: _write_gastemperatureinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "temperature.inp" input file for RADLite.
        Inputs: N/A
        Outputs: 1 (written, not returned)
            > "temperature.inp" is written to the inp_path directory specified
              ...during initialization.
        Notes: N/A
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


    def _write_radliteinp(self):
        """
        Method: _write_radliteinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "radlite.inp" input file for RADLite.
        Inputs: N/A
        Outputs: 1 (written, not returned)
            > "radlite.inp" is written to the inp_path directory specified
              ...during initialization.
        Notes: N/A
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


    def _write_turbulenceinp(self):
        """
        Method: _write_turbulenceinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "turbulence.inp" input file for RADLite.
        Inputs: N/A
        Outputs: 1 (written, not returned)
            > "turbulence.inp" is written to the inp_path directory specified
              ...during initialization.
        Notes: N/A
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


    def _write_velocityinp(self):
        """
        Method: _write_velocityinp
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Writes the "velocity.inp" input file for RADLite.
        Inputs: N/A
        Outputs: 1 (written, not returned)
            > "velocity.inp" is written to the inp_path directory specified
              ...during initialization.
        Notes: N/A
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
class RadliteSpectrum():
    @func_timer
    def __init__(self, infilename):
        """
        Method: __init__
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Initializes an instance of the RadliteSpectrum() class.
            > Checks that input file (infilename) has minimum required
              ...parameters to run the gen_spec() method.
        Inputs: 1 required
            > infilename (required)
                - Type: string
                - Example: "/User/path/to/file/input_spectrum.json"
                - Description: A .json file containing ALL input parameters
                  ...for this class instance.
        Outputs: N/A
        Notes:
            > Creates an underlying dictionary to hold all input parameters...
              ...and future attributes.
        """
        ##Below Section: READ IN + STORE input file
        #Read in input spectrum data
        with open(infilename, 'r') as openfile:
            inputdict = json.load(openfile)
        #Store in secret dictionary (stripping out comments)
        self._attrdict = {}
        for key in inputdict:
            self._attrdict[key] = inputdict[key]["value"]
        self._attrdict["units"] = {} #Inner dictionary to hold units
        #NOTE: Unit dictionary for calculated quantities only! NOT inputs


        ##Below Section: PRINT the welcome message, if so desired
        if self.get_attr("verbose"):
            print(welcome_message)
            print("")


        ##Below Section: EXIT
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Welcome!  You have successfully initialized an instance "
                    +"of RadliteSpectrum(). You can use this instance to "
                    +"process and plot RADLite output. Start by running the "
                    +"gen_spec() method to process RADLite spectra.\n")
        return
    #


    ##CHECK METHODS
    def _check_inputs(self):
        """
        Method: _check_inputs
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Checks that the user has initialized this class instance with
              ...sensible inputs.  Generally we assume that users will use
              ...sensible inputs, so this method does not actually check all
              ...inputs at the moment.  Can easily be updated/changed as
              ...absolutely needed.
        Inputs: N/A
        Outputs: N/A
        Notes:
            > Will raise an error if any of the necessary inputs are missing.
        """
        ##Below Section: CHECK user inputs; make sure they are valid
        #Truncate number of cores, if too many requested
        maxcores = mp.cpu_count()
        if self.get_attr("numcores") > (maxcores - 1):
            newnumcores = max([(maxcores - 1), 1])
            print("Oh no, looks like you requested too many cores (you "
                    +"requested "+str(self.get_attr("numcores"))+")!")
            print("Looks like you only have "+str(maxcores)+" available in "
                    +"total, so we're reducing the number of cores down to "
                    +str(newnumcores)+".")
            self._set_attr(attrname="numcores", attrval=newnumcores)
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
                - Example: "dist"
                - Description: Name of the given attribute (such as "dist" for
                  ...distance to the source).
        Outputs: 1
            > <Attribute value>
                - Type: Varies
                - Example: 140
                - Description: The value of the given attribute.
        Notes:
            > If the desired attribute has not been stored before, then the code
              ...will try to read in its value.  If that doesn't work, then
              ...the code will stop trying (see note below).
            > If an invalid name is given, then the code will return an error
              ...listing all available attributes.
        """
        ##Below Section: DETERMINE requested attribute
        #Try accessing it
        try:
            return self._attrdict[attrname]
        except KeyError: #If attribute not yet recorded...
            pass

        #Otherwise, raise an error
        raise AttributeError("'"+attrname+"' doesn't seem to be a valid "
                            +"attribute.  Valid attributes (sorted) are:\n"
                            +str(np.sort([key for key in self._attrdict]))+".\n"
                            +"Run the method gen_spec() (if you haven't "
                            +"yet) to automatically populate more "
                            +"attributes.\n")
    #


    def get_unit(self, attrname):
        """
        Method: get_unit
        Purpose:
            > Fetches the automatic unit of the given attribute.  Will not
              ...return any units for inputs from initialization.
        Inputs: 1 required
            > attrname (required)
                - Type: string
                - Example: "spectrum"
                - Description: Name of the given attribute (such as "spectrum"
                  ...for the spectrum).
        Outputs: 1
            > <Attribute value>
                - Type: string
                - Example: "Jy"
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
                            +"the method gen_spec() (if you haven't "
                            +"yet) to automatically populate more "
                            +"attributes.\n")
    #


    def _set_attr(self, attrname, attrval, attrunit=""):
        """
        Method: _set_attr
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Sets the given value for the given attribute in the underlying
              ...attribute database.
        Inputs: 2 required, 1 optional
            > attrname (required)
                - Type: string
                - Example: "mstar"
                - Description: Name of the desired attribute (such as "mstar"
                  ...for the stellar mass).
            > attrval (required)
                - Type: Varies
                - Example: 1.0E33
                - Description: Desired value for the given attribute.
            > attrunit (optional; default="")
                - Type: string
                - Example: "g"
                - Description: Unit of the given attribute value.
        Outputs: N/A
        Notes: N/A
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
        Method: gen_spec
        Purpose:
            > Reads in and processes data from RADLite runs specified during
              ...initialization.
            > Puts the molecular lines together into a single spectrum
              ...(removing any duplicate lines).
        Inputs: N/A
        Outputs: N/A
        Notes:
            > After gen_spec() completes, the spectra can be accessed using the
              ...get_attr() method, like so:
              - get_attr('wavelength') for the wavelengths
              - get_attr('frequency') for the frequencies
              - get_attr('spectrum') for the full spectrum
              - get_attr('emission') for the emission only
              - get_attr('continuum') for the continuum
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
            print("get_attr('frequency') for the frequencies")
            print("You can also plot them using the plot_spec() method, "
                    +"using the same keywords just written above.\n")
        return
    #



    ##OUTPUT DISPLAY METHODS
    @func_timer
    def plot_spec(self, yattrname, xattrname="wavelength", fig=None, figsize=(10,10), linewidth=3, linestyle="-", markersize=0, markerstyle="o", markercolor="blue", linecolor="black", xlog=False, ylog=False, xscaler=1.0, yscaler=1.0, alpha=1.0, xlim=None, ylim=None, xunit=None, yunit=None, xlabel=None, ylabel=None, axisfontsize=16, titlefontsize=18, legfontsize=16, tickfontsize=14, title="", dolegend=False, leglabel="", legloc="best", dopart=False, dosave=False, savename="testing.png"):

        """
        Method: plot_attr
        Purpose:
            > Plots the given attribute(s).
        Inputs: 1 required, 24 optional
            > yattrname (required)
                - Type: string
                - Example: "spectrum"
                - Description: Name of the attribute (such as "spectrum" for
                  ...the full spectrum) to plot as a line+scatter plot.
            > xattrname (optional; default="wavelength")
                - Type: string
                - Example: "frequency"
                - Description: Name of the attribute (such as "frequency" for
                  ...the frequencies) to plot along the x-axis.
            > alpha (optional; default=1.0)
                - Type: integer OR float; in [0,1]
                - Example: 0.5
                - Description: Measure of transparency of the line+scatter plot.
                  ...1 is fully opaque and 0 is fully transparent.
            > axisfontsize (optional; default=16)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the x-axis and y-axis labels and
                  ...for the colorbar label (if yattrname attribute is 2D).
            > dolegend (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will show a legend on the plot.
                - Notes: The value for leglabel will be shown on the legend.
            > dopart (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will exit the function without saving
                  ...or displaying the plot.
                - Notes: Useful for overplotting different attributes and/or
                  ...different class instances onto the same figure, which can
                  ...be done using multiple plot_spec() calls.  Requires that a
                  ...<matplotlib.pyplot Figure instance> be passed for fig.
            > dosave (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will save the plot as savename.
            > fig (optional; default=None)
                - Type: <matplotlib.pyplot Figure instance> OR None
                - Example: matplotlib.pyplot.figure()
                - Description: A figure upon which to plot the given attribute.
                  ...See matplotlib.pyplot documentation for Figure.
            > figsize (optional; default=(10,10))
                - Type: list-like, with 2 values
                - Example: (10,10)
                - Description: 2 values indicating the figure size along the
                  ...x-axis and y-axis.
            > legfontsize (optional; default=16)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the legend (if dolegend=True).
            > leglabel (optional; default="")
                - Type: string OR None
                - Example: "Y-axis Values"
                - Description: The label for the line+scatter plot.  Will be
                  ...displayed only if dolegend=True.
            > legloc (optional; default="best")
                - Type: string
                - Example: "lower left"
                - Description: The location for the plot legend.  It should
                  ...should be a location supported by matplotlib.pyplot.legend.
                  ...Will be used only if dolegend=True.
            > linecolor (optional; default="black")
                - Type: string
                - Example: "black"
                - Description: This should be the name of a color available
                  ...from matplotlib.pyplot.
            > linestyle (optional; default="-")
                - Type: string
                - Example: "-"
                - Description: Style of the line in the line+scatter plot.
                - Notes: All values accepted by matplotlib.pyplot.plot() are
                ...accepted here.
            > linewidth (optional; default=3)
                - Type: integer OR float
                - Example: 3
                - Description: Thickness of the line in the line+scatter plot.
                - Notes: Set to 0 if no line is desired.
            > markercolor (optional; default="blue")
                - Type: string
                - Example: "blue"
                - Description: Color of the markers for the line+scatter plot.
                - Notes: All values accepted by matplotlib.pyplot are accepted
                  ...here.
            > markersize (optional; default=30)
                - Type: integer OR float
                - Example: 30
                - Description: Size of the markers for the line+scatter plot.
                - Notes: Set to 0 if no markers are desired.
            > markerstyle (optional; default="o")
                - Type: string
                - Example: "o"
                - Description: Style of the markers for the line+scatter plot.
                - Notes: All values accepted by matplotlib.pyplot.scatter() are
                ...accepted here.
            > savename (optional; default="testing.png")
                - Type: string
                - Example: "plot_folder/snazzy_plot.png"
                - Description: Filepath+filename for saving the figure.  Only
                  ...used if dosave=True.
            > tickfontsize (optional; default=14)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the x-axis and y-axis ticks.
            > title (optional; default="")
                - Type: string
                - Example: "A Snazzy Plot"
                - Description: The title for the overall plot.
            > titlefontsize (optional; default=18)
                - Type: integer OR float
                - Example: 10
                - Description: The fontsize for the figure title.
            > xlabel (optional; default=None)
                - Type: string OR None
                - Example: "X-axis Values"
                - Description: The label for the x-axis.  If None, then
                  ...the capitalized name of the x-axis attribute (xattrname)
                  ...will be used for the label.
            > xlim (optional; default=None)
                - Type: list-like with 2 values OR None
                - Example: [50,200]
                - Description: The x-axis range for the plot.  If None, then
                  ...the x-axis range will not be changed.
            > xlog (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will set the x-axis to log-scale.  If
                  ...False, will leave the x-axis as is.
            > xscaler (optional; default=1.0)
                - Type: integer OR float
                - Example: 0.01
                - Description: A factor by which to multiply the x-axis values.
                - Notes: Useful for changing the unit of the x-axis (e.g., a
                ...value of 0.01 to change from cm to m).
            > xunit (optional; default=None)
                - Type: string OR None
                - Example: "GHz"
                - Description: If xunit is None, then the code will use the
                  ...default unit (if it exists) for the requested x-axis
                  ...attribute.  If xunit is a string, then the code will use
                  ...the given string as part of the label for the x-axis.
                  ...It will be wrapped in square brackets.
                - Notes: Useful if, for example, the user scaled the x-axis with
                  ...xscaler and thus changed the unit of the x-axis values.
                  ...If the user does not want a unit shown at all, then the
                  ...user should set xunit to an empty string ("").
            > ylabel (optional; default=None)
                - Type: string OR None
                - Example: "Y-axis Values"
                - Description: The label for the y-axis.  If None, then
                  ...the capitalized name of the y-axis attribute (yattrname)
                  ...will be used for the label.
            > ylim (optional; default=None)
                - Type: List-like with 2 values OR None
                - Example: [50,200]
                - Description: The y-axis range for the plot.  If None, then
                  ...the y-axis range will not be changed.
            > ylog (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will set the y-axis to log-scale.  If
                  ...False, will leave the y-axis as is.
            > yscaler (optional; default=1.0)
                - Type: integer OR float
                - Example: 1000
                - Description: A factor by which to multiply the y-axis values.
                - Notes: Useful for changing the unit of the y-axis
                  ...(e.g., a value of 1000 to change from Jy to mJy).
            > yunit (optional; default=None)
                - Type: string OR None
                - Example: "microJanskies"
                - Description: If yunit is None, then the code will use the
                  ...default unit (if it exists) for the requested y-axis
                  ...attribute.  If yunit is a string, then the code will use
                  ...the given string as part of the label for the y-axis.
                  ...It will be wrapped in square brackets.
                  ...If the user does not want a unit shown at all, then the
                  ...user should set xunit to an empty string ("").
                - Notes: Useful if, for example, the user scaled the y-axis with
                  ...yscaler and thus changed the unit of the y-axis values.
        """
        ##Below Section: INITIALIZE empty plot, if no existing plot given
        if fig is None:
            fig = plt.figure(figsize=figsize)


        ##Below Section: FETCH x and y-axis values
        yvals = self.get_attr(yattrname)
        xvals = self.get_attr(xattrname)


        ##Below Section: PLOT the desired attributes
        plt.plot(xvals*xscaler, yvals*yscaler, color=linecolor,
                    markerfacecolor=markercolor, markeredgecolor=markercolor,
                    markersize=markersize, marker=markerstyle,
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
        Method: _read_core_radliteoutput
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Reads in and processes a single set of molecular line files output
              ...by RADLite.
        Inputs: 1 required
            > filedict (required)
                - Type: dictionary
                - Key-Value Pairs:
                    o "spec": "linespectrum_moldata.dat" file produced by one
                      ...core during one RADLite run.
                    o "mol": "moldata.dat" file produced by one core during one
                      ...RADLite run.
                - Description: Dictionary containing molecular line and
                  ...data fileset produced by one core during one RADLite run.
        Outputs: 1
            > <dictionary>
                - Key-Value Pairs:
                    o "line": Data extracted from "linespectrum_moldata.dat" file produced by one core during one RADLite run.
                    o "mol": Data extracted from "moldata.dat" file produced by one core during one RADLite run.
                - Description: Dictionary containing processed molecular line
                  ...and molecular info from one fileset produced by one core
                  ...during one RADLite run.
        Notes: N/A
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
        molname = moldata[1].split("\n")[0] #Name of molecule (newline stripped)
        molweight = float(moldata[3]) #Weight of molecule
        numlevels = int(moldata[5]) #Number of energy levels
        numtrans = int(moldata[7+numlevels+1]) #Number of transitions

        #Extract level information
        levhere = moldata[7:7+numlevels] #Section containing all levels
        gbank = np.array([lhere[17:(17+7)] for lhere in levhere]
                            ).astype(float) #All possible degeneracies
        vbank = np.array([lhere[24:(24+15)]
                            for lhere in levhere]) #All possible vib. levels
        qbank = np.array([lhere[39:(39+15)]
                            for lhere in levhere]) #All possible rot. levels

        #Extract transitions
        iloc = 7 + numlevels + 3 #Starting index of transitions
        sechere = moldata[iloc:(iloc+numtrans)] #Section containing all trans.
        sechere = [lhere.split() for lhere in sechere] #Split out spaces
        #For upper and lower degeneracies
        gup_arr = np.array([gbank[int(lh[1])-1]
                            for lh in sechere]).astype(float)
        glow_arr = np.array([gbank[int(lh[2])-1]
                            for lh in sechere]).astype(float)
        #For einstein coefficients, central wavenumbers, and upper energies
        A_arr = np.array([lh[3] for lh in sechere]).astype(float)
        wavenum_arr = np.array([lh[4] for lh in sechere]).astype(float)
        Eup_arr = np.array([lh[5] for lh in sechere]).astype(float) #wavenum.
        EupK_arr = Eup_arr*h0*c0/1.0/kB0 #Kelvin
        #For vibrational and rotational levels
        vup_list = np.array([vbank[int(lh[1])-1] for lh in sechere])
        vlow_list = np.array([vbank[int(lh[2])-1] for lh in sechere])
        qup_list = np.array([qbank[int(lh[1])-1] for lh in sechere])
        qlow_list = np.array([qbank[int(lh[2])-1] for lh in sechere])


        #FOR SPECTRUM FILE
        #Extract all spectral traits from spectrum file
        numlines = int(specdata[4]) #Number of mol. lines
        incl = float(specdata[6].split()[2]) #Velocity info
        maxnumpoints = int(specdata[5]) #Max. number of data points per line

        #Iterate through molecular lines in subrun
        linedict_list = [{} for ai in range(0, numlines)] #To hold all fluxes
        moldict_list = [{} for ai in range(0, numlines)] #To hold all mol. info
        iloc = 6 #Starting index within line files
        for ai in range(0, numlines):
            freq = c0/(1.0/wavenum_arr[ai]) #Central frequency [Hz]
            wavelength = c0/1.0/freq #Central wavelength [cm]
            numpoints = int(specdata[iloc+5]) #Num of freq.

            #Extract section of data for current molecular line
            sechere = specdata[(iloc+7):(iloc+7+numpoints)]
            #Convert into 2D array of floats
            sechere = np.array([lhere.split()
                                for lhere in sechere]).astype(float)

            #Record data for this mol. line
            linedict_list[ai]["vel"] = sechere[:,0] #Velocities
            linedict_list[ai]["flux"] = sechere[:,1] #Fluxes
            linedict_list[ai]["freq"] = freq #Central frequency
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
            #For einstein coeff., central wavenumber, central frequency,
            #...central wavelength, upper energy
            moldict_list[ai]["A"] = A_arr[ai]
            moldict_list[ai]["wavenum"] = wavenum_arr[ai]
            moldict_list[ai]["wavelength"] = wavelength
            moldict_list[ai]["freq"] = freq
            moldict_list[ai]["Eup"] = Eup_arr[ai] #In wavenumber
            moldict_list[ai]["Eup_K"] = EupK_arr[ai] #In Kelvin
            #For vibrational and rotational levels
            moldict_list[ai]["vup"] = vup_list[ai]
            moldict_list[ai]["vlow"] = vlow_list[ai]
            moldict_list[ai]["qup"] = qup_list[ai]
            moldict_list[ai]["qlow"] = qlow_list[ai]


        ##Below Section: RETURN mol. line data + EXIT
        return {"line":linedict_list, "mol":moldict_list}
    #


    @func_timer
    def _read_radliteoutput(self):
        """
        Method: _read_radliteoutput
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Reads in and processes sets of molecular line files produced by
              ...RADLite for all runs specified by user during initialization.
            > Removes any duplicate line occurrences encountered from the given
              ...files.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
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
                raise ValueError("Oh no! "+pathhere+" does not contain an "
                        +"equal # of spectrum ('linespectrum_moldata_*.dat') "
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
                            +dhere["vup"]+dhere["vlow"]
                            +dhere["qup"]+dhere["qlow"]
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
            raise ValueError("Oh no! The RADLite runs you specified don't all "
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
        Method: _process_spectrum
        WARNING: THIS METHOD IS NOT INTENDED FOR DIRECT USE BY USER.
        Purpose:
            > Synthesizes spectrum from processed RADLite output.
            > Interpolates the continuum from the molecular lines.
            > Records the full spectrum, the emission-only spectrum, and the
              ...continuum.
        Inputs: N/A
        Outputs: N/A
        Notes: N/A
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
        freq_list = np.array([dhere["freq"]
                                    for dhere in linedict_list]) #Center freqs.
        freq_list = freq_list[freq_list != 0] #Cut out any freq=0
        #For min. mu
        muminraw = np.min(cinmu0/1.0/freq_list) #Unbroadened min. mu value
        mumin = muminraw - (boxwidth*1.0E9*muminraw/cinmu0) #Broad. by box width
        #For max. mu
        mumaxraw = np.max(cinmu0/1.0/freq_list) #Unbroadened max. mu value
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
                    *1.0E9/dhere["freq"]) + (cinmu0/1.0/dhere["freq"]))
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
        tempindsort = np.argsort([(cinmu0/1.0/dhere["freq"])
                                for dhere in linedict_list])
        tempmusort = np.array([(cinmu0/1.0/dhere["freq"])
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
        #For frequencies
        outfreq_arr = c0/(outmu_arr*1.0E-4) #Hz


        ##Below Section: STORE spectra and molecule information + EXIT
        self._set_attr(attrname="spectrum", attrval=outy_arr,
                        attrunit="Jy") #Line-spec.
        self._set_attr(attrname="emission", attrval=outem_arr,
                        attrunit="Jy") #Em-spec.
        self._set_attr(attrname="continuum", attrval=outcont_arr,
                        attrunit="Jy") #Continuum
        self._set_attr(attrname="wavelength", attrval=outmu_arr,
                        attrunit=r"$\mu$m") #Wavelength [mu]
        self._set_attr(attrname="frequency", attrval=outfreq_arr,
                        attrunit="Hz") #Frequency [Hz]
        if self.get_attr("verbose"): #Verbal output, if so desired
            print("Done processing RADLite spectrum!\n")
        return
    #



    ##WRITE METHODS
    @func_timer
    def write_fits(self, fitsname, overwrite=False):
        """
        Method: write_fits
        Purpose:
            > Writes the spectrum data generated by the gen_spec() method to
              ...a .fits file.
        Inputs: 1 required, 1 optional
            > fitsname (required)
                - Type: string
                - Example: "folder/testfits.fits"
                - Description: Name of the output .fits file.
            > overwrite (optional; default=False)
                - Type: boolean
                - Example: True
                - Description: If True, will overwrite any previous files that
                  ...have the same file name as that given by fitsname.  If
                  ...False, will not overwrite any previous files.
        Outputs: 1 (written, not returned)
            > The .fits file will be written to the given fitsname.
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
            wavelen_arr = self.get_attr("wavelength") #Wavelength [mu]
            freq_arr = self.get_attr("frequency") #Frequency [Hz]
            moldict_list = self.get_attr("molinfo") #All mol. info
        except AttributeError: #Throw an error if hasn't been processed yet
            raise AttributeError("Oh no! Looks like you haven't processed "
                        +"any RADLite output data yet. You can do so by "
                        +"running the gen_spec() method for this class.")


        ##Below Section: #ASSEMBLE data for .fits file
        #Create primary header for .fits file
        hdr = fitter.Header()
        hdr["WAVEUNIT"] = "micron"
        hdr["FREQUNIT"] = "Hz"
        hdr["FLUXUNIT"] = "Jy"
        hdr["INCL_deg"] = moldict_list[0]["incl"]
        hdr["DIST_pc"] = self.get_attr("dist")
        hdu_hdr = fitter.PrimaryHDU(header=hdr)

        #Create and fill table container with flux data
        c1 = fitter.Column(name="wavelength", array=wavelen_arr, format='D')
        c2 = fitter.Column(name="frequency", array=freq_arr, format='D')
        c3 = fitter.Column(name="emission", array=fluxem_arr, format='D')
        c4 = fitter.Column(name="spectrum", array=fluxspec_arr, format='D')
        c5 = fitter.Column(name="continuum", array=fluxcont_arr, format='D')
        hdu_flux = fitter.BinTableHDU.from_columns([c1, c2, c3, c4, c5])

        #Determine format types for all molecular data
        namelist = [key for key in moldict_list[0]] #Key-names of mol. data
        formatlist = ['D']*len(namelist) #Init. formats as floats
        #Change format to strings as applicable
        for ai in range(0, len(namelist)):
            key = namelist[ai] #Current name of molecular data
            if isinstance(moldict_list[0][key], str): #If string, change format
                formatlist[ai] = 'A' #String format (instead of float format)

        #Create and fill table container with molecular data
        arrlist = [np.array([dhere[key] for dhere in moldict_list])
                    for key in namelist] #Arrays corresponding to data names
        collist = [fitter.Column(name=namelist[ai], array=arrlist[ai],
                                    format=formatlist[ai])
                    for ai in range(0, len(namelist))] #Data columns
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
