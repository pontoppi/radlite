##FILE:
##PURPOSE:


##Below Section: IMPORT necessary functions
import subprocess
import multiprocessing as mp
import numpy as np


##
class Radlite():
    def __init__(verbose=True):
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
        ##Below Section: PRINT a welcome message, if so desired
        if verbose:
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


        ##Below Section: Check input parameters for any user error
        #!!!Variables to add: turbmode = ["kep", "sou"]
        #Make sure that desired image cube output is valid
        validimage = [0, 2]
        if image not in validimage:
            raise ValueError("Sorry, the image you have chosen ("
                    +str(image)+") is not a valid image.  The value "
                    +"of image must be one of the following: "
                    +str(validimage)+".")
        ##Below Section: CHOOSE RADLite exec. based on desired output
        #Prepare executable for desired image cube output
        if self.image == 0: #If desired cube output is spectrum
            if verbose: #Verbal output, if so desired
                print("Will prepare a spectrum-formatted image cube...")
                print("")
            self.executable = self.exe_path+"RADlite"
        elif self.image == 2: #Else, if desired cube output is circular
            if verbose: #Verbal output, if so desired
                print("Will prepare a circular-formatted image cube...")
                print("")
            self.executable = self.exe_path+"RADlite_imcir"
        else: #Throw error, otherwise
            pass #Fix later


        #Make sure that desired LTE format is valid
        if not isinstance(lte, bool):
            raise ValueError("Sorry, the lte you have chosen ("
                    +str(lte)+") is not a valid lte.  The value "
                    +"of lte must be a boolean (True or False).")


        ##Below Section: RECORD inputs as attributes


    def set_numprocessors(self, numprocessors):
        pass

    def get_abundance():
        pass

    def get_velocity():
        pass

    def get_gastemperature():
        ##Below Section: RECORD physical structure of underlying model
        radius = self.radmc.radius #Radius
        theta = self.radmc.theta #Theta
        dustdens = self.radmc.dustdensity #Dust density
        dusttemp = self.radmc.dusttemperature #Dust temperature

    def get_moldata():
        pass

    def get_levelpop():
        pass

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
        ##Below Section: Throws an error if invalid number of processors
        #!!!


        ##Below Section: SET UP result directory
        printstamp = time.datetoprint?() #???
        printtime = time.timetoprint?() #???
        if self.nodate: #If date should not be appended to directory name
            rundir = self.run_name + printstamp
        else: #If date should be appended to directory name
            rundir = self.run_name + printstamp + printtime

        #Make the desired directory, if nonexistent
            try:
                comm = subprocess.call(["mkdir", rundir]) #Create directory
            except (comm != 0): #Will override files within, otherwise
                pass
        if verbose: #Verbal output, if so desired
            print("All input files and final data will be saved to the "
                    +"following directory: "+rundir)

        #Copy initial files into final directory
        if verbose: #Verbal output, if so desired
            print("Copying over initial data files into "+rundir+"...")
        initfilelist = ["problem_params.pro", "line_params.ini", "radius.inp",
                        "theta.inp", "frequency.inp", "density.inp",
                        "dustdens.inp", "dusttemp_final.dat",
                        "dustopac.inp", "dustopac_*.inp", "abundance.inp",
                        "temperature.inp"]
        for iname in initfilelist:
            comm = subprocess.call(["cp", "./"+iname, rundir+"/"])


        ##Below Section: PROCESS data from the LTE or NLTE data file
        #Read in either LTE or NLTE data
        if self.lte: #If LTE treatment desired
            if verbose: #Verbal output, if so desired
                print("Extracting LTE molecular data...")
                print("")
            self.lineset = self._read_hitran(numprocessors)
        else: #Else, if non-LTE treatment desired
            if verbose: #Verbal output, if so desired
                print("Extracting NLTE molecular data...")
                print("")
            molfile = "molfile.dat"
            self.lineset = self._read_lambda(numprocessors)


        ##Below Section: DECIDE whether or not to generate NLTE files
        #!!!! - UNFINISHED - ASSUMING run_nlte=True FOR NOW
        #FILL THIS IN FROM LINE_RUN LINES 145-162
        run_nlte = True #!!!ASSUMING FOR NOW


        ##Below Section: RUN RADLITE on single/multiple processors in subfolders
        if verbose: #Verbal output, if so desired
            print("Running RADLite on "+str(numprocessors)+" processor(s)...")
        #Prepare pool of processors
        plist = []
        for ai in range(0, numprocessors):
            if verbose: #Verbal output, if so desired
                print("Prepping "+str(ai)+"th processor...")

            #Make a directory for this run
            cpudir = "./workingdir_cpu"+str(ai)+"/" #Processor-specific dir.
            if verbose: #Verbal output, if so desired
                print("Generating directory: "+cpudir)
            try:
                comm = subprocess.call(["mkdir", cpudir]) #Create subdirectory
            except (comm != 0): #If directory already exists, replace it
                comm = subprocess.call(["rm", "-r", cpudir]) #Erase prev. dir.
                comm = subprocess.call(["mkdir", cpudir]) #New empty dir.
            #Copy initial files into this new processor subdirectory
            for iname in initfilelist:
                comm = subprocess.call(["cp", rundir+"/"+iname, cpudir])

            #Call processor routine
            phere = mp.Process(target=self._run_processor,
                                    args=(pind, cpudir, rundir,))
            plist.append(phere)
            #Start process
            if verbose: #Verbal output, if so desired
                print("Starting "+str(ai)+"th processor in "+cpudir+"...")
            phere.start()

        #Close pool of processors
        if verbose: #Verbal output, if so desired
            print("Done running processor(s)!")
        for ai in range(0, numprocessors):
            if verbose: #Verbal output, if so desired
                print("Closing "+str(ai)+"th processor...")
            plist[ai].join()


        ##Below Section: FINISH and EXIT
        if verbose: #Verbal output, if so desired
            print("Done running run_lines()!")
        return
    #


    def _run_processor(self, pind, cpudir, rundir):
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
        ##NOTE: SOME OF THESE ARE PROBABLY GENERAL ENOUGH TO NOT BE DONE PER PROCESSOR
        self._write_radliteinp(cpudir) #Direct radlite input file
        self._write_velocityinp(cpudir) #Velocity input file
        self._write_gasdensityinp(cpudir) #Gas density input file
        self._write_abundanceinp(cpudir) #Abundance input file
        self._write_turbulenceinp(cpudir) #Turbulence input file
        self._write_levelpopinp(cpudir) #Level population input file
        self._write_moldatadat(filepathandname=cpudir, pind=pind,
            outfilename=(cpudir+"moldata_"+str(pind)+".dat")) #Mol. datafiles


        ##Below Section: CHECK that all files are in order
        #BELOW FROM IDL !!! - read_vel, 'velocity.inp', vel
        #Check passband width
        #velmax = np.max(np.abs(vel.vphi))
        #if verbose: #Verbal output, if so desired
        #    print("Max. velocity = {0:.2f}km/s".format(velmax/1E5))


        ##Below Section: RUN RADLITE
        logfile = open(cpudir+"RADLITE_core"+str(pind)+".log", 'w')
        comm = subprocess.call([self.executable], #Call RADLite
                                cwd=cpudir+"/", #Call within processor subdir.
                                stdout=logfile) #Send output to log file


        ##Below Section: EXTRACT output files
        comm = subprocess.call(["mv",
                            cpudir+"/moldata_"+str(pind)+".dat",
                            rundir+"/"]) #Molecular data used by this processor
        comm = subprocess.call(["mv",
                            cpudir+"/linespectrum_moldata_"+str(pind)+".dat",
                            rundir+"/"]) #Line spectrum output
        if self.image == 2:
            comm = subprocess.call(["mv",
                            cpudir+"/lineposvelcirc_moldata_"+str(pind)+".dat",
                            rundir+"/"]) #Circular 3D image cube output


        ##Below Section: DELETE subdir. + EXIT
        comm = subprocess.call(["rm", "-r", cpudir]) #Erase processor dir.
        return
    #


    ##READ METHODS
    def _read_hitran(self, filepathandname):
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
        if verbose: #Verbal output, if so desired
            print("Extracting molecular lines from file "+filepathandname+"...")
        #Read in all filelines
        with open(filepathandname, 'r') as openfile:
            alllines = openfile.readlines()
        if verbose: #Verbal output, if so desired
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
                                    for linehere in alllines, #All entries
                                    dtype=hitrandtypes[ai]]) #Entry datatype
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
        if verbose: #Verbal output, if so desired
            print("Removing molecular lines outside of specified criteria...")
        #Extract indices of lines that fall within criteria
        keepinds = ~np.where(
                    #If incorrect isotopologue
                    (hitrandict["isonum"] != self.isotop)
                    #If wavenumber outside of desired wavenumber range
                    | (hitrandict["wavenum"] < self.wavenumrange[0])
                    | (hitrandict["wavenum"] > self.wavenumrange[1])
                    #If intensity, level, or upper energy beyond given cutoffs
                    | (hitrandict["ins"] < self.cutoff)
                    | (hitrandict["Eu"] > self.Eupmax)
                    | (float(hitrandict["vup"]) > self.vupmax))[0]
        numlines = np.sum(keepinds) #Number of lines that fall within criteria

        #Keep only those lines that fall within criteria; delete all other lines
        for key in hitrandict:
            hitrandict[key] = hitrandict[key][keepinds]
        self.hitrandict = hitrandict #Record final set of molecular lines
        self.numlines = numlines #Record final count of lines
        if verbose: #Verbal output, if so desired
            print("There are "+str(numlines)+" molecular lines "
                    +"that fall within specified criteria.")


        ##Below Section: SPLIT data across given number of processors
        if verbose: #Verbal output, if so desired
            print("Dividing up the lines for "
                        +self.numprocessors+" processors...")
        #Determine indices for splitting up the data
        numpersplit = self.numlines // self.numprocessors #No remainder
        splitinds = [[(ai*numpersplit),((ai+1)*numpersplit)]
                        for ai in range(0, self.numprocessors)] #Divide indices
        splitinds[-1][1] = self.numlines #Tack leftovers onto last processor
        if verbose: #Verbal output, if so desired
            print("Here are the chosen line intervals per processor:")
            print([("Processor "+str(ehere[0])+": Interval "+str(ehere[1]))
                                    for ehere in enumerate(splitinds)])
        self.splitinds = splitinds #Record lines per processor
        ##Below Section: RETURN compiled lists of line data
        if verbose: #Verbal output, if so desired
            print("Done extracting molecular data!")
    #


    #












    ##READ METHODS
    BELOW STUFF NOT UPDATED YET!!!




    ##WRITE METHODS
    def _write_linespectrum_inp(self, numlines, molfilename, outfilename):
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


    def _write_molefile(self, outfilename, hitrandict):
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
        if verbose: #Verbal output, if so desired
            print("Writing "+outfilename+"...")
        ##Below Section: COMBINE levels + REMOVE duplicates to get unique levels
        if verbose: #Verbal output, if so desired
            print("Counting up unique levels for "+outfilename+"...")
        #Combine energies, transitions, and degeneracies
        Ealllist = np.concatenate((hitrandict["Elow"], hitrandict["Eup"]))
        valllist = np.concatenate((hitrandict["vlow"], hitrandict["vup"]))
        qalllist = np.concatenate((hitrandict["qlow"], hitrandict["qup"]))
        galllist = np.concatenate((hitrandict["glow"], hitrandict["gup"]))
        #Sort combined lists by energies
        sortinds = np.argsort(Ealllist)
        Ealllist = Eallist[sortinds]
        valllist = vallist[sortinds]
        qalllist = qallist[sortinds]
        galllist = gallist[sortinds]

        #Extract indices for unique levels only
        uniqEinds = ~np.where( #Indices of unique Elow, Eu value combinations
                    (np.abs(Ealllist[0:-1]-Ealllist[1:len(Ealllist)])
                            /1.0/(Ealllist[1:len(Ealllist)]+0.1)) < 0.0001)
        uniqvinds = ~np.where( #Indices of unique vlow, vup value combinations
                        vuplist[0:-1] == vuplist[1:len(Ealllist)])
        uniqginds = ~np.where( #Indices of unique glow, gup value combinations
                        glowlist[0:-1] == glowlist[1:len(Ealllist)])

        #Combine and apply indices
        if verbose: #Verbal output, if so desired
            print("Removing any duplicate levels for "+outfilename
                    +"...  To start, there are "+str(len(Ealllist))+" levels.")
        uniqinds = np.unique(np.concatenate((uniqEinds, uniqvinds, uniqginds)))
        Ealllist = Ealllist[uniqinds]
        valllist = valllist[uniqinds]
        qalllist = qalllist[uniqinds]
        galllist = galllist[uniqinds]
        if verbose: #Verbal output, if so desired
            print(str(np.sum(~uniqinds))+" levels have been removed for "+
            print(outfilename+", leaving "+str(len(Eallist))+" unique levels.")


        ##Below Section: BUILD string to form the molecular data file
        writestr = ""
        writestr += "!MOLECULE\n{0:s}\n".format(self.molname)
        writestr += "!MOLECULAR WEIGHT\n{0:4.1f}\n".format(self.molweight))
        writestr += "!NUMBER OF ENERGY LEVELS\n{0:6d}\n".format(
                            len(Ealllist))
        #Tack on unique levels, energies, and degeneracies
        writestr += "!LEVEL + ENERGIES(cm^-1) + WEIGHT + v + Q\n"
        for ai in range(0, len(Ealllist)):
            writestr += "{0:5d}{1:12.4f}{2:7.1f}{3:>15s}{4:>15s}\n".format(
                            (ai+1), Ealllist[ai], galllist[ai], valllist[ai],
                            qalllist[ai])
        #Tack on transitions
        writestr += "!NUMBER OF RADIATIVE TRANSITIONS\n"
        writestr += "{0:6d\n}".format(len(hitrandict["wavenum"]))
        writestr += ("!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(cm^-1) + "
                        +"E_u(cm^-1) + v_l + Q_p + Q_pp\n")
        for ai in range(0, len(hitrandict["wavenum"])):
            levu = np.where((np.abs(Ealllist - hitrandict["Eup"][ai])
                                /1.0/hitrandict["Eup"][ai]) < 1E-4)[0]
            if hitrandict["Elow"][ai] != 0: #If not down to 0-level
                levl = np.where(np.abs(Ealllist - hitrandict["Elow"][ai])
                                /1.0/hitrandict["Elow"][ai]) < 1E-4)[0]
            else:
                levl = np.where(Ealllist == 0)[0]
            writestr += ("{0:5d}{1:5d}{2:5d}{3:12.3e}{4:16.7f}".format(
                                (ai+1), levu[0]+1, levl[0]+1,
                                hitrandict["A"][ai], hitrandict["wavenum"][ai])
                        +"{5:12.5f}{6:>15}{7:>15}{8:>15}{9:>15}\n".format(
                                hitrandict["Eup"][ai], hitrandict["vup"][ai],
                                hitrandict["vlow"][ai], hitrandict["qup"][ai],
                                hitrandict["qlow"][ai]))


        ##Below Section: WRITE the results to file + EXIT function
        with openfile as open(outputfilename, 'w'):
            f.write(writestr)
        return
    #


        ##



    ##WRITE METHODS
    _write_radliteinp(cpudir) #Direct radlite input file
    _write_velocityinp(cpudir) #Velocity input file
    _write_gasdensityinp(cpudir) #Gas density input file
    _write_abundanceinp(cpudir) #Abundance input file
    _write_turbulenceinp(cpudir) #Turbulence input file
    _write_levelpopinp(cpudir) #Level population input file
    _write_moldatadat(filepathandname=cpudir,




#
#









#
