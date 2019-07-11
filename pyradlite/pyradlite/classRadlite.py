##FILE:
##PURPOSE:


##
class Radlite():
    def __init__():
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


        ##Below Section: CHECK input parameters for any user error
        ##Below Section: RECORD inputs as attributes



    def make_abundance():

    def make_velocity():

    def make_gastemperature():

    def make_moldata():

    def make_levelpop():

    def run_lines(self):
        """
        DOCSTRING
        Function:
        Purpose:
        Inputs:
        Variables:
        Outputs:
        Notes:
        """
        ##Below Section: SET UP result directory
        printstamp = time.datetoprint?() #???
        printtime = time.timetoprint?() #???
        if nodate: #If date should not be appended to directory name
            rundir = run_name + printstamp
        else: #If date should be appended to directory name
            rundir = run_name + printstamp + printtime
        #Make the desired directory, if nonexistent
        ???
        if verbose: #Verbal output, if so desired
            print("All data will be saved to the following "
                    +"directory: "+rundir)


        ##Below Section: CHOOSE RADLite exec. based on desired output
        #Prepare executable for desired image cube output
        if self.image == 0: #If desired cube output is spectrum
            if verbose: #Verbal output, if so desired
                print("Preparing a spectrum-formatted image cube...")
                print("")
            executable = exe_path+"RADlite"
        elif self.image == 2: #Else, if desired cube output is circular
            if verbose: #Verbal output, if so desired
                print("Preparing a circular-formatted image cube...")
                print("")
            executable = exe_path+"RADlite_imcir"


        ##Below Section: PROCESS data from the LTE or NLTE data file
        #Read in either LTE or NLTE data from other file
        if self.lte == 1: #If LTE treatment desired
            if verbose: #Verbal output, if so desired
                print("Extracting LTE molecular data...")
                print("")
            linedict = self._read_hitran()
        elif self.lte == 0: #Else, if NLTE treatment desired
            if verbose: #Verbal output, if so desired
                print("Extracting NLTE molecular data...")
                print("")
            molfile = "molfile.dat"
            linedict = self._read_lambda()


        ##Below Section: DECIDE whether or not to generate NLTE files
        #!!!! - UNFINISHED - ASSUMING run_nlte=True FOR NOW
        #FILL THIS IN FROM LINE_RUN LINES 145-162
        run_nlte = True #!!!ASSUMING FOR NOW


        ##Below Section:
        #RUN - SET_MODEL (IF NOT ALREADY DONE SO (I.E., SELF.MODEL=NONE))
        #RUN - CALC_LINES
        #RUN - RADLITE
        #SAVE+MOVE RESULTS





        #Below assembles/converts line wavelengths
        #wave = np.sort(1E4/1.0/freq) #Sorted line wavelengths [um]
        #nlines = len(wave) #Number of lines


        ##Below Section:
        #Below ensures that

    #


    ##READ METHODS
    def _read_hitran(self, filepathandname, numprocessors):
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
        ##Below Section: ITERATE through lines in file
        if verbose: #Verbal output, if so desired
            print("Extracting molecular lines from file "+filepathandname+".")
        #Set up lists to hold data
        mollist = []
        isolist = []
        wavenumlist = []
        Alist = []
        Elowlist = []
        Euplist = []
        vuplist = []
        vlowlist = []
        guplist = []
        glowlist = []
        #Iterate through each line
        openfile = open(filepathandname, 'r')
        for linehere in openfile:
            #Extract molecular info from current line (160 tot. char.)
            molinehere = int(linehere[0:1+1]) #2 char.; molecule
            isohere = int(linehere[2:2+1]) #1 char.; isotopologue
            wavenumhere = float(linehere[3:14+1]) #12 char.; wavenumber [cm^-1]
            inshere = float(linehere[15:24+1]) #10 char.; line intensity
            Ahere = float(linehere[25:34+1]) #10 char.; Einstein A coeff.
            airhere = float(linehere[35:39+1]) #5 char. #!!!FIX/FINISH COMMENTS
            selfhere = float(linehere[40:44+1]) #5 char.
            Elowhere = float(linehere[45:54+1]) #10 char.
            Euphere = Elowhere + wavehere
            temphere = float(linehere[55:58+1]) #4 char.
            preshere = float(linehere[59:66+1]) #8 char.
            vuptext = linehere[67:81+1] #15 char.
            vlowtext = linehere[82:96+1] #15 char.
            quptext = linehere[97:111+1] #15 char.
            qlowtext = linehere[112:126+1] #15 char.
            errhere = int(linehere[127:132+1]) #6 char.
            refhere = int(linehere[133:144+1]) #12 char.
            flaghere = linehere[145:145+1] #1 char.
            guphere = float(linehere[146:152+1]) #7 char.
            glowhere = float(linehere[153:159+1]) #7 char.
pandas, astropy - lookinto; read into pandas dataframe
            #Throw an error if line contains signs of incomplete data
            if guphere == 0:
                raise ValueError("Something's wrong!  Incomplete molecular "
                        +"HITRAN data at wavenum="+str(wavenumhere)+".")

            #Skip this line if falls outside of desired range
            if ((isonum != self.isotop) #If incorrect isotopologue
                    or (wavenumhere < self.wavenumrange[0]) #If < wavenum range
                    or (wavenumhere > self.wavenumrange[1]) #If > wavenum range
                    or (inshere < self.cutoff) #If intensity < cutoff
                    or (Euphere > self.Eupmax) #If E_up > E_up cutoff
                    or (float(vuphere) > self.vupmax): #If v_up > v_up cutoff
                continue

            #Record the parsed data
            mollist.append(molinehere)
            isolist.append(isohere)
            wavenumlist.append(wavenumhere)
            Alist.append(Ahere)
            Elowlist.append(Elowhere)
            Euplist.append(Euphere)
            vuplist.append(vuphere)
            vlowlist.append(vlowhere)
            quplist.append(quphere)
            qlowlist.append(qlowhere)
            guplist.append(guphere)
            glowlist.append(glowhere)
        #Politely close the file
        openfile.close()

        #Convert the lines into sorted arrays (for easier parsing later)
        sortinds = np.argsort(Elowlist) #Indices, sorted by E_low
        mollist = np.asarray(mollist)[sortinds]
        isolist = np.asarray(isolist)[sortinds]
        wavenumlist = np.asarray(wavenumlist)[sortinds]
        Alist = np.asarray(Alist)[sortinds]
        Elowlist = np.asarray(Elowlist)[sortinds]
        Euplist = np.asarray(Euplist)[sortinds]
        vlowlist = np.asarray(vlowlist)[sortinds]
        vuplist = np.asarray(vuplist)[sortinds]
        qlowlist = np.asarray(qlowlist)[sortinds]
        quplist = np.asarray(quplist)[sortinds]
        glowlist = np.asarray(glowlist)[sortinds]
        guplist = np.asarray(guplist)[sortinds]
        numlines = len(wavenumlist)
        if verbose: #Verbal output, if so desired
            print("There are "+str(numlines)+" molecular lines.")


        ##Below Section: SPLIT data across the given number of processors
        if verbose: #Verbal output, if so desired
            print("Dividing up the lines for "+numprocessors+" processors...")
        #Determine indices for splitting up the data
        numpersplit = numlines // numprocessors #Number of lines, no remainder
        splitinds = [[(ai*numpersplit),((ai+1)*numpersplit)]
                            for ai in range(0, numprocessors)] #Dividing indices
        splitinds[-1][1] = numlines #Tack leftovers onto last processor
        #Divide up data using indices
        mollist_split = [mollist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        isolist_split = [isolist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        Alist_split = [Alist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        wavenumlist_split = [wavenumlist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        Elowlist_split = [Elowlist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        Euplist_split = [Euplist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        vuplist_split = [vuplist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        vlowlist_split = [vlowlist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        guplist_split = [guplist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        glowlist_split = [glowlist[splitinds[ai][0]:splitinds[ai][1]]
                                for ai in range(0, numprocessors)]
        if verbose: #Verbal output, if so desired
            print("Here are the chosen line intervals per processor:")
            print([("Processor "+str(ehere[0])+": Interval "+str(ehere[1]))
                                    for ehere in enumerate(splitinds)])


        ##Below Section: RETURN compiled lists of line data
        if verbose: #Verbal output, if so desired
            print("Done extracting molecular data!")
        return {"mol":mollist_split, "iso":isolist_split, "A":Alist_split,
                "wavenum":wavenumlist_split, "Elow":Elowlist_split,
                "Eup":Euplist_split, "vup":vuplist_split,
                "vlow":vlowlist_split, "gup":guplist_split,
                "glow":glowlist_split,
                "Eall_uniq":Ealllist, "vall_uniq":valllist,
                "qall_uniq":qalllist, "gall_uniq":galllist} - attri.
    #


    def _read_model(self, modelpathandname):
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
        ##Below Section: LOAD the model data
        self.model = np.load(modelpathandname).item()
    #


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




#
#









#
