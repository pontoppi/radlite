import os
from itertools import islice
import numpy as np
import matplotlib.pylab as plt
import astropy.io.ascii as at

from scipy.integrate import simps
import scipy.constants as cst

from scipy.integrate import cumtrapz

class radmc_model():
    def __init__(self,path):
        self.path = path
        self.radius = self.read_radius()
        self.theta = self.read_theta()
        self.frequency = self.read_freq()
        self.wavelength = cst.c*1e6/self.frequency # micron
        self.gasdensity = self.read_gasdensity()
        self.dustdensity = self.read_dustdensity()
        self.dusttemperature = self.read_dusttemperature()
        self.dustopac = self.read_opacity()
        self.pars = self.read_parameters()

        self.nr = self.radius.shape[0]
        self.nt = self.theta.shape[0]
        self.nf = self.frequency[0]
        self.nspecies = self.dustdensity.shape[2]

        self.x,self.y = self.cartesian()

    def cartesian(self):
        x = np.zeros((self.nr,self.nt))
        y = np.zeros((self.nr,self.nt))
        for i in np.arange(self.nr):
            x[i,:] = self.radius[i]*np.sin(self.theta)
            y[i,:] = self.radius[i]*np.cos(self.theta)
        return x, y

    def gasmass(self):
        cols = np.zeros(self.nr)
        for i in np.arange(self.nr):
            cols[i] = simps(self.gasdensity[i,:],self.theta*self.radius[i])

        mass = simps(2*np.pi*self.radius*cols,self.radius)
        return mass

    def dustmass(self):
        masses = []
        for j in np.arange(nspecies):
            cols = np.zeros(self.nr)
            for i in np.arange(self.nr):
                cols[i] = simps(self.dustdensity[i,:,j],self.theta*self.radius[i])

            mass = simps(2*np.pi*self.radius*cols,self.radius)

        masses.append(mass)

        return masses

    def read_radius(self):
        filename = 'radius.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        nr = int(file.readline())

        radius = np.zeros(nr)
        i = 0
        for l in file:
            if not l.isspace():
                radius[i] = float(l)
                i += 1

        file.close()

        return radius

    def read_theta(self):
        filename = 'theta.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        nt,mirror = file.readline().split()
        nt = int(nt)
        mirror = int(mirror)

        theta = np.zeros(nt)
        i = 0
        for l in file:
            if not l.isspace():
                theta[i] = float(l)
                i += 1

        if mirror==1:
            theta = np.append(theta,np.pi-theta[nt-1::-1])

        file.close()

        return theta

    def read_freq(self):
        filename = 'frequency.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        nnu = int(file.readline())

        freq = np.zeros(nnu)
        i = 0
        for l in file:
            if not l.isspace():
                freq[i] = float(l)
                i += 1

        file.close()

        return freq

    def read_opacity(self):
        filename = 'dustopac.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')
        file_format = int(file.readline().split()[0])
        nspecies = int(file.readline().split()[0])
        if nspecies>1:
            raise ValueError('Currently only models with one dust species are supported for this method')
        file.close()

        filename = 'dustopac_1.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')
        nfreq, ispecies = file.readline().split()
        nfreq = int(nfreq)

        dustopac = {'cabs':np.zeros(int(nfreq)),'csca':np.zeros(int(nfreq))}

        empty = file.readline()

        lines = list(islice(file, nfreq))
        for i,l in enumerate(lines):
            dustopac['cabs'][i] = float(l)

        empty = file.readline()

        lines = list(islice(file, nfreq))
        for i,l in enumerate(lines):
            dustopac['csca'][i] = float(l)

        file.close()
        return dustopac

    def read_gasdensity(self):
        filename = 'density.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        nr,nt,mirror = file.readline().split()
        nr = int(nr)
        nt = int(nt)
        mirror = int(mirror)

        density = np.zeros((nr,nt*(1+mirror)))
        i = 0
        j = 0
        for l in file:
            if not l.isspace():
                density[j,i] = float(l)
                if i==nt-1:
                    density[j,i+1:] = density[j,nt-1::-1]
                    i = 0
                    j += 1
                else:
                    i += 1

        file.close()

        return density

    def read_dustdensity(self):
        #multiple densities not yet supported
        filename = 'dustdens.inp'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        nspecies,nr,nt,mirror = file.readline().split()
        nspecies = int(nspecies)
        nr = int(nr)
        nt = int(nt)
        mirror = int(mirror)

        density = np.zeros((nr,nt*(1+mirror),nspecies))
        i = 0
        j = 0
        k = 0
        for l in file:
            if not l.isspace():
                density[j,i,k] = float(l)

                if j==nr:
                    i = 0
                    j = 0
                    k += 1
                else:
                    if i==nt-1:
                        density[j,i+1:,k] = density[j,nt-1::-1,k]
                        i = 0
                        j += 1
                    else:
                        i += 1


        file.close()

        return density

    def read_dusttemperature(self):
        #multiple densities not yet supported
        filename = 'dusttemp_final.dat'
        fullpath = os.path.join(self.path,filename)
        file = open(fullpath,'r')

        nspecies,nr,nt,mirror = file.readline().split()
        empty = file.readline()
        idum = file.readline()

        nspecies = int(nspecies)
        nr = int(nr)
        nt = int(nt)
        mirror = int(mirror)

        temperature = np.zeros((nr,nt*(1+mirror),nspecies))
        i = 0
        j = 0
        k = 0
        for l in file:
            if not l.isspace():
                temperature[j,i,k] = float(l)

                if j==nr:
                    i = 0
                    j = 0
                    k += 1
                else:
                    if i==nt-1:
                        temperature[j,i+1:,k] = temperature[j,nt-1::-1,k]
                        i = 0
                        j += 1
                    else:
                        i += 1

        file.close()

        return temperature

    def read_parameters(self):
        filename = 'problem_params.pro'
        fullpath = os.path.join(self.path,filename)
        f = open(fullpath)
        content = [x.strip('\n') for x in f.readlines() if '=' in x]

        pars = {}
        for line in content:
            elements = line.split('=')
            try:
                key = elements[0].strip()
                value = float(elements[1])
            except:
                key = elements[0].strip()
                value = elements[1]
            pars[key] = value

        return pars

    def surface_at_tau(self,tau,wave=0.55,au=False):

        freq = (cst.c*1e6)/wave
        cabs_at_wave = np.interp(freq,self.frequency,self.dustopac['cabs'])
        csca_at_wave = np.interp(freq,self.frequency,self.dustopac['csca'])
        opac_at_wave = cabs_at_wave+csca_at_wave

        cum = cumtrapz(self.dustdensity[:,:,0],-self.y,initial=0.)*opac_at_wave  # only using first dust species
        ys = np.zeros(self.nr)
        for i in np.arange(self.nr):
            csub = np.argmin(np.abs(cum[i,:].flatten()-tau))
            ys[i] = self.y[i,csub]

        if au:
            return (self.radius/self.au,ys/self.au)
        else:
            return (self.radius,ys)
