#!/usr/bin/env python

# EEEEEE   XX    XX         AAA        TTTTTTTTTT
# EE        XX  XX         AA AA       TT  TT  TT
# EE         XXXX         AA   AA          TT
# EEEEEE      XX         AAA   AAA         TT
# EE         XXXX       AAAAAAAAAAA        TT
# EE        XX  XX     AA         AA       TT
# EEEEEE   XX    XX   AA           AA      TT   
#
# EXcitonic Analysis Tool         @MoLECoLab 
# https://molecolab.dcci.unipi.it/tools/
#

#
# *************************************
# EXAT - EXcitonic Analysis Tool
# specview.py VISUALIZE SPECTRA MODULE
# *************************************
#
# Copyright (C) 2014-2017 
#   S. Jurinovich, L. Cupellini, C.A. Guido, and B. Mennucci
#
# This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE or at
#   <http://www.gnu.org/licenses/>.
#

from __future__ import print_function,division
import sys,os
import matplotlib.pyplot as plt
import numpy    as np
import common   as c
import argparse as arg


# *****************************************************************************
#
# Funtion to visualize the spectrum 
#
# x     ... X-Axis (value in cm-1)
# E     ... Transition energies (cm-1)
# OD    ... Y-Axis data of Linear Apsortion
# LD    ... Y-Axis data of Linear Dichroism
# CD    ... Y-Axis data of Circular Dichroism
# mu2od ... Square of transition dipole moments (Debye**2) for OD
# mu2ld ... Square of transition dipole moments (Debye**2) for LD
# rotst ... Rotational strength (10^-40 c.g.s.)
#
#
def view(OPT,x,E,OD,LD,CD,mu2od,mu2ld,rotst):
  ODsticks  = mu2od*E # intensity = mu^2*energy
  CDsticks  = rotst*E # intensity = R*energy
  
  Plots = ['uv','cd']
  NPlot = len(Plots)

  # Set subplots
  f,ec = plt.subplots(NPlot,sharex=True)
  plt.subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=None, hspace=None)

  # Define X-Axis units
  if   OPT['UAxi'] == "nm" :
    print(" Spectrum will be shown in nm")
    x = 1.0E7/x
    E = 1.0E7/E
    plt.xlim(1.0E7/OPT['XMax'],1.0E7/OPT['XMin'])
    E1 = np.where(E>plt.xlim()[0])[0]
    E2 = np.where(E<plt.xlim()[1])[0]
    IndEne = np.intersect1d(E1,E2)

  elif OPT['UAxi'] == "cm-1" :
    print(" Spectrum will be shown in cm-1")
    plt.xlim(OPT['XMin'],OPT['XMax'])
    E1 = np.where(E>plt.xlim()[0])[0]
    E2 = np.where(E<plt.xlim()[1])[0]
    IndEne = np.intersect1d(E1,E2)

  
  elif OPT['UAxi'] == "eV" :
    print(" Spectrum will be shown in eV")
    x = x/c.PhyCon['eV2wn']
    E = E/c.PhyCon['eV2wn']
    plt.xlim(OPT['XMin']/c.PhyCon['eV2wn'],OPT['XMax']/c.PhyCon['eV2wn'])
    E1 = np.where(E>plt.xlim()[0])[0]
    E2 = np.where(E<plt.xlim()[1])[0]
    IndEne = np.intersect1d(E1,E2)


  k = 0

  # Plot UV
  if 'uv' in Plots:
    TOD = np.sum(OD,axis=1)
    MOD = np.max(TOD)
    # Normalize sticks
    ODsticks *= MOD/np.max(ODsticks[IndEne])
    ec[k].set_ylim(0,MOD+MOD*0.1)
    ec[k].set_title('Absorption Spectrum')
    ec[k].plot(x,TOD,linewidth=2.5,linestyle="-",color="red")
    ec[k].vlines(E,ODsticks,0,linewidth=2.0)  
    k+=1

  # Plot LD
  if 'ld' in Plots:
    TLD = np.sum(LD,axis=1)
    MLD = np.max(TLD)
    mLD = np.min(TLD)
    ec[k].set_ylim(mLD+mLD*0.1,MLD+MLD*0.1)
    ec[k].set_title('Linear Dichroism Spectrum')
    ec[k].plot(x,TLD,linewidth=2.5,linestyle="-",color="red")
    ec[k].plot(x,np.zeros(len(x)),linewidth=1.0,linestyle="-",color="black")
    ec[k].vlines(E,mu2ld*MLD/np.max(np.abs(mu2ld[IndEne])),0)  
    k+=1
 
  # Plot CD
  if 'cd' in Plots:
    TCD = np.sum(CD,axis=1)
    MCD = np.max(TCD)
    mCD = np.min(TCD)
    # Normalize sticks
    CDmax = max(abs(mCD),MCD)
    STmax = max(np.max(CDsticks[IndEne]),abs(np.min(CDsticks[IndEne])))
    CDsticks *= CDmax/STmax

    ec[k].set_ylim(mCD+mCD*0.1,MCD+MCD*0.1)
    ec[k].set_title('Circular Dichroism Spectrum')
    ec[k].plot(x,TCD,linewidth=2.5,linestyle="-",color="red")
    ec[k].plot(x,np.zeros(len(x)),linewidth=1.0,linestyle="-",color="black")
    ec[k].vlines(E,CDsticks,0,linewidth=2.0)  
    k+=1
  
  # Show plot
  plt.show()


# *****************************************************************************

#
# If called from command line
#
if __name__ == "__main__":

  # Set default options
  OPT = { 'pltformat' : 'png',
          'pltname'   : 'spectrum'}

  # Read input line
  parser = arg.ArgumentParser(description="Plot spectra using MatPlotLib library")
  parser.add_argument('-i',help='Data for stiks',default="ec.dat")
  parser.add_argument('-uv',help='Absorption data file',default="spec.uv.dat")
  parser.add_argument('-cd',help='CD data file',default="spec.cd.dat")
  parser.add_argument('--max',help='High energy limit of the spectrum (nm)', type=float,default=200)
  parser.add_argument('--min',help='Low energy limit of the spectrum (nm)', type=float,default=800)
  parser.add_argument('-save',help='File name for the picture to save',default="spectrum")
  parser.add_argument('--format',help='Format for the figure',choices=['eps','ps','png','pdf'],default='png')

  # Version
  parser.add_argument('--version','-V',action='version',version=c.PROGVERS,
                      help="Show program's full version number, and exit")  

  args = parser.parse_args()

  # Set options
  InFile = args.i
  UVFile = args.uv
  CDFile = args.cd
  OPT['Emax'] = args.max
  OPT['Emin'] = args.min
  OPT['pltname'] = args.save
  OPT['pltformat'] = args.format

  # Print the welcome message
  c.welcome()

  # Read the input files
  c.checkfile(InFile)
  c.checkfile(UVFile)
  c.checkfile(CDFile)

  StData = np.loadtxt(InFile)
  Energy = 1E7/(StData[:,1]*c.PhyCon['eV2wn'])
  Dipo2  = StData[:,2]
  RotSt  = StData[:,3]

  UVData = np.loadtxt(UVFile)
  CDData = np.loadtxt(CDFile)

  # View the spectrum
  view(OPT,UVData,CDData,Energy,Dipo2,RotSt)

  # END
  print("\n Done! \n")

