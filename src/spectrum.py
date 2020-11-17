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
# spectrum.py COMPUTE SPECTRA
# *************************************
#
#
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

import sys,os
import argparse as arg
import math     as m
import numpy    as np
import common   as c
from scipy import signal

OPT = {}

#
# =============================================================================
# BUILD THE EXCITONIC SPECTRA:
# =============================================================================
#
# INPUT :  energy  : site enrgies (cm-1)                       [array]
#          dipoles : square of electric tr. dip. moment (a.u.) [array]
#          rotstr  : rotational strength (10^-40 esu2cm2)      [array]
#          sigma   : bandwith of gaussian lineshape (eV)       [array]
#          OPT     : options (see help)                        [dictionary]
#
# OUTPUT : UV      : Absorption spectrum (eps vs nm)           [array]
#          CD      : Circular Dichroism spectrum (Deps vs nm)  [array]
#

#   The Input data files must contains 5 coloumns:  
#   (1) Id number of the transition;   
#   (2) Energy of the transition (ev);     
#   (3) Square of dipole moment (a.u.);      
#   (4) Linear Dichro;
#   (5) Rotational strenght (10^-40 esu2cm2...).
#   An additional coloumn can be specified for bandwidth specification.""")
#


# *****************************************************************************
factOD = 108.86039     # w*|mu|^2 to epsilon         if mu^2 is in Debye = 10^-36 esu^2 cm^2
factCD = factOD*4e-4   # w*R      to Delta(epsilon)  if R    is in 10^-40 esu^2 cm^2

def specalc(NTran,energy,dipo,dipold,rotstr,broad):

  # Convert square dipoles from a.u. to Debeye
  dipo *= c.PhyCon['ToDeb']**2

  # Define X-Axis (in wavenumbers)
  x=np.arange(OPT['XMin'],OPT['XMax'],OPT['XStp'],dtype=float)
  NStep = len(x) 
 
  # Initilize the intensity
  OD = np.zeros((NStep,NTran),dtype=float)
  CD = np.zeros((NStep,NTran),dtype=float)
  LD = np.zeros((NStep,NTran),dtype=float)

  # ILLEGGIBILE! (controllare formule)
  if OPT['LShp'] == "gaussian":
    LineShape = lambda p, x: (x/(p[1]*np.sqrt(2*m.pi))*np.exp(-0.5*((x-p[0])/p[1])**2)) 

  elif OPT['LShp'] == "lorentzian":
    LineShape = lambda p, x: (x/(np.pi*p[1]))*(p[1]**2/((x-p[0])**2+p[1]**2))

  for j in range(NTran):	
    p = [energy[j],broad[j]]
    OD[:,j] = LineShape(p,x)*dipo[j]*factOD
    LD[:,j] = LineShape(p,x)*dipold[j]*factOD
    CD[:,j] = LineShape(p,x)*rotstr[j]*factCD 
 
  return(x,OD,LD,CD)

# *****************************************************************************
#
# Compute maximum and min of spectra
#

def findpeaks(typespec,x,y,broad):


#broad = np.mean(broad)
#width = np.arange(1,10)*broad/3
#peakind = signal.find_peaks_cwt(y, width)
 peakind = signal.argrelmax(y)[0]
 if typespec is not 'OD': peakind = np.append(peakind,signal.argrelmin(y)[0])
 peakind = np.sort(peakind)

 NPeaks = len(peakind)
 print "\n Found %2d peaks in %s spectrum%s" % (NPeaks,typespec,":"*np.sign(NPeaks))
 for i in range(NPeaks):
   Ind = peakind[i]
   print " > %2d  X =   %8.0f cm-1    %6.1f nm   %8.4f eV   Y= % 8.4E" % (i+1,x[Ind],1.0E7/x[Ind],x[Ind]/c.PhyCon['eV2wn'],y[Ind])

 return


# *****************************************************************************
#
# Compute auto xmin and xmax 
#
def autowindow(EMin,EMax,broad=None): 
  if broad is not None:
    OPT['XMin'] = EMin-4*broad[0]
    OPT['XMax'] = EMax+4*broad[-1]
  else:
    OPT['XMin'] = EMin-1000.0
    OPT['XMax'] = EMax+1000.0


# *****************************************************************************
#
# Compute auto xmin and xmax 
#
def printlimits():
  if OPT['UAxi'] == "nm"     : 
    print(" Spectrum Limits: <-------------- lambda (nm) -----------") 
    print("                    %4.0f                           %4.0f" % (1.E7/OPT['XMax'],1.E7/OPT['XMin']) ) 
  elif OPT['UAxi'] == "cm-1" : 
    print(" Spectrum Limits: ----------- wavenumber (cm-1) -------->") 
    print("                  %8.0f                       %8.0f"       % (OPT['XMin'],OPT['XMax']) )
  elif OPT['UAxi'] == "eV" : 
    print(" Spectrum Limits: ---------------- E (eV) -------------->") 
    print("                  %8.4f                       %8.4f"       % (OPT['XMin']/c.PhyCon['eV2wn'],OPT['XMax']/c.PhyCon['eV2wn']) )



# *****************************************************************************
#
# Save Output Spectra
#

def savespec(typespec,x_in,y_in,OPT):
  x = x_in.copy()
  y = y_in.copy()

  # Set to zero very small numbers
  y[np.abs(y)<1e-99] = 0.0

  if OPT['bin']:
    OutFileName = "%s.%s.npy" % (OPT['FOut'],typespec)
  else:
    OutFileName = "%s.%s.dat" % (OPT['FOut'],typespec)
  
    OutFile = open(OutFileName,'w')
    OutFile.write("# %s spectrum generated with spectrum.py\n" % typespec)
    OutFile.write("# VERSION: %s\n" % c.VERSION)
    OutFile.write("# Initial command: \n")
    OutFile.write("# %s \n" % ' '.join(sys.argv))
  
  N = 1

  # Coumpute the total spectrum
  sum_y = np.sum(y,axis=1)


  # If required normalize the spectrum at the maximum (positive) value
  if OPT['INrm'] == True :
    max_y = np.max(sum_y)
    if typespec in ('CD','LD'):
      max_y = max(max_y,abs(np.min(sum_y)))
    if max_y > 0.0:
      y     /= max_y
      sum_y /= max_y

  # If requested save the contributions
  if OPT['ICnt'] == True : 
    sum_y = np.column_stack((sum_y, y))
    N += np.shape(y)[1]

  # Add the x coulumn to the output matrix
  if OPT['UAxi'] == "nm":
    outfmt = '%12.4f'+' %18.10E'*N
    x = 1.0E7/x
  elif OPT['UAxi'] == "eV":
    outfmt = '%12.6f'+' %18.10E'*N
    x /= c.PhyCon['eV2wn']
  elif OPT['UAxi'] == "cm-1":
    outfmt = '%12.2f'+' %18.10E'*N

  spec  = np.column_stack((x, sum_y))

  if OPT['bin']:
    np.save(OutFileName,spec)
  else:
    np.savetxt(OutFile,spec,fmt=outfmt,newline="\n")
    OutFile.close()
  pass



# *****************************************************************************
#
# If called from command line
#
if __name__ == "__main__":

  # Read input line
  parser = arg.ArgumentParser(prog="spectrum.py",description=""" 

    This module plot the absorption OD LD and CD spectrum starting from the
    the transition energies and the relative intensities. A Gaussian or Lorentzian
    function can be used to obtain the spectra.""")

  parser.add_argument('input',   help='Input file (default = results.out)',default="results.out",nargs='?')
  parser.add_argument('-o','--out',  help='Prefix of the output file containing the data' ,default="spec")
  parser.add_argument('-p','--plot', help='Plot the spectra on screen using matlibplot library',action="store_true")
  parser.add_argument('-n','--norm', help='Normalize the spectra to the maximum',action="store_true")
  parser.add_argument('--unitaxis',  help='Define the unit for the limit of the spectrum',choices=["nm","eV","cm-1"],default="nm")
  parser.add_argument('--unitbroad', help='Define the unit for lineshape broadening',choices=["eV","cm-1"],default="cm-1")
  parser.add_argument('--max',   help='High energy limit',default=180.0, type=float)
  parser.add_argument('--min',   help='Low energy limit',default=900.0,type=float)
  parser.add_argument('--autorange', help='Set automatic energy range',action="store_true")
  parser.add_argument('--step',  help='Frequency step (cm-1)',default=2.0,type=float)
  parser.add_argument('--bin', help='Save binary npy files instead of text files',action="store_true")
  parser.add_argument('--shape', help='Define the spectral lineshape',choices=["gaussian","lorentzian"],default="gaussian")
  parser.add_argument('--lorentzian',help='Equivalent to --shape lorentzian',action="store_true")
  parser.add_argument('--contr', help='Save the contributions of single transitions to the total spectrum',action="store_true")
  parser.add_argument('-s','--shift', help='Shift the spectrum of given quantity',type=float)

  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--HWHH','--gamma','--broadening', help="Half Width at Half Heigth for any lineshape, if HWHH=read: read the values from input file")
  group.add_argument('--sigma', help="sigma for the Gaussian lineshape,if sigma=read : read the values from input file")

  # Version
  parser.add_argument('--version','-V',action='version',version=c.PROGVERS,
                      help="Show program's full version number, and exit")  

  args = parser.parse_args()

  InFileName   = args.input
  OPT['FOut']  = args.out
  OPT['XMax']  = args.max
  OPT['XMin']  = args.min
  OPT['XStp']  = args.step
  OPT['LShp']  = args.shape
  OPT['IPlt']  = args.plot
  OPT['INrm']  = args.norm
  OPT['ICnt']  = args.contr
  OPT['UBrd']  = args.unitbroad
  OPT['UAxi']  = args.unitaxis
  OPT['shift'] = args.shift
  OPT['bin']   = args.bin
  autorange    = args.autorange
  if autorange is True:
    autorange = True
    OPT['XMax'] = 'Auto'
    OPT['XMin'] = 'Auto'

  if args.lorentzian:
    OPT['LShp'] = 'lorentzian'

  # Parse and process broeadening type
  OPT['BRed'] = False
  OPT['BVal'] = 0.0
  if args.HWHH  != None : 
    OPT['BTyp'] = "HWHH"
    if args.HWHH == "read"  : OPT['BRed'] = True
    else: OPT['BVal'] = float(args.HWHH)
  if args.sigma != None : 
    OPT['BTyp'] = "sigma"
    if args.sigma == "read" : OPT['BRed'] = True
    else: OPT['BVal'] = float(args.sigma)

  #print OPT['BTyp']
  #print OPT['BRed']
  #print OPT['BVal']

  # Print the welcome message
  c.welcome()

  # Check if the input file exists
  c.checkfile(InFileName)


  # Print Initial Information
  print(" XAxis will be saved in %s" % OPT['UAxi'])

  # Print user-defined limits and chec and check
  if not autorange: 
    # Convert UAxi in wavenumbers/eV
    if OPT['UAxi'] == "nm" :
      OPT['XMax'] = 1.0E+7/OPT['XMax']
      OPT['XMin'] = 1.0E+7/OPT['XMin']
    elif OPT['UAxi'] == "eV" :
      OPT['XMax'] *= c.PhyCon['eV2wn']
      OPT['XMin'] *= c.PhyCon['eV2wn']
  
    # Check extrema
    if OPT['XMin']>OPT['XMax']:
      # Swap min and max
      tmp = OPT['XMin']
      OPT['XMin'] = OPT['XMax']
      OPT['XMax'] = tmp
    elif OPT['XMin'] == OPT['XMax'] :
      c.error("XMin is equal to XMax!","XAxis")


  data = np.loadtxt(InFileName,comments="#").T      # Read data from input file
  # ^ Transpose to tweak case w. only one line
  energy = np.atleast_1d(data[1])*c.PhyCon['eV2wn']    # Convert the energy coulomn from eV to cm-1
  dipo   = np.atleast_1d(data[2])                    
  dipold = np.atleast_1d(data[3])
  rotstr = np.atleast_1d(data[4])
  NTran  = energy.size
  NCols  = np.shape(data)[0] # data is transposed

  # Check if the correct broadening is assigned to the lineshape
  if OPT['LShp'] == "lorentzian":
    if OPT['BTyp'] != "HWHH":
      c.error(" You can't specify a sigma with Lorentzian LineShape: you have to specify HWHH!","Broadening")

  if OPT['BRed'] == True:    # Read brodening from input file or assign the parsed value
    if NCols < 6 : c.error(" You asked to read the broadening but you forgot to put the data into the table!","Broadening")
    broad = np.atleast_1d(data[5])
  else:
    broad = np.array([OPT['BVal']]*NTran)

  # If user request Gaussian with HWHH broadeining we report it to sigma
  if OPT['LShp'] == "gaussian" and OPT['BTyp'] == "HWHH":
    broad = broad/np.sqrt(2*np.log(2))

  if OPT['UBrd'] == "eV"  : broad *= c.PhyCon['eV2wn']
  
  ETrMax = np.max(energy)
  ETrMin = np.min(energy)
  # Automatic spectrum range
  if autorange is True:  autowindow(ETrMin,ETrMax,broad)

  # Print spectrum limits
  printlimits()

  # Number of points
  OPT['NPts'] = np.round((OPT['XMax']-OPT['XMin'])/OPT['XStp'])
  
  print 
  print " E min  : %8.0f cm-1" % OPT['XMin']
  print " E max  : %8.0f cm-1" % OPT['XMax']
  print " Step   : %8.0f cm-1" % OPT['XStp']
  print " NPts   : %8d       " % OPT['NPts']
  print(" Prefix Output files %s " % OPT['FOut'])
  print

  print " Number of electronic transitions : %3d" % NTran


  if ETrMax > OPT['XMax'] or ETrMin < OPT['XMin'] :
    print " WARNING! The selected spectrum range is too small to plot all transtions!"
    print " WARNING! --min < %8.0f" % ETrMin
    print " WARNING! --max > %8.0f" % ETrMax

  # Line shape broadening
  print " Lineshape function      : %s"  % OPT['LShp']
  print " Broadening is given in  : %s"  % OPT['UBrd']


  # Shift energies
  if OPT['shift'] != None:
    print " Excitation energies will be shifted by %8.4f eV (%10.2f cm^-1)" % (OPT['shift'],OPT['shift']*c.PhyCon['eV2wn'])
    print energy
    print OPT['shift']*c.PhyCon['eV2wn']
    energy += OPT['shift']*c.PhyCon['eV2wn']
    print energy

  # Call the function to write the spectrum into a file
  x,OD,LD,CD = specalc(NTran,energy,dipo,dipold,rotstr,broad)

  findpeaks("OD",x,OD.sum(axis=1),broad)
  findpeaks("CD",x,CD.sum(axis=1),broad)
  findpeaks("LD",x,LD.sum(axis=1),broad)

  savespec("OD",x,OD,OPT)
  savespec("LD",x,LD,OPT)
  savespec("CD",x,CD,OPT)

 
  # Visualize the spectrum
  if OPT['IPlt'] == True :
    print("\n ... Calling specview.py to visulize and save the spectra")
    try: 
      import specview
    except:
      c.error('The specview module could not be called!','Spectrum')
    specview.view(OPT,x,energy,OD,LD,CD,dipo,dipold,rotstr)
 
  print("\n Done!\n")
