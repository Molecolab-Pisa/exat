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
# loganalizer.py UTIL
# *************************************
#

#
# loganalizer.py is a program to analize your gaussian output
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



# Standard Python Modules
from __future__ import print_function
import sys,os,glob
import numpy    as np
import argparse as arg
import common   as c
import util     as u
import readdata as r
from excsystem import ExcSystem


# EEN,Del,RxDel can be floats or np.arrays
def RVelIso(EEN,Del,RxDel):

  FactRV = 471.4436078822227
  Fact   = -FactRV/EEN
  RVel   = np.sum(Del*RxDel,axis=-1)
  Angle  = np.arccos(RVel/(np.linalg.norm(Del,axis=-1)*np.linalg.norm(RxDel,axis=-1)))
  RVel   = RVel*Fact/2
  Angle  = np.degrees(Angle)

  return RVel,Angle



def checkdipovel(dipolen,dipovel,site):

  dipovel1 = -dipovel/(site[:,None]/c.PhyCon['ToeV']) 
  
  delta = np.linalg.norm(dipovel1 - dipolen, axis=-1)

  norm1 = np.linalg.norm(dipovel1,axis=-1)
  norm2 = np.linalg.norm(dipolen ,axis=-1)
  cos   = np.sum(dipolen*dipovel1,axis=-1)/(norm1*norm2)
  angle = np.degrees(np.arccos(cos))

  deltaperc = delta/norm2*100
  
  print()
  for i in range(s.NTran[0]):
    print(" Nabla %2d ----- delta = %6.2f%%  angle = %8.4f " % (i+1,deltaperc[i],angle[i]))
  
  pass


# *****************************************************************************

#
# If called by command line, read the Gaussian log file provided
#

if __name__ == "__main__" : 

  # Read input line
  parser = arg.ArgumentParser(description="Read Gaussian .log file and extract some useful data...")
  parser.add_argument('--version','-V',action='version',version=c.PROGVERS,
                      help="Show program's full version number, and exit")  

  parser.add_argument('-v',help='Increase the verbosity of the output',action="count",default=0)
  parser.add_argument('logfile',help='Gaussian output log file to process')
  parser.add_argument('-o','--out',help='Prefix for output files')
  parser.add_argument('--seltran',help='Activate seltran keyword',action="store_true")
  parser.add_argument('--cent',help='How to compute the chormophore center',nargs="*",default="geom")
  parser.add_argument('--anadipo',help='Analysis of the transition dipole orientation with respect to a certain axis')
  parser.add_argument('--mu',help='Prefix for output files',action="store_true")
  parser.add_argument('--ang',help='Compute angle',action="store_true")
  parser.add_argument('--spec',help='Extract data for produce uv and cd spectrum',action="store_true")
  parser.add_argument('--reorient',help='Indicate 2 atoms as reference direction',
          nargs=2,type=int,default=None)
  parser.add_argument('--checkvel',action="store_true",
          help='Check velocity dipole moments')
  args = parser.parse_args()

  c.welcome()

  # Set options
  InLogFile   = args.logfile
  if args.out:
    c.OPT['OutPrefix'] = args.out
    c.setoutfiles()

  c.OPT['verbosity'] = args.v
  c.OPT['seltran']   = args.seltran
  c.OPT['reorient']  = args.reorient
  c.OPT['Cent']      = args.cent
  if args.anadipo:
    c.OPT['anadipo'] = True
    c.ExtFiles['refaxis'] = args.anadipo

  fixang = True#args.fixang
  velcheck = args.checkvel

  # Check if the .log files
  c.checkfile(InLogFile)

  # Read Gaussian Log files
  #out = r.readgaulog(InLogFile)
  anum,xyz,site,dipo,dipovel,mag,NAtom,NTran,NChrom = r.readgaulog(InLogFile)


  dipo     = np.array(dipo)
  dipovel  = np.array(dipovel)
  mag      = np.array(mag)
  site     = np.array(site)


  NAtom    = [NAtom]


  # Compute the chromophore center
  Cent = u.calchromcent(NAtom,anum,xyz,[NTran])
  
  s = ExcSystem([(1,list(range(1,NTran+1)))],site=site,Cent=np.asarray(Cent),
          DipoLen=dipo,DipoVel=dipovel,Mag=mag)
  s.add_geom(anum,xyz,NAtom)
  # Tweak for couplings
  s.coup = np.zeros(0)
  s.buildmatrix()


  # Select transtion (if requested by user)
  #  if c.OPT['seltran'] == True:
  #    Coup = 0.0
  #    Site,Dipo,DipoVel,Mag,Cent,NewCoup = readdata.seltran(site,dipo,dipovel,mag,Cent,Coup)

  # Analyze electric transition dipole moment
  c.ChromList = [InLogFile.split('.')[0]]
  if c.OPT['reorient'] is not None:
    u.reorientdipo(dipo,dipovel,mag,coup)


# Compute internal magnetic moment (DipoMag is gauge-dependent)
  Cent = np.array(Cent)
  RxDel = np.cross(Cent/c.PhyCon['ToAng'],dipovel)
  MagInt = mag - RxDel
  MagTot = mag
  MagExt = MagTot-MagInt


  print('\n')
  print(" Electric transition dipoles: ")
  for i in range(s.NTran[0]):
    NormDip =  np.linalg.norm(dipo[i])
    print(" Trans %2d  %8.4f  mu  = %8.4f %8.4f %8.4f  ->  %8.4f " % (i+1,site[i],dipo[i][0],dipo[i][1],dipo[i][2],NormDip))

  #c.OPT['verbosity'] = -1
  if c.OPT['anadipo'] is not None:
    u.dipoanalysis(s)

  if velcheck:
    checkdipovel(dipo,dipovel,site)

  print('\n')
  print(" Magnetic transition dipoles (intrinsic): ")
  for i in range(s.NTran[0]):
    NormMagInt =  np.linalg.norm(MagInt[i])
    print(" Trans %2d  %8.4f  mag = %8.4f %8.4f %8.4f  ->  %8.4f " % (i+1,site[i],MagInt[i][0],MagInt[i][1],MagInt[i][2],NormMagInt))

    # Dipole analysis with internal magnetic moment
  if c.OPT['anadipo'] is not None:
    s1 = s.copy()
    s1.DipoLen = MagInt
    u.dipoanalysis(s1,True)

# Save visudipo
  u.savegeom(anum,xyz)
  u.savevisudipo(s,s.DipoLen,s.Mag)

# Recompute rotatory strenght
  R,angle = RVelIso(site/c.PhyCon['ToeV'],dipovel,-MagTot)
  R = np.array(R)

# Save results.out
  SqDip = np.sum(dipo**2,axis=1)
  u.resout(site*c.PhyCon['eV2wn'],SqDip,site*0.0,R)

  print("\n Done! \n")
  

