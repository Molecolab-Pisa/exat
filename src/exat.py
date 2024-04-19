#!/usr/bin/env python
#
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
# exat.py MAIN PROGRAM
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

# Import python modules
from __future__ import print_function
import sys, os, glob 
import numpy as np

# Import COMMON modules
import common as c   # Common
import util   as u   # Utilities

# Import WORKING modules
import readdata      # Read data (from Gaussian output)
import trans         # Compute the excitonic dip and R
from excsystem import ExcSystem, ChromTranList

# *****************************************************************************
#
# (0) Read the options from input line and set the OPT dictionary
#    

c.welcome()
u.GetOpts()

# *****************************************************************************
#
# (1) Read the input parameters from external files provided
#     by the user or extract data from different Gaussian versions
#

print("\n (1) Read data") 

# System-independent call (including seltran)
system,Sel = readdata.Read()


# *****************************************************************************
#
# 2. Build the Excitonic Matrix
#
print("\n (2) Build the excitonic Hamiltonian") 

if c.OPT['ScaleCoup'] != 1.0:
  if c.v():
    print("   ... applying coupling scaling factor"+ \
  " of %4.1f" % c.OPT['ScaleCoup'])
  system.coup *= c.OPT['ScaleCoup']

if c.OPT['CleanCoup'] > 0.0:
  thr = c.OPT['CleanCoup']
  if c.v():
    print("   ... applying treshold of %4.1f" % thr)
  system.coup[abs(system.coup) < thr] = 0.0

# If requested, scale dipole (Length)
if  c.OPT['ScaleDipo'] != 1.0 : 
  print(" ... All transition dipoles (length) will be scaled by factor %8.4f"\
     % c.OPT['ScaleDipo'])
  system.DipoLen  *= c.OPT['ScaleDipo']  # Scale transition dipole moments 

system.buildmatrix()

# Apply selection of transitions
if c.OPT['seltran']:
  system = system.seltran(Sel)

# All other modifications are AFTER seltran!

# If requested reorient TrDipo and change the sign to corresponding couplings
# This option works only with 1 transition per chromophore
if ( c.OPT['reorient'] is not None ) : 
  if c.OPT['read'] == 'external': 
    c.error(' Reorient Dipoles is not possible with read external')
  else:
    print(" ... The direction of transition dipole moments will be changed according to selected axis")
    if c.OPT['forcedipo'] == True : 
      print(" ... The orientation of transition dipole moments will forced parallel to selected axis")
    # Reorient DipoLen and DipoVel, on the basis of DipoLen (!!)
    system = u.reorientdipo(system)

if c.OPT['read'] != 'external':
  # If requested, modify transition centers TODO
  if c.OPT['ModCent']: Cent = u.modcent(Cent)
  # If requested, modify electric transition dipoles
  if c.OPT['ModDipoLen']: DipoLen = u.moddipo(DipoLen,'dipo')
  if c.OPT['ModDipoMag']: Mag     = u.moddipo(Mag,'magdipo')
  # If requested, modify site energies
  if c.OPT['ModSite']: 
    system = u.modsite(system)
  # Possibly modify couplings here
  if c.OPT['ModCoup']: 
    system = u.modcoup(system)
      


  # If requested performs the transition dipole angle analysis
  # This should be done AFTER coupling modification
  # meaning that the input dipoles and couplings are consistent
  if c.OPT['anadipo']: 
    system, _ = u.dipoanalysis(system)

  # If requested, scale dipoles and couplings
  if c.OPT['ScaleTran']: 
    system = u.scaletran(system) 

      
# Save the exciton Hamiltonian to file
system.savematrix(c.OutFiles['matrix'])
# If verbosity is > 1 , print couplings and distances and save coup.out file
if c.v(1): u.prtcoup(system)
# If verbosity is > 1 , print site energies and save site.out file
if c.OPT['verbosity'] > 1 : u.prtsite(system)

# *****************************************************************************
#
# 3. Diagonalize the Excitonic Hamiltonian
#

print("\n (3) Diagonalize the excitonic Hamiltonian") 
system.diagonalize()
u.savediag(system.energy,system.coeff,system.coef2)

# Convert some quantities to A.U.
EEN  =  system.energy/c.PhyCon['Town']  # Excitonic Energies (Hartree)

# *****************************************************************************
#
# 4. Perform excitonic calculations
#

RxDel  = np.cross(system.Cent/c.PhyCon['ToAng'],system.DipoVel)

# Compute internal magnetic moment (DipoMag is gauge-dependent)
MagInt = system.Mag - RxDel

print("\n (4) Compute the excitonic properties") 

EXCDipoLen = trans.EXCalc(system.coeff,system.DipoLen)

# Compute Linear Absorption Spectrum
print("\n ... Compute the Linear Absorption Intensities") 
EXCDipo2   = np.sum(EXCDipoLen**2,axis=1)

# Compute Linear Dichroism Spectrum
print("\n ... Compute the Linear Dichroism Intensities") 
LD = trans.LinDichro(EXCDipoLen)

# Compute Rotational Strength ...
print("\n ... Compute the Circular Dichroism Intensities") 
EXCRot = trans.RotStrength(EEN,system.Cent,system.coeff,system.DipoLen,
                             EXCDipoLen,
                             system.DipoVel,MagInt,RxDel,system.Site)

#
# 5. Print out the reuslts
#

# Visualization options

# Save file for DIPOLE visualization in VMD
# This is here to print also excitonic dipo
# Print intrinsic magnetic moments
if system.has_geom: 
  u.savegeom(system.anum,system.xyz.tolist())
  u.savevisudipo(system,EXCDipoLen,-MagInt)
 #np.savetxt(c.OutFiles['cent'],Cent,fmt='%12.5f',
 #  header='Coordinates of centers (ang)')

# Print out results.out
u.resout(system.energy,EXCDipo2,LD,EXCRot)

# Save the exciton system to file 
system.save(c.OutFiles['exatdata'])

print("\n Done! \n")

