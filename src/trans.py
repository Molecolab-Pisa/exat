#!/usr/bin/env python
#
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
# trans.py TRANSITION PROPERTIES MODULE
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

from __future__ import print_function
import os,sys
import readdata, math, util
import numpy    as np
import argparse as arg
import common   as c

#****************************************************

# Compute Rotational Strength

def RotStrength(EEN,Cent,coeff,DipoLen,EXCDipoLen,DipoVel,MagInt,RxDel,Site):


  if c.OPT['RCalc'] == "mu":
    if c.v():
      print(" ... using approximate treatment (electric dipoles only)") 
    EXCRot = Rapprox(EEN*c.PhyCon['Town'],coeff,DipoLen,Cent)

  elif c.OPT['RCalc'] in ("mag","velocity"):
    if c.v():
      print(" ... Using full treatment (electric and magnetic dipoles)\n") 

    # Compute some excitonic properties
    EXCDipoVel = EXCalc(coeff,DipoVel) 
    EXCMagInt  = EXCalc(coeff,MagInt)
    EXCMagExt  = EXCalc(coeff,RxDel)

    EXCMagInt  = EXCalc(coeff,MagInt)
    EXCMagExt  = EXCalc(coeff,RxDel)
    EXCMagTot  = EXCMagInt+EXCMagExt 
  
    if c.v(0):
      print("     Splitting mu-mu and mu-mag contributions") 
 
    if c.OPT['RCalc'] == "mag":  
      # Compute extrinsic R Strength (mu-mu) 
      RExt,AngExt = RVelIso(EEN,EXCDipoVel,-EXCMagExt)
      # Compute intrinsic R Strength (mu-mag)
      RInt,AngInt = RVelIso(EEN,EXCDipoVel,-EXCMagInt)
      # Compute    approx R Strength (Length)
      EXCRotMu    = Rapprox(EEN*c.PhyCon['Town'],coeff,DipoLen,Cent)

    # Compute     total R Strength 
    EXCRot,Angles = RVelIso(EEN,EXCDipoVel,-EXCMagTot)
  

    if c.OPT['RCalc'] == "mag":
      Out = np.column_stack((EXCRotMu,RInt,RExt,EXCRot))
      np.savetxt(c.OutFiles['rstrength'],Out,fmt="%10.4f")
    
      if c.v(0):
        print(" "*5+"-"*58) 
        print("  "+(4*"%15s") % ("Mu (Approx)","R Int","R Ext","RTot"))
        print(" "*5+"-"*58) 
        for i in range(len(EXCRot)):
          print("  "+" %14.2f"*4 % (EXCRotMu[i],RInt[i],RExt[i],EXCRot[i]))
        print()
  
  elif c.OPT['RCalc'] == "length":
    if c.v():
      print(" ... using full treatment in Length formulation")
    # Extrinsic part
    EXCRot = Rapprox(EEN*c.PhyCon['Town'],coeff,DipoLen,Cent)
    # Add Intrinsic part
    EXCRot += RLenIso(EXCDipoLen,-EXCMagInt)[0]

  elif c.OPT['RCalc'] == "newlen":
    print(" ... PROVA!")
    EXCMagInt  = EXCalc(coeff,MagInt)
    EXCMagExt  = EXCalc(coeff,RxDel)
    EXCMagTot  = EXCMagInt+EXCMagExt 

    SiteAU     = Site/c.PhyCon['ToeV'] 
    DipoVel1   = (DipoVel.T/SiteAU).T
    EXCDipoV1  = EXCalc(coeff,DipoVel1)

    EXCRot     = RLenIso(EXCDipoV1,EXCMagTot)[0]

  if c.v(1): 
    print(' '*5+'Sum of rotatory strengths: %10.2f\n' % (EXCRot.sum()))
  return EXCRot  

# Compute Linear Dichroism

def LinDichro(Dipo):

  if c.OPT['LDAxis']   == "x" : Axis = np.array([1.0,0.0,0.0])
  elif c.OPT['LDAxis'] == "y" : Axis = np.array([0.0,1.0,0.0]) 
  elif c.OPT['LDAxis'] == "z" : Axis = np.array([0.0,0.0,1.0]) 
  else:
    Axis = np.array(c.OPT['LDAxis'],dtype=float)
    c.OPT['LDAxis'] = "[%4.1f,%4.1f,%4.1f]" % tuple(Axis)
  if c.v():
    print("   ... LD Spectrum computed with respect to the %s axis"\
          % c.OPT['LDAxis'])

  LD = []
  
  Axis /= np.linalg.norm(Axis)
  for i in range(len(Dipo)):
    DotProd = np.dot(Dipo[i,:],Axis)
    Dip2    = np.dot(Dipo[i,:],Dipo[i,:])
    LD.append(1.5*(3*DotProd**2-Dip2))
  return LD

#****************************************************

def EXCalc(coeff,V):
  EXV = np.dot(coeff,V)
  return EXV

#****************************************************
# Compute the rotational strenght in the velocity form
# Gaussian util : ecdten.F
#
# EEN,Del,RxDel can be floats or np.arrays
def RVelIso(EEN,Del,RxDel):

  FactRV = 471.4436078822227
  Fact   = -FactRV/EEN
  RVel   = np.sum(Del*RxDel,axis=-1)
  Angle  = np.arccos(RVel/(np.linalg.norm(Del,axis=-1)*np.linalg.norm(RxDel,axis=-1)))
  RVel   = RVel*Fact/2
  Angle  = np.degrees(Angle)

  return RVel,Angle


#****************************************************
# Compute the rotational strenght in the velocity form
# Gaussian util : ecdten.F
#
# r,RxDel can be floats or np.arrays
def RLenIso(r,RxDel):

  FactRL = 471.4436078822227
  Fact   = FactRL
  RLen   = np.sum(r*RxDel,axis=-1)
  Angle  = np.arccos(RLen/(np.linalg.norm(r,axis=-1)*np.linalg.norm(RxDel,axis=-1)))
  RLen   = RLen*Fact/2

  return RLen,Angle

#****************************************************

def calc_tprod(dipo,cent):
  #Returns upper triangular matrix of products (mu_i X mu_j)*(r_j-r_i)
  n = dipo.shape[0]
  A = np.zeros((n,n))

  for i in range(n):
    for j in range(i+1,n):  
     A[i,j] = np.dot(cent[i] - cent[j],np.cross(dipo[i],dipo[j]))
  return A 

def Rapprox(energy,coeff,dipo,cent):

  ntot = len(energy)

  #compute triple products
  tprod = calc_tprod(dipo,cent)
  if c.OPT['savetprod']:
    if c.v(0):
      print(' ... save triple product matrix in %s'\
        % (c.OutFiles['tprod']))
    np.savetxt(c.OutFiles['tprod'],tprod,fmt='%10.5f')

  #R_k = coeff_k * tprod * coeff_k.T
  R = (np.dot(coeff,tprod)*coeff).sum(axis=1)

#
# Unit of R at this ponints are dipole^2 (a.u.) (units of centers (Ang) are 
# cancelled with the unit of energy (cm^-1) by inserting the factor 1E8 in the
# prefactor.
#
  R *= energy
  R *= 1E-8*np.pi
#
# Returns the Excitonic Rotational Strenghts in 10^-40 esu^2 cm^2 (same unit as
# GaussView CD Spectrum)
#
  DipAu2cgs = ((2.541746E-18)**2)*1E+40
  R *= DipAu2cgs

  return R

# -----------------------------------------------------------------------------
