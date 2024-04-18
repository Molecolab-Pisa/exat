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
# anacoeff.py 
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
import sys, os, glob
import numpy as np
import argparse as arg
import readdata 
import common as c 


#
# Compute inverse participation ratio
#
def calc_IPR(Coeff):
  IPR = 1/ np.sum(Coeff**4,axis=1)
  return IPR


def buildformat(NSt,C2,Thresh,NTranChrom):
  CgtT = np.where( C2 > Thresh)
  fmt = ('%2d) %10.4f %5.1f %8.1f %5.2f')
  fmt += ' | '
  k = 0
  for j in NTranChrom:
    for i in range(j):
      if C2[k] > Thresh : 
        fmt += '% 0.2f'
        #print C2[i],fmt
      else:
        fmt += '    %0.0f'
        #print C2[i],fmt
      k += 1
    fmt += ' | '
  return fmt


#
# Print resume table
# ------------------------------------------------------------------

def print_summary(ExcEne,Dip2,RotStr,Coeff2,Thresh,NTran):

  NSt = len(ExcEne)
  IPR = calc_IPR(Coeff)

  #fmthdr = '                          '
  NTranChrom = []
  for i in NTran:
    NTranChrom.append(len(i))
    #fmthdr += '|'+'%6s'+'|'
  header = ' #      E (ev)   mu2   RotStr   IPR '

  header += '|  G(1)      |      G(2)  |       C(3)      |        A(4)     |       C(5)      |       U(6)      |     G(9)   |    G(10)   |     U(11)       |     G(12)  |       C(13)     |      C(14)      |'
  #  header = ['G(1)','G(2)','C(3)','A(4)','C(5)','U(6)','G(9)','G(10)','U(11)','G(12)','C(13)','C(14)']
  #print ( fmthdr % tuple(header))

  print('\n'+230*('-'))
  print(header)
  print(230*('-'))

  for i in range(NSt):
    fmt    = buildformat(NSt,Coeff2[i],Thresh,NTranChrom)
    OutRow = (i+1,ExcEne[i],Dip2[i],RotStr[i],IPR[i])+(tuple(Coeff2[i]))
    print(( fmt % OutRow))

  print("Minimum IPR = %5.1f on ExcState %2d" % (np.min(IPR),np.argmin(IPR)+1))
  print("Maximum IPR = %5.1f on ExcState %2d" % (np.max(IPR),np.argmax(IPR)+1))
  print("Average IPR = %5.1f (sigma = %5.1f)" % (np.average(IPR),np.std(IPR)))

  pass


#
# Read results file
# ------------------------------------------------------------------

def read_results(InResFile):

  # Check if the external file exists
  c.checkfile(InResFile)

  # Read file and store data
  Data = np.loadtxt(InResFile,comments='#')

  ExcEne = Data[:,1]
  Dip2   = Data[:,2]
  RotStr = Data[:,4]

  return ExcEne,Dip2,RotStr

#
# Read coefficients file
# -------------------------------------------------------------------

def read_coeff(InCoeffFile):

  # Check if the external file exists
  c.checkfile(InCoeffFile)

  # Initialize lists
  ExcEne = [] ; Coeff  = [] ; Coeff2 = []

  # Read external file and store excitation energies and coefficients
  with open(InCoeffFile,'r') as f:
    while True:
      Tmp = list(map(float,f.readline().split()))
      try:
        ExcEne.append(Tmp[0])  # Excitonic energy in eV
        Coeff.append(Tmp[2:])
      except:
        break

  # Transform lists into arrays
  ExcEne = np.array(ExcEne)
  Coeff  = np.array(Coeff)

  # Compute the Coeff2
  Coeff2 = Coeff**2

  return ExcEne,Coeff,Coeff2



if __name__ == "__main__":

  parser = arg.ArgumentParser(description="Analyze the coefficients matrix")
  parser.add_argument('-i',help='File Containing the coefficient matrix as in the output of EXAT program',default='diag.dat')
  parser.add_argument('-c','--chromlist',help='Chromlist File as in EXAT calculations',default='chromlist.in')
  parser.add_argument('-r','--results',help='File containing the results',default='results.out')
  parser.add_argument('-t','--threshold',help='Threshold for printing coeffs',default=0.1,type=float)
  args = parser.parse_args()

  InCoeffFile = args.i
  InChromFile = args.chromlist
  InResFile = args.results
  Thresh = args.threshold

  c.OPT['seltran'] = True
  c.ExtFiles['crlist'] = InChromFile
  
  
  ChromList,NChrom,NTran = readdata.ReadChromList()
  ExcEne,Coeff,Coeff2 = read_coeff(InCoeffFile) 

  NSt = len(ExcEne)
 
  ExcEne,Dip2,RotStr  = read_results(InResFile) 

  print_summary(ExcEne,Dip2,RotStr,Coeff2,Thresh,NTran)
