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
# viewlevel.py UTIL
# *************************************
#
# ------------------------------------------------------------------------
#
#   This module allows the visualization of energy level diagram for
#   site energies and excitonic energies
#
# ------------------------------------------------------------------------
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
import os, sys
import numpy as np
import common as c
from excsystem import ExcSystem, load_npz
import argparse as arg
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt

CMAP = mpl.cm.afmhot

#
# Plot energy levels
#
def plot_levels(thresh,site,energy,coeff,labels=None):

  N = len(energy)

  if labels is None:
    labels = ['%2d' % (I+1) for I in range(len(site))]

  # Hide x-axis and set reasonable range
  ax = plt.gca()
  ax.set_ylabel('Energy (cm$^{-1}$)')
  ax.xaxis.tick_top()
  ax.tick_params(axis='x',which='both',bottom='off',top='off',pad=15)

  xtics  = [1.0]
  xlabel = ['Site']
  for i in range(N):
    xtics  += [xtics[i]+1.0]
    xlabel += ['Exciton %d' % (i+1)]


  ax.xaxis.set_ticks(xtics)
  ax.xaxis.set_ticklabels(xlabel)
  bfont = mpl.font_manager.FontProperties(size=14, weight='bold')

  for lab in ax.xaxis.get_ticklabels():
    lab.set_fontproperties(bfont)
 
  ax.set_xlim([xtics[0]-0.5,xtics[-1]+0.5])
  # Margin
  plt.subplots_adjust(left=0.2,right=0.9)

  # Plot site energy levels
  for I,E in enumerate(site):
    ax.plot([0.8,1.2],[E,E],'b-',linewidth=2.0,picker=5)
    # Labels site
    ax.text(0.6,E,labels[I],fontsize=16,horizontalalignment='right',verticalalignment='bottom')

  # Plot excitonic levels
  thresh = thresh/1e2 
  if Dipo is not None:
    DipoK = Dipo.copy()/c.PhyCon['ToDeb']
    D2sum = np.sum(DipoK**2)
    D2 = []
    for K in range(N):
      # Transform dipoles to new/intermediate basis
      DipoK  = np.dot(coeff[K],DipoK)
      D2K    = np.sum(DipoK**2,axis=1)
      D2.append(D2K)
    D2    = np.asarray(D2)
    D2max = np.max(D2)

  for K in range(N):
    C2     = coeff[K]**2
    for I,E in enumerate(energy[K]):
      if Dipo is None:
        ax.plot([xtics[K+1]-0.2,xtics[K+1]+0.2],[E,E],'b-',linewidth=2.0,picker=5)
      else:
        #lw = D2[I]/D2sum*20.0
        #ax.plot([xtics[K+1]-0.2,xtics[K+1]+0.2],[E,E],'b-',linewidth=lw,picker=5)
        #ax.plot([xtics[K+1]-0.2,xtics[K+1]+0.2],[E,E],'-',color=CMAP(intens),linewidth=2.0,picker=5)
        intens = D2[K,I]/D2max
        pos = xtics[K+1]-0.2 + intens*0.4
        ax.plot([xtics[K+1]-0.2,xtics[K+1]+0.2],[E,E],'k-',
                        linewidth=1.0,alpha=0.3,picker=5)
        ax.plot([xtics[K+1]-0.2,pos           ],[E,E],'b-',
                        linewidth=2.0,picker=5)

      contribs = np.where(C2[I] > thresh)[0]

      if K == 0:
        for J in contribs:
          ax.plot([xtics[0]+0.2,xtics[1]-0.2],[site[J],E],'k--')
      else:
        for J in contribs:
          ax.plot([xtics[K]+0.2,xtics[K+1]-0.2],[energy[K-1][J],E],'k--')

  plt.draw()
  plt.show()

  return


#
# Diagonalization routines
#

def eighsort(matrix):
  e,v = np.linalg.eigh(matrix)
  index = e.argsort()
  e = e[index]
  v = v[:,index]
  return e,v

def diagonalize(M):

  # Diagonaliztion
  eigenval, eigenvec = eighsort(M)
  eigenvec = eigenvec.T
  dim = len(eigenval)

  # Convert the results into array-type objects
  eigenvec = np.asarray(eigenvec).reshape((dim,dim))
  eigenval = np.asarray(eigenval).reshape(-1)

  return(eigenval,eigenvec)



# *****************************************************************************
#
# M A I N  P R O G R A M
#

if __name__ == "__main__" :

  # Read input line
  parser = arg.ArgumentParser(description="Tool for energy level visualization")
  parser.add_argument('-v',help='Increase the verbosity of the output',action="count")
  parser.add_argument('-t','--threshold',help='Threshold for connecting levels (%%)',type=float,default=10.0)
  parser.add_argument('-s','--strongcoupling',help='Define threshold for strong couplings (cm^-1)',type=float,default=None)
  parser.add_argument('--partition',help='Select different sets of states to partition the Hamiltonian',
      nargs='+',default=None)
  #parser.add_argument('--dipo',help='Use dipoles to estimate dipole strength of excitons',default=None)
  #parser.add_argument('matrix',help='Matrix file containing site energies and couplings',default='matrix.dat')
  parser.add_argument('exatfile',default='exat.npz',
      help='Exat npz file containing the exciton system')
  args = parser.parse_args()

  # Set the options
  #InMatFile  = args.matrix
  exatfile   = args.exatfile
  threshold  = args.threshold
  strongcoup = args.strongcoupling
  partition  = args.partition

  c.welcome()
  print("\n > viewlevel.py module")
  print("   A tool for visualizing the energy levels \n")

  if args.v > 0  : OPT['verbosity'] = args.v

  # Read the excitonic matrix
  c.checkfile(exatfile)
  print(("\n Reading the excitonic system from %s" % (exatfile) ))

  excsys = load_npz(exatfile)

  excsys.diagonalize()



  M    = excsys.H
  Dipo = excsys.DipoLen
  dim  = np.shape(M)
  print((" ... matrix dimension are: %s" % str(dim)))

  labels = []
  for Chrom in excsys.ChromList:
    for j in excsys.ChromList[Chrom]:
      labels.append("%s-%d" % (Chrom,j))
    

  if strongcoup is not None:

    print("\n Advanced visualization with two step diagonalization")
    print((" strong coupling threshold = %6.1f cm-1" % strongcoup))

    # Build H0 Hamiltonian
    print(" ... build H0 Hamiltonian (contains only stron couplings)")
    mat = np.ndarray.flatten(M)
    for i,x in enumerate(mat):
      if abs(x) < strongcoup :
        mat[i] = 0.0
    M0 = np.reshape(mat,dim)

    # Diagonalize H0 Hamiltonian
    print(" ... diagonalize H0 Hamiltonian")
    e0,v0 = diagonalize(M0)

    # Transform H in the basis of v0
    print(" ... project original Hamiltonian in the basis of H0 to obtain H1")
    M1 = np.dot(np.dot(v0,M),v0.T)

    # Diagonalize H1
    print(" ... diagonalize the H1 Hamiltonian")
    e1,v1 = diagonalize(M1)

    # Plot energy level
    print(" ... print energy levels")
    plot_levels(threshold,M.diagonal(),[e0,e1],[v0,v1],labels)

  elif partition is not None:
    sets = [np.array(c.stringconverter(x))-1 for x in partition]
    
    print("\n Advanced visualization with two step diagonalization")
    print(" Hamiltonian partitioning")
    print(sets)
    print()
    #
    print(" ... build H0 Hamiltonian (only intra-set)")
    M0 = np.diag(M.diagonal()) # Set diagonal part
    # For every set, add couplings within the set
    for SS in sets:
      nn = len(SS)
      for I in range(nn):
	for J in range(nn): M0[SS[I],SS[J]] = M[SS[I],SS[J]]
    
    # Diagonalize H0 Hamiltonian
    print(" ... diagonalize H0 Hamiltonian")
    e0,v0 = diagonalize(M0)

    # Transform H in the basis of v0
    print(" ... project original Hamiltonian in the basis of H0 to obtain H1")
    M1 = np.dot(np.dot(v0,M),v0.T)

    # Diagonalize H1
    print(" ... diagonalize the H1 Hamiltonian")
    e1,v1 = diagonalize(M1)

    # Plot energy level
    print(" ... print energy levels")
    plot_levels(threshold,M.diagonal(),[e0,e1],[v0,v1],labels)
     

  else:                 

    print("\n Diagonalizing the exciton matrix" )

    # Diagonalize M
    e,v = diagonalize(M)

    # Plot energy level
    print(" ... print energy levels")
    plot_levels(threshold,M.diagonal(),[e],[v],labels)


  print("\n Done! \n")

