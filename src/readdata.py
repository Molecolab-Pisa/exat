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
# readdata.py INPUT READING MODULE
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

# Standard Python Modules
import sys,os,glob,fileinput,re,copy
import subprocess as sp
import numpy  as np
import argparse as arg

# Import COMMON Modules
from excsystem import ChromTranList,ExcSystem,load_npz
import common as c
import util   as u

# *****************************************************************************
# 
# MAIN function to read data
#

def Read():
  
  # system-dependent read
  SelChromList = None
  if c.OPT['seltran']:
    ChromList = ReadChromList()
    # Save selection
    SelChromList = ChromTranList(ChromList.iteritems())
  elif c.OPT['read'] not in ('g16','load'):
    ChromList = ReadChromList()

  if c.OPT['read'] == 'load':
    return load_npz(c.ExtFiles['load']),SelChromList

  if c.OPT['read'] == 'g16':
    if c.v():
      print  " ... excpected gaussian version: g16"
      print  " ... reading in %s .log file" % c.OPT['logfile']
    Site,Dipo,DipoVel,Mag,Coup,NTran,anum,xyz,NAtom = readgaulog36(c.OPT['logfile'])
    ChromList = ChromTranList( ((i+1),[]) for i in range(len(NTran)) )
    ChromList.set_NTran(NTran)
  
  elif c.OPT['read'] == 'external':
  
    if c.v(): 
      print " ... all parameters are read from external files" 
      print " ... reading %s file to set the chromophore list" % c.ExtFiles['crlist'] 
    Cent,Site,Dipo,DipoVel,Mag,Coup = ReadExternal(ChromList)
  

  if c.OPT['read'] != 'external':
    if c.v():
      print " ... compute the center of each transtion" 
    Cent = u.calchromcent(NAtom,anum,xyz,ChromList.NTran)
  
  # Compute dipole couplings
  if c.OPT['coup'] == 'PDA' :
    if c.v():
      print " ... using dipole-dipole couplings based on electric transition dipole moments" 
      print "     REFRACTION INDEX = %5.2f" % c.OPT['refrind'] 
    Coup,Kappa = u.coupforster(Cent,Dipo,ChromList.NChrom,ChromList.NTran)
    Kappa    = np.array(Kappa)
  else:
    Kappa = None
  
  
  Cent     =  np.array(Cent)
  DipoLen  =  np.array(Dipo)                    # Electric tranistion dipole moments (a.u.)
  DipoVel  =  np.array(DipoVel)                 # Nabla (a.u.)
  Mag      =  np.array(Mag)
  Site     = np.array(Site)*c.PhyCon['eV2wn']
  Coup     = np.array(Coup)

  system   = ExcSystem(ChromList,Site,Coup,Cent,DipoLen,
                        DipoVel,Mag,Kappa)

  if c.OPT['read'] != 'external':
    system.add_geom(anum,xyz,NAtom)

  return system,SelChromList


# *****************************************************************************
# 
# Read data from external files
#
def ReadExternal(ChromList):

  # Read External files:
  if c.v():
    print "     > CHROMLIST         : %s" % c.ExtFiles['crlist']
    print "     > site energies     : %s" % c.ExtFiles['insite']
    print "     > coupling          : %s" % c.ExtFiles['incoup']
    print "     > elec. tr. dipoles : %s" % c.ExtFiles['dipo']
    print "     > chrom centers     : %s" % c.ExtFiles['incent']

  NChrom,NTran  =  ChromList.NChrom, ChromList.NTran

  if c.v():  
    print " ... number of chromophores         : %3d" % NChrom 
    print " ... N. transitions per chromophore : %s"  % NTran 

  # Read Site Energies
  Site = np.zeros(sum(NTran))
  Site = np.loadtxt(c.ExtFiles['insite'])[:,1:].flatten()

  # Check if the file exists and read the files
  c.checkfile(c.ExtFiles['incoup'])
  c.checkfile(c.ExtFiles['incent'])
  c.checkfile(c.ExtFiles['dipo'])
  Coup = np.loadtxt(c.ExtFiles['incoup'],dtype="float",ndmin=1)
  Cent = np.loadtxt(c.ExtFiles['incent'],dtype="float")
  Dipo = np.loadtxt(c.ExtFiles['dipo'],dtype="float")

  DipoVel = np.copy(Dipo)
  Mag     = np.copy(Dipo)
  Cent    = np.array(Cent)
  Site    = np.array(Site)
  Dipo    = np.array(Dipo)/c.PhyCon['ToDeb']
  DipoVel = np.array(DipoVel)
  Mag     = np.array(Mag)
  Coup    = np.array(Coup)


  return Cent,Site,Dipo,DipoVel,Mag,Coup



# *****************************************************************************
# 
# Read ChromList file
#
def ReadChromList():
  " Reads the ChromList file"        

  InFile = c.ExtFiles['crlist']
  c.checkfile(InFile)

  ChromList = ChromTranList.from_crfile(InFile)

  if c.OPT['seltran']:
      nonsel = [n == 0 for n in ChromList.NTran]
      if any(nonsel):
        Tmp = ChromList.Chrom[nonsel.index(True)]
        ErrMsg = "Chromophore %s has no transition selected" % Tmp
        c.error(ErrMsg,"ReadChromList")

  return ChromList



# *****************************************************************************
# 
# Prepare the file lists for gdvh23 version
#

def ChromListBuilder(chromlist):
  chromfilelist = []
  if c.OPT['env'] == 'vac'     : suffix = "_vac"
  elif c.OPT['env'] == 'mmpol' : suffix = ""
  for chrom in chromlist:
    chromfilelist.append(chrom+"/"+chrom+suffix+".log")
  return chromfilelist

def CoupListBuilder(ChromList):
  nchrom = ChromList.NChrom  

  crlist = ChromList.Chrom

  CoupList = []
  CoupFileList = [] ; Prefix = ""
  if c.OPT['env'] == 'vac'   : suffix = "_vac"
  if c.OPT['env'] == 'mmpol' : suffix = ""
  for i in range(nchrom):
    for j in range(i+1,nchrom):
      coup = crlist[i]+"."+crlist[j]
      coupfileij = Prefix+"V_"+coup+"/V_"+coup+suffix+".log"
      if c.OPT['coup'] != 'forster' :
        #c.checkfile(coupfileij)
        CoupFileList.append(coupfileij)
        CoupList.append([i,j])
  return CoupList,CoupFileList

def preplists(ChromList):
  ChromFileList = ChromListBuilder(ChromList.Chrom)
  CoupList,CoupFileList = CoupListBuilder(ChromList)
  return ChromFileList,CoupList,CoupFileList

# *****************************************************************************
#
# Read into a property gdvh23 .log file to extract site energies, dipoles ...
# should be called for each prop file.
#

def readgaulog(logfile):

  FINDNAt = "NAtoms="

  # Check if the file exists
  c.checkfile(logfile)

  # If the file exists then open it, read all data and close it
  infile = open(logfile,'r')
  data = infile.read().split("\n")
  infile.close()

  # Inizialize lists
  xyz  = [] ; anum = [] ; site = []; dipo = [] ; dipovel = [] ; mag = [] ; NAtoms = []

  # Read the loaded file line by line
  i = 80
  while i < len(data):
    
    # Looks for NTran:
    if "9/41=" in data[i]:
      NTran = int(data[i].split("9/41=")[1].split(",")[0])

    # Looks for NChrom:
    if "62=" in data[i]:
      NChrom = int(data[i].split("62=")[1].split(",")[0])
      if c.OPT['read'] == 'gdvh36': NAtoms = [0]*NChrom 


    # Looks for the atomic coordinates 
    String = "Input orientation:"
    if String in data[i]:
      atom = True
      j = 0
      while atom == True :
        tempcen = data[i+5+j].split()
        if len(tempcen) != 6:
          atom = False
        else:
          xyz.append(map(float,tempcen[3:6]))
          anum.append(int(tempcen[1]))
          j = j + 1

    # Extracts the dipole moments (length)
    elif "electric dipole" in data[i]:
      atom = True
      j = 0
      while atom == True :
        tempdip = data[i+2+j].split()
        if len(tempdip) != 6 :
          atom = False
        else:
          dipo.append(map(float,tempdip[1:4]))
          j = j + 1

    # Extracts the dipole moments (velocity)
    elif "transition velocity dipole" in data[i]:
      atom = True
      j = 0
      while atom == True :
        tempdip = data[i+2+j].split()
        if len(tempdip) != 6 :
          atom = False
        else:
          dipovel.append(map(float,tempdip[1:4]))
          j = j + 1

    # Extracts the magnetic moment
    elif "transition magnetic dipole" in data[i]:
      atom = True
      j = 0
      while atom == True :
        tempdip = data[i+2+j].split()
        if len(tempdip) != 4 :
          atom = False
        else:
          mag.append(map(float,tempdip[1:4]))
          j = j + 1

    # Extracts the site energy values
    elif "nm " in data[i]:
      site.append(float(data[i].split()[4]))

    # Compute the center of the chromophore
    elif FINDNAt in data[i]:
      NAtoms.append(data[i].split(FINDNAt)[1].split()[0])

    i += 1


  NChrom = 1
  NAtoms = int(NAtoms[0])
  NTran  = len(site)

  return (anum,xyz,site,dipo,dipovel,mag,NAtoms,NTran,NChrom)


# *****************************************************************************
#
# Read into gdvh36 .log files to extract site energies, dipoles ...
# Also extract couplings
#

def readgaulog36(logfile):

  # Check if the file exists
  c.checkfile(logfile)

  # Inizialize lists
  xyz  = [] ; anum = [] ; site = []; dipo = [] ; dipovel = [] ; mag = [] ; NAtoms = []

  # Initialize NTran to Gaussian default
  # to avoid errors later
  NTran = 3

  # What coupling to look for
  if   c.OPT['coup'] == 'trden' or c.OPT['coup'] == 'total':  
    if c.OPT['read'] == 'gdvh36': 
      FindCoup = "> TOTAL" 
    else:
      FindCoup = "Total coupling"
    DoCoup   = True 
  elif c.OPT['coup'] == 'coulomb':  
    if c.OPT['read'] == 'gdvh36': 
      FindCoup = "> Coulomb"
    else: 
      FindCoup = "Coulomb"
    DoCoup   = True 
  elif c.OPT['coup'] == 'PDA':
    print("   ... skipping couplings reading: coupling will be computed using point dipole approximation\n")
    DoCoup   = False 
  else:
    c.error("Confused in coupling choice!","readgaulog36")


  # If the file exists then open it and read line by line
  with open(logfile,'r') as f:
    # We split the file reading in 2 loops to improve performance
    # Read properties
    while True:  
      line = f.readline()
      if not line: break
      if 'Electronic Coupling for Excitation Energy Tranfer' in line: break

      # Looks for NTran:
      if "9/41=" in line:
        try:
          NTran = int(line.split("9/41=")[1].split(",")[0])
        except:
          NTran = int(line.split("9/41=")[1].split("/")[0])
          

      # Looks for NChrom:
      if "62=" in line:
        NChrom = int(line.split("62=")[1].split(",")[0])
        NAtoms = [0]*NChrom 
        FragAt = [None]*NChrom


      # Extract number of atoms in gdvh36
      if line[1:19] == 'Symbolic Z-matrix:':
        kk = 0
        while True:
          line = f.readline()
          if line[0:9] == ' Charge =':
            pass 
          elif '(fragment=' in line.lower():
            d = line.split()[0]
            IFrag = int(d.split('=')[1].split(')')[0].split(',')[0])
            NAtoms[IFrag-1] += 1
            # Assign atom to fragment
            try:    FragAt[IFrag-1].append(kk)
            except: FragAt[IFrag-1] = [kk]
            kk += 1
          elif "NAtoms=" in line:
            NAtmTot = int(line[11:16])
            break

      # Looks for the atomic coordinates 
      elif "Input orientation:"  in line:
        atom = True
        j = 0
        while atom == True :
          line = f.readline()
          # Skip 4 lines
          if j < 4: j += 1;  continue
          tempcen = line.split()
          if len(tempcen) != 6:
            atom = False
          else:
            xyz.append(map(float,tempcen[3:6]))
            anum.append(int(tempcen[1]))
            j = j + 1

      # Extracts the dipole moments (length)
      elif "electric dipole" in line:
        atom = True
        j = 0
        while atom == True :
          pos  = f.tell()
          line = f.readline()
          # Skip 1 line
          if j < 1: j += 1;  continue
          tempdip = line.split()
          if len(tempdip) != 6 : 
            f.seek(pos) #go back one line
            break
          else:
            dipo.append(map(float,tempdip[1:4]))
            j = j + 1

      # Extracts the dipole moments (velocity)
      elif "transition velocity dipole" in line:
        atom = True
        j = 0
        while atom == True :
          pos  = f.tell()
          line = f.readline()
          if j < 1: j += 1;  continue
          tempdip = line.split()
          if len(tempdip) != 6 : 
            f.seek(pos) #go back one line
            break
          else:
            dipovel.append(map(float,tempdip[1:4]))
            j = j + 1

      # Extracts the magnetic moment
      elif "transition magnetic dipole" in line:
        atom = True
        j = 0
        while atom == True :
          pos  = f.tell()
          line = f.readline()
          if j < 1: j += 1;  continue
          tempdip = line.split()
          if len(tempdip) != 4 : 
            f.seek(pos) #go back one line
            break
          else:
            mag.append(map(float,tempdip[1:4]))
            j = j + 1

      # Extracts the site energy values
      elif "nm " in line:
        site.append(float(line.split()[4]))
    #end while

    # read couplings
    if DoCoup: 
      if c.v(): print  " ... reading couplings in %s file using %s values" % (c.OPT['logfile'],c.OPT['coup'])
      Cdtyp = [('Ch1',int),('Tr1',int),('Ch2',int),('Tr2',int),('Coup',float)]
      Coup = []

      if c.OPT['read'] == 'gdvh36':
        # read couplings for gdvh36-molecolab 
        if c.OPT['coup'] == 'total': ExplCoup = []
        while True:  
          line = f.readline()
          if not line: break
          elif FindCoup in line: 
            Coup += [(int(line[25:29]),int(line[31:35]),int(line[43:47]),int(line[49:53]),float(line[55:69]))]
          elif c.OPT['coup'] == 'total' and '> Explicit MMPol'  in line:
            ExplCoup += [(int(line[25:29]),int(line[31:35]),int(line[43:47]),int(line[49:53]),float(line[55:69]))]
      else:
        # read couplings for gdvh36 (plain)
        while True:
          line = f.readline()
          if not line: break
          elif 'Frag=' in line and 'State=' in line:
            Ch1,St1,Ch2,St2 = (int(line[6:9]),int(line[16:19]),int(line[45:48]),int(line[55:58]))
            # Read other lines
            while True:
              line = f.readline()
              if not line: break
              if FindCoup in line: 
                CoupValue = float(line[31:43])*c.PhyCon['eV2wn']
                break
            # swap chrom 1 and 2 for compatibility with order of gdvh36-molecolab
            Coup += [(Ch2,St2,Ch1,St1,CoupValue)]

      # Transform in np.array and sort 
      Coup =  np.array(Coup,dtype=Cdtyp)
      Coup.sort(order=['Ch1','Ch2','Tr1','Tr2','Coup'],axis=0)
      Coup = Coup['Coup']# Get rid of indices
      if c.OPT['coup'] == 'total' and len(ExplCoup) == len(Coup): 
      # Compute and subtract explicit term
        ExplCoup =  np.array(ExplCoup,dtype=Cdtyp)
        ExplCoup.sort(order=['Ch1','Ch2','Tr1','Tr2','Coup'],axis=0)
        Coup -= ExplCoup['Coup']# Get rid of indices
         
    # End read couplings
    else:
      print("   ... skipping couplings reading: coupling will be computed using point dipole approximation\n")
      Coup = None
  # File is closed
  
  #Transform NTran in list for compatibility
  NTran = [NTran]*NChrom

  # Reorder anum and xyz frag-wise order
  anum1 = np.array(anum)
  xyz1  = np.array(xyz)
  anum = []; xyz = []
  for i in range(NChrom):
    anum = anum + anum1[FragAt[i]].tolist()
    xyz  = xyz  +  xyz1[FragAt[i]].tolist()

  return site,dipo,dipovel,mag,Coup,NTran,anum,xyz,NAtoms



def readcouph36(logfile): 

  c.NChrom = NChrom
  c.NTran = NTran # Save NTran to Common 

  return Site,Dipo,DipoVel,Mag,Cent,NewCoup,NewKappa

# *****************************************************************************
 
def printlocal(Site,Dipo,Dip2,Cent):

  for n in range(c.OPT['NChrom']):
    print("\nChromophore     : %10s" % n )
#    print("Center of trans : %10.4f %10.4f %10.4f" % (Cent[n][0],Cent[n][1],Cent[n][2]))
    print'----------------------------------------------------------'
    print' #   E (eV)      mx      my     mz        Dip2   (a.u.) '
    print'----------------------------------------------------------'
#
    for i in range(c.OPT['NTran'][n]):
      k = sum(c.OPT['NTran'][0:n])+i
      print("%2d %8.4f   %6.3f  %6.3f  %6.3f    %6.3f "%(i+1,Site[k],Dipo[k][0],Dipo[k][1],Dipo[k][2],Dip2[k]))

# *****************************************************************************

