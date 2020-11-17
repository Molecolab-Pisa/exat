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
# util.py UTILITIES MODULE
# *************************************
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
#
import os,sys,glob
import numpy  as np
from scipy.spatial import distance
import common as c
import argparse as arg


# *****************************************************************************
#
# Parse input line
#

class CustomHelpFormatter(arg.HelpFormatter):
  # 
  def _split_lines(self, text, width):
  #return all but include blanklines 
    return ['']+arg.HelpFormatter._split_lines(self, text, width)+['']

# Get Options with argparse
def GetOpts():

  # Read input line
  parser = arg.ArgumentParser(description="""EXAT - EXcitonic Analysis Tool: 
  Performs an excitonic calculation from the 
  output of Gaussian gdvh23dev or gdvh36eet log file(s). """,
  formatter_class=lambda prog: CustomHelpFormatter(prog,
   max_help_position=36,indent_increment=3))

  # Specify which chromlist file should be used
  parser.add_argument('chromlist',metavar='ChromListFile',
    help='''File containing the chromophore list (Default: chromlist.in).
    The file should contain one chromophore per line.''',
    default="chromlist.in",nargs='?')

  # Version
  parser.add_argument('--version','-V',action='version',version=c.PROGVERS,
     help="Show program's full version number, and exit")  

  # Input options
  grpinp = parser.add_argument_group(' Input options')
  read_group = grpinp.add_mutually_exclusive_group()
  read_group.add_argument('--readexternal','-e',action="store_true",
     help='Read all quantities from external files')
  read_group.add_argument('--load','-l',default=False,
     help='''Load an exat .npz file instead of 
            reading from Gaussian/External files''')
  grpinp.add_argument('--log','-i',metavar='LogFile',
     help='Gaussian logfile name (for gdvh36 only)')

  # Output options
  grpout = parser.add_argument_group(' Output options')
  grpout.add_argument('--out','-o',metavar='Prefix',
    help='Prefix for all the output files.')
  grpout.add_argument('--savecoeff',action='store_true',
    help='Save excitonic coefficient in a numpy file for later use.')
  grpout.add_argument('--savetprod','--tprod',action='store_true',
    help='Save matrix of triple products R_ij*(mu_i x mu_j).')

  # Site Energies
  grpsit = parser.add_argument_group(' Options for site energies')
  grpsit.add_argument('--insite','--modsite','--site',metavar='SiteEnergyFile',
    help='Modify site energies according to SiteEnergyFile',default=None)

  # Centers options
  grpcnt = parser.add_argument_group(' Center options')
  read_grpcnt = grpcnt.add_mutually_exclusive_group()
  read_grpcnt.add_argument('--cent',metavar='Center',nargs="*",
    help='How to compute the chromophore center. \
    Default is the geometric center of the chromophore.',default="geom")
  read_grpcnt.add_argument('--incent',metavar='CenterFile',
    help='Modify chromophore centers according to CenterFile',default=None)


  # Couplings
  grpcou = parser.add_argument_group(' Options for couplings') 
  grpcou.add_argument('--coup', 
   choices=['trden','PDA','coulomb'],    
   help='''Choose the coupling calculation type. 
   trden includes all terms. 
   Coulomb includes only the coulomb term. 
   PDA computes coupling with point-dipole approxmation.
   Default is trden.''')
  grpcou.add_argument('--cleancoup',metavar='Threshold',type=float,
    help='Set to zero all the couplings below Threshold')
  grpcou.add_argument('--seltran',action="store_true",
    help="""Activate selection of transitions on the basis of ChromListFile. 
    ChromListFile should contain, for each line, the transitions to consider. """)
  grpcou.add_argument('--refrind',  metavar='n',type=float,
    help='''Specify the refraction index n for PDA coup calculations. 
    Default is 1.0 (vacuum).''')
  grpcou.add_argument('--scalecoup',metavar='Factor',type=float,
    help='Scale all couplings by Factor')
  grpcou.add_argument('--modcoup',metavar='CoupFile',
    help='Modify electronic couplings according to the CoupFile',default=None)
  grpcou.add_argument('--incoup',metavar='CoupFile',
    help='Read all ectronic couplings from CoupFile',default=None)

  # Transition Dipoles
  grpdip = parser.add_argument_group(' Options for transition dipoles ')
  grpdip.add_argument('--scaledipo',type=float,
    help='Scale all transition dipole moments')
  grpdip.add_argument('--reorient',default=None,nargs=2,
    type=int,metavar=('ID1','ID2'),
    help='''Id of two atoms that define the axis 
    with respect to which the TrDip will be reoriented''')
  grpdip.add_argument('--forcedipo',action="store_true",
    help='Force TrDip to be aligned along the axis specified by --reorient keyword')
  grpdip.add_argument('--anadipo',metavar='FILE',
    help='''Analysis of the transition dipole orientation with respect to the
    reference frame indicated in FILE.''')
  grpdip.add_argument('--scaletran',metavar='ScaleFile',
    help='Request scaling of all dipoles and couplings as indicated in ScaleFile.')
  grpdip.add_argument('--indipo','--dipolen',metavar='DipoFile',default=None,
    help='Modify electric dipole moments according to DipoFile')
  grpdip.add_argument('--inmag','--dipomag',metavar='DipoFile',default=None,
    help='Modify magnetic dipole moments according to DipoFile')

  # CD related
  rsttitle = ' CD related options: \n  Options to compute the rotational strength'
  grprst = parser.add_argument_group(rsttitle)
  rcalc_group = grprst.add_mutually_exclusive_group()
  rcalc_group.add_argument('--mu',action='store_true',
    help='Use approximate treatment (electric dipoles only). This is the default.')
  rcalc_group.add_argument('--mag',action='store_true',
    help='''Use full treatment (electric and magnetic dipoles), 
    splitting the electric and magnetic contributions''')
  
  # Advanced options
  grpadv = parser.add_argument_group(' Advanced and expert options ')
  grpadv.add_argument('--ldaxis',
    help='Choose the reference axis for LD calculation')

  # Verobisity of the output and printing options
  grpprt = parser.add_argument_group(' Printout and verbosity options')
  grpprt.add_argument('-v',action='count',
    help='Increase the verbosity of the output')
  grpprt.add_argument('-q',action="store_true",
    help='Quiet output: minimum printout')
  grpprt.add_argument('--prtcoup-threshold',metavar='Threshold',
    type=float,default=100.0,
    help='Threshold for coupling printout with -vv option enabled (default: 100cm^-1)')

  # Parse args
  args = parser.parse_args()

  # Set OPT dictionary
  if args.v > 0   : c.OPT['verbosity'] = args.v
  elif args.q     : c.OPT['verbosity'] = -1

  c.ExtFiles['crlist'] = args.chromlist

  # Check I/O options
  if args.readexternal:
    c.OPT['read'] = 'external'
    # With --readexternal, --seltran option is activated by default
    c.OPT['seltran'] = True  
    if args.insite is not None: c.ExtFiles['insite'] = args.insite
    if args.incoup is not None: c.ExtFiles['incoup'] = args.incoup
    if args.indipo is not None: c.ExtFiles['dipo']   = args.indipo
    if args.incent is not None: c.ExtFiles['incent'] = args.incent

  elif args.load:
    c.OPT['read'] = 'load'
    c.ExtFiles['load'] = args.load

  else:
    if not args.log:
      c.error("You have to indicate the Gaussian .log file!","Input parse")
    else:
      c.OPT['logfile']  = args.log
      c.OPT['read'] = 'g16'

  if args.mu == True       : c.OPT['RCalc'] = "mu"
  elif args.mag == True    : c.OPT['RCalc'] = "mag"
  else                     : c.OPT['RCalc'] = "mu"  # Default 

  if args.cent : c.OPT['Cent'] =  args.cent

  if args.out      : 
    c.OPT['OutPrefix'] = args.out
    c.setoutfiles()
  if args.savecoeff:
    c.OPT['savecoeff'] = True
  if args.savetprod:
    c.OPT['savetprod'] = True

  if args.seltran  : c.OPT['seltran']  = args.seltran
  if args.scalecoup is not None: c.OPT['ScaleCoup']  = args.scalecoup
  if args.scaledipo : c.OPT['ScaleDipo']  = args.scaledipo
  if args.reorient  : c.OPT['reorient']  = args.reorient

  if args.forcedipo  : 
    if c.OPT['reorient'] != None:
      c.OPT['forcedipo']  = args.forcedipo
    else:
      c.error("To force a dipole you have also to specify --reorient keyword",
      "argparse")

  if args.anadipo  : 
    c.OPT['anadipo'] = True
    c.ExtFiles['refaxis'] = args.anadipo
  if args.scaletran : 
    c.OPT['ScaleTran']  = True 
    c.ExtFiles['ScaleTran']  = args.scaletran

  if args.insite is not None:
    c.OPT['ModSite']   = True
    c.ExtFiles['insite'] = args.insite

  if args.indipo is not None:
    c.OPT['ModDipoLen']  = True
    c.ExtFiles['dipo']   = args.indipo

  if args.inmag is not None:
    c.OPT['ModDipoMag']  = True
    c.ExtFiles['magdipo'] = args.inmag

  if args.modcoup is not None:
    c.OPT['ModCoup']  = True
    c.ExtFiles['modcoup'] = args.modcoup

  if args.incent is not None:
    c.OPT['ModCent']  = True
    c.ExtFiles['incent'] = args.incent

  if args.coup      : c.OPT['coup']   = args.coup
  if args.refrind   : c.OPT['refrind']   = args.refrind
  if args.cleancoup : c.OPT['CleanCoup']   = args.cleancoup
  if args.ldaxis    : c.OPT['LDAxis'] = args.ldaxis
  setldaxis()

  c.OPT['CouPrtThr'] = args.prtcoup_threshold

  if c.v(1): c.PrintOpt()



# *****************************************************************************
#
# Compute the chromophore centers
#
def calchromcent(NAtom,anum,xyz,NTran):

  Cent = []
  Init = 0

  if "geom" in c.OPT['Cent'] : ICalc = "geom"
  elif "mass" in c.OPT['Cent'] : ICalc = "mass"
  else: 
    ICalc = "select"
    ISel  = c.stringconverter(c.OPT['Cent'][0])
    ISel  = (np.array(ISel)-1).tolist()
    if c.v():
      print " The following atoms will be selected to compute chromophore center: "
      print ISel

  for i in range(len(NAtom)):
    End  = Init + NAtom[i]
    IAnu = anum[Init:End]
    IXYZ = xyz[Init:End]
    Init = End

    if    ICalc == "geom" : 
      ICent = np.average(IXYZ,axis=0)
    elif  ICalc == "mass" : 
      IMas = anu2awg(IAnu)
      ICent = np.average(IXYZ,axis=0,weights=IMas)
    elif  ICalc == "select" :
      SAnu = list(IAnu[i] for i in ISel)
      SXYZ = list(IXYZ[i] for i in ISel)
      SMas = anu2awg(SAnu)
      ICent = np.average(SXYZ,axis=0,weights=SMas)

    # Replicate the center for every transiton...
    for j in range(NTran[i]):
      #Cent.extend(calcent(IAnu,IXYZ))
      Cent.extend(ICent)

  Cent = np.reshape(Cent,(sum(NTran),3)).tolist()
  return Cent

# *****************************************************************************
#
# Reorient transition dipoles and change the coupling sign
#
#def reorientdipo(DipoLen,DipoVel,Mag,Coup):
def reorientdipo(system):

  exc = system.copy()
    
  Init = 0
  At1,At2 = c.OPT['reorient'] # direction: from 1 to 2
  K=0
  for I in range(exc.NChrom):
    # Extract xyz of Chrom I
    End  = Init + exc.NAtom[I]
    IXYZ = np.array(exc.xyz[Init:End])
    Init = End
    # Define axis 
    axis  = IXYZ[At2-1]-IXYZ[At1-1]
    # Normalize
    axis /= np.linalg.norm(axis)
    for J in range(exc.NTran[I]):
      sign = np.sign(np.dot(axis,exc.DipoLen[K]))

      if sign < 0 : 
        if c.v():
          print "Chrom: %3d  Trans: %3d  Transition electric dipole moment will be reoriented"\
             % (I+1,J+1)
        exc.DipoLen[K] *= -1 
        # Change signs to DipoVel and Mag as well
        exc.DipoVel[K] *= -1 
        exc.Mag[K]     *= -1 
        #exc.coup = changesign(I,J,exc.NChrom,exc.NTran,exc.Coup)
        exc.H[K]       *= -1
        exc.H[:,K]     *= -1

      if c.OPT['forcedipo'] == True:
        if exc.NTran[I] > 1 : c.error("NTran > 1: you cannot force all dipoles to be parallel to a specific axis!","reorientdipo")
        exc.DipoLen[K] = axis*np.linalg.norm(exc.DipoLen[I])

      K += 1

    exc.update_sitecoup()
    exc.H[exc.H == 0.0] = 0.0

  return exc

# *****************************************************************************
#
# Transition dipole analysis
#
def buildrefframe(A,B,C):
  # Creates a reference frame from three points.
  # x is the A->B axis, and y is in the ABC plane
  # z = x cross y is defined by the couterclockwise ABC rotation

  ux = (B-A)/np.linalg.norm(B-A)  # Versor defined A ----> B
  AC = (C-A)
  ACprojx = np.dot(AC,ux) # projection on ux
  # check quasi-linear dependence
  if abs(ACprojx/np.linalg.norm(AC)) > np.cos(np.radians(10.0)) and c.v():
    print " WARNING: choosen atoms are almost aligned!"
    print " Angle between vectors: %8.2f degrees; threshold: %8.2f degrees" \
    % (np.degrees(np.arccos(ACprojx/np.linalg.norm(AC))),10)  
  # subtract the projection on x and normalize
  AC -= ACprojx*ux
  uy  = AC/np.linalg.norm(AC)
  # create z axis (x cross y)
  uz = np.cross(ux,uy)
  return ux,uy,uz

def changerefframe(oldframe,newframe,V):
  # Gives components of V in the new reference frame
  # works only for orthonormal ref frames
  transform = np.dot(oldframe,newframe.T)
  U = np.dot(transform,V)
  return U

def diporient(ChromFrame,crlist,DipoLen,Theta,Phi,RefTheta,RefPhi,ForceDipo):
  # Orient dipoles for all chromophores, return new dipoles
  # ChromFrame: list of molecular ref. frames, one per chromophore
  
  if c.v(): 
    print " ... Dipoles of the following transitions will be oriented: "

  # Mask reference with ForceDipo
  mask = np.logical_not(ForceDipo)
  RefTheta[mask[:,0]] = Theta[mask[:,0]]
  RefPhi[mask[:,1]]   = Phi[mask[:,1]]

  # Reference axis from theta/phi in spherical coords
  rft = np.radians(90. - RefTheta)
  rfp = np.radians(RefPhi)
  RefAxes = [np.sin(rft)*np.cos(rfp),np.sin(rft)*np.sin(rfp),np.cos(rft)]
  RefAxes = np.asarray(RefAxes).T # ref axes in molecule frame

  for i in range(crlist.NChrom): 
    k0 = sum(crlist.NTran[:i])
    Chrom = crlist.Chrom[i]
    for j in range(crlist.NTran[i]):
      k = k0+j
      if not any(ForceDipo[k]): continue
      # Compute new ref axis in lab frame
      ref   = np.dot(RefAxes[k],ChromFrame[i]).T
      # Dipole magnitude
      MDipo    = np.linalg.norm(DipoLen[k])
      DipoLen[k] = ref*MDipo # Force orientation
  
      if c.v():
        print " Chrom %7s  Tran %3d   theta = %10.4f  phi = %10.4f   Dipole = %8.4f %8.4f %8.4f" \
        % (Chrom,j+1,RefTheta[k],RefPhi[k],DipoLen[k,0],DipoLen[k,1],DipoLen[k,2])

  return DipoLen
 
def dipoanalysis(system,nochange=False):
  
  exc = system.copy()  
  NChrom = exc.NChrom
    
  if c.v():
    print "\n ==== Transition Dipole Analysis Function ==== \n"
    print " Computing angle between transition dipole moments and a given reference frame"
    print " The reference frame is defined from the positions of three (groups of) atoms "
    print " Note that the sign of the angle is related to the third (group of) atoms" 
    print 

  c.checkfile(c.ExtFiles['refaxis'])  
  if c.v():
    print " > Reference file: %s \n" % c.ExtFiles['refaxis']

  # Read reference file
  with open(c.ExtFiles['refaxis'],'r') as f:
    lines = [line.split('#')[0] for line in f.readlines()]
  InData  = [l for l in lines if l.strip()]

  # Process reference file
  NLines = len(InData)
  if NLines != NChrom: 
    c.error("Reference file must contain the same number of chromphores of the input","dipoanalysis")
  RefChrom = []
  AA = [] ; BB = [] ; CC = []

  NTranTot = sum(exc.NTran)
  RefTheta = np.zeros(NTranTot)#[]
  RefPhi   = np.zeros(NTranTot)#[]
  ForceDipo  = np.zeros((NTranTot,2),dtype=bool) #[]

  Init = 0
  for i in range(NChrom):
    # Line for chromophore i
    Data = InData[i].split()

    if len(Data) < 4 or ( (len(Data) % 2) != 0) :
      c.error("Check reference input! Number of coulmns for at least one data is wrong!","dipoanalysis")

    RefChrom.append(Data[0])
    RefGrpA = c.stringconverter(Data[1])
    RefGrpA = (np.array(RefGrpA)-1).tolist()
    RefGrpB = c.stringconverter(Data[2])
    RefGrpB = (np.array(RefGrpB)-1).tolist()
    RefGrpC = c.stringconverter(Data[3])
    RefGrpC = (np.array(RefGrpC)-1).tolist()
    End  = Init + exc.NAtom[i]
    IXYZ = exc.xyz[Init:End]
    AA.append(np.average(list(IXYZ[j] for j in RefGrpA),axis=0))
    BB.append(np.average(list(IXYZ[j] for j in RefGrpB),axis=0))
    CC.append(np.average(list(IXYZ[j] for j in RefGrpC),axis=0))

    # Read whether to change the dipole orientation
    for j in range(exc.NTran[i]):
      # theta[Chrom i, Tr j] in Data[col]
      # phi    "        "    in Data[col+1] 
      col = 4 + 2*j 
       
      if len(Data) <= col: break

      jj = sum(exc.NTran[:i])+j

      # theta
      if Data[col].lower().endswith('f'): 
        RefTheta[jj] = float(Data[col][:-1])
        ForceDipo[jj] = [True,False]
      else:  
        RefTheta[jj] = float(Data[col])
        ForceDipo[jj] = [False,False]

      # phi   
      if Data[col+1].lower().endswith('f'):
        RefPhi[jj] = float(Data[col+1][:-1])
        ForceDipo[jj][1] = True
      else:
        RefPhi[jj] = float(Data[col+1])
        ForceDipo[jj][1] = False
      
    Init = End

  # Reference axis from theta/phi in spherical coords
  rft = np.radians(90. - RefTheta)
  rfp = np.radians(RefPhi)
  RefAxes = [np.sin(rft)*np.cos(rfp),np.sin(rft)*np.sin(rfp),np.cos(rft)]
  RefAxes = np.asarray(RefAxes).T # ref axes in molecule frame

  Theta = []; Phi = []

  ChromFrame = np.zeros((NChrom,3,3))

  k = 0
  for i in range(NChrom):
    # Build new reference frame
    ux,uy,uz = buildrefframe(AA[i],BB[i],CC[i])
    ChromFrame[i] = np.asarray([ux,uy,uz])
    for j in range(exc.NTran[i]):
      # Evaluate the angle
      Dipo  = exc.DipoLen[k]
  
      ref   = np.dot(RefAxes[k],ChromFrame[i]).T # ref axis in lab frame
      sign  = np.dot(ref,Dipo)
      
      if sign < 0 and not nochange: 
        if c.v():
          print " Chrom: %3d  Trans: %3d  Transition sign will be changed "\
             % (i+1,j+1)
        exc.DipoLen[k] *= -1 
        # Change signs to DipoVel and Mag as well
        exc.DipoVel[k] *= -1 
        exc.Mag[k]     *= -1 
        exc.H[k]       *= -1
        exc.H[:,k]     *= -1

      UDipo = Dipo/np.linalg.norm(Dipo) 
      theta = np.arccos(np.dot(UDipo,uz)) 
      theta = 90.0 - np.degrees(theta) # angle wrt the plane
      phi   = np.arctan2(np.dot(UDipo,uy),np.dot(UDipo,ux))
      phi   = np.degrees(phi) # angle wrt x axis (AB)

      # save values
      Theta.append(theta)
      Phi.append(phi)
      k += 1

  Theta = np.array(Theta)
  Phi = np.array(Phi)

  # Print out the results of the dipo analysis
  if c.v():
    print 
    k = 0
    for i in range(NChrom):
      for j in range(exc.NTran[i]): 
        print " Chrom %7s  Tran %3d   theta = %10.4f  phi = %10.4f"\
          % (exc.ChromList.Chrom[i],j+1,Theta[k],Phi[k])
        if c.v(0): 
          print "   Dipole: " + "  %10.4f  %10.4f  %10.4f " % tuple(exc.DipoLen[k])
        k += 1   

  
  if c.v(): print 


  if ForceDipo.any(): 
    # Note: diporient masks RefTheta/Phi with Theta/Phi  
    #       where ForceDipo is False
    exc.DipoLen = diporient(ChromFrame,exc.ChromList,exc.DipoLen,
            Theta,Phi,RefTheta,RefPhi,ForceDipo)
  elif c.v():
    print " ... Dipoles will not be changed"

  # End
  if c.v():
    print "\n ====  End of Transition Dipole Analysis  ==== \n"

  exc.H[exc.H == 0.0] = 0.0  
  return exc,(np.array(Theta),np.array(Phi))


# *****************************************************************************
#
# Transition dipoles and coupling scaling (ScaleTran)
#
def scaletran(system):
# Scales electric and magnetic dipoles of transition as defined in c.ExtFiles['ScaleTran']
# Scales all couplings accordingly

  exc    = system.copy()
  crlist = exc.ChromList

  reffile = c.ExtFiles['ScaleTran']
  c.checkfile(reffile)  
  if c.v():
    print 
    print " ... scaling transition densities using scaling factors from file: %s" % reffile 

  # Read reference file
  with open(reffile,'r') as f:
    lines = f.readlines()

  for line in lines: 
    data  = line.split()

    if len(data) < 3:
      c.error("Wrong format for scaletran")

    chrom = data[0]

    try:
      itran = int(data[1])
    except:
      c.error("Transition number not understood: %s" % data[1])

    #k = c.TrIdx(chrom,itran) - 1
    k = crlist.TrIdx(chrom,itran) - 1
  
    try:
      fact = float(data[2])
    except:  
      c.error('Scale factor not understood: %s' % data[2])
    if fact != 1.0:
      if c.v():
        print "     Chrom: %10s   Tran: %3d  scale transition density by %10.6f" \
               % (chrom,itran,fact)
        exc.DipoLen[k] *= fact
        exc.DipoVel[k] *= fact
        exc.Mag[k]     *= fact
        exc.H[:,k]     *= fact
        exc.H[k,:]     *= fact

        exc.H[k,k] = system.H[k,k]

  return exc

# *****************************************************************************
#
# Modify site energies 
#
def modsite(system):
  exc = system.copy()
  crlist = system.ChromList

  reffile = c.ExtFiles['insite']
  c.checkfile(reffile)  

  # Read reference file
  with open(reffile,'r') as InFile:
    lines = InFile.readlines()


  for line in lines: 
    data  = line.split()

    if len(data) < 3:
      c.error("Wrong format for modsite")

    chrom = data[0]

    try:
      itran = int(data[1])
    except:
      c.error("Transition number not understood: %s" % data[1])

    k = crlist.TrIdx(chrom,itran) - 1
  
    if "+" in data[2] or "-" in data[2]: shift = True  # check if we want to shift or not
    else: shift = False

    try:
      Ene = float(data[2])
    except:  
      c.error('Site energy value not understood: %s' % data[2])
    if Ene is not None:
      if shift == False:
        if c.v():
          print "     Chrom: %10s   Tran: %3d  site energy set to %8.5f eV"\
                % (chrom,itran,Ene)
        exc.H[k,k] = Ene*c.PhyCon['eV2wn']
      elif shift == True:
        if c.v():
          print "     Chrom: %10s   Tran: %3d  site energy sfhit of %8.5f eV"\
                % (chrom,itran,Ene)
        exc.H[k,k] += Ene*c.PhyCon['eV2wn']
 
  exc.update_sitecoup()

  return exc

# *****************************************************************************
#
# Modify electric transition dipole moments
#
def moddipo(Dipo,name):

  # Expected format:
  # Chrom Tran mu_x mu_y mu_z

  reffile = c.ExtFiles[name]
  c.checkfile(reffile)
  if c.v():
    print " Dipoles (Debye) will be read from %s" % (reffile)

  fmt = " Chrom: %5s Tran: %3d  Orig dipo = %8.4f %8.4f %8.4f    New dipo = %8.4f %8.4f %8.4f"

  with open(reffile,'r') as InFile:
    lines = InFile.readlines()
  for line in lines:
    data  = line.split()

    if len(data) < 5:
      c.error("Wrong format for moddipo")

    chrom = data[0]
    try:
      itran = int(data[1])
    except:
      c.error("Transition number not understood: %s" % data[1])

    k = c.TrIdx(chrom,itran) - 1
    newdip = np.array(data[2:5],dtype=float)

    # Possibly go back to atomic units
    if name == 'dipo': newdip /= c.PhyCon['ToDeb']

    if c.v():
      print fmt % (chrom,k+1,Dipo[k][0],Dipo[k][1],Dipo[k][2],newdip[0],newdip[1],newdip[2])
    Dipo[k] = newdip
        
  return Dipo


# *****************************************************************************
#
# Modify centers 
#
def modcent(Cent):
  # Expected format:
  # Chrom [Tran] X Y Z 

  reffile = c.ExtFiles['incent']
  c.checkfile(reffile)
  if c.v():
    print " Centers will be read from %s" % (reffile)

  fmt = "   Chrom: %5s Tran: %3d -- center moved by %12.4f Ang"

  with open(reffile,'r') as InFile:
    lines = InFile.readlines()
  for line in lines:
    data  = line.split()

    if len(data) < 4 or len(data) > 5:
      c.error("Wrong format for modcent")
    
    chrom = data[0]

    if len(data) == 5: 
      # Look for transition number
      try:
        itran = int(data[1])
      except:
        c.error("Transition number not understood: %s" % data[1])
      NewCent = np.array(data[2:5],dtype=float)
      trlist  = [itran-1]
    else: 
      trlist = range(c.NTran[c.ChromList.index(chrom)])
      NewCent = np.array(data[1:4],dtype=float)

    for j in trlist:
      k = c.TrIdx(chrom,j) 
      if c.v():
        print fmt % (chrom,j+1,np.linalg.norm((NewCent-Cent[k])))
      Cent[k] = NewCent

  return Cent

# *****************************************************************************
#
# Modify electronic couplings
#
def modcoup(system):
  # Note this function is called AFTER seltran
  # and AFTER building of exciton matrix, in order
  # to be able to add intra-site couplings

  exc = system.copy()

  reffile = c.ExtFiles['modcoup']
  c.checkfile(reffile)
  if c.v():
    print " ... Read coupling changes in file %s " % reffile

  with open(reffile) as f:
    while True:
      # strip comments
      line = f.readline().split('#')[0]
      if not line: break
      data = line.split()
      try: 
        iChrom,jChrom,iTran,jTran = data[:4]
        iTran,jTran = int(iTran),int(jTran)        
        newcoup     = float(data[4])
      except ValueError: 
        c.error('Incorrect format in coupling file')
      ii = exc.ChromList.TrIdx(iChrom,iTran) - 1 
      jj = exc.ChromList.TrIdx(jChrom,jTran) - 1
      if ii == jj:
        c.error("Trying to change a site energy in modcoup")
      if c.v():
        print "     Coupling Chrom %5s (%3d) -- %5s (%3d) set to %10.2f cm^-1"\
              % (iChrom,iTran,jChrom,jTran,newcoup)
      exc.H[ii,jj] = newcoup
      exc.H[jj,ii] = newcoup
      
  return exc



# *****************************************************************************
#
# Assign the atomic weigths
#

def anu2awg(In):
  # Return atomic masses 
  Out = []
  for I in In:
    Out.append(c.AtMass[I])
  return Out

# Save the complete geometry in xyz file
def savegeom(anum,xyz):
  if c.v(1): print " ... save geometry to file %s " %  c.OutFiles['xyz'] 
  NAtoms = len(anum)
  OutF = open(c.OutFiles['xyz'],'w')
  OutF.write("%d\n" %  NAtoms)
  OutF.write("Generated by %s\n" %  c.PROGVERS)
  for i in range(NAtoms):
    OutF.write("%3d %10.5f %10.5f %10.5f \n" % tuple([anum[i]]+xyz[i])) 
  return


# *****************************************************************************
#
# Save a .vmd file to visualize dipoles in VMD
# Usage: vmd -e file.vmd
#
def savevisudipo(system,ExcDipo,MagDipo=None):

  C = system.Cent
  Dipo = system.DipoLen
  ChromList = system.ChromList
  NChrom = ChromList.NChrom
  NTran  = ChromList.NTran

  if c.v(1): 
    print " ... save dipole visualization script to %s " %  c.OutFiles['visudipo'] 
    print "     Edit the file following the instructions therein, then type "
    print "     vmd -e %s " % c.OutFiles['visudipo']

  OutFN = c.OutFiles['visudipo']
  with  open(OutFN,'w') as OutF:
    OutF.write(c.VMDHeader())
  
    Scale    = 2.0*c.PhyCon['ToDeb']
    ScaleMag = 2.0
    # Loads molecule
    OutF.write("mol new %s type xyz\n" % c.OutFiles['xyz'])
    # Changes default rep
    OutF.write("mol modstyle    0 0 CPK 0.5 0.2 24.0 24.0 \n")
    OutF.write("mol modmaterial 0 0 EdgyShiny\n")
    # Sets resids
    for i in range(NChrom):
      idx1 = sum(system.NAtom[:i])
      idx2 = idx1 + system.NAtom[i]
      OutF.write("[atomselect top {index %d to %d}] set resid %d\n" % (idx1,idx2,i+1))
  
    # Write down site dipoles
    OutF.write("\n## SITE DIPOLES ...  \n\n")
    
    k = 0
    for i in range(NChrom):
      OutF.write(" \n")
      for j in range(NTran[i]):
        OutF.write(" \n")
        OutF.write("## CHROM %-10s TRAN %3d \n" % (ChromList.Chrom[i],j+1))
        # Default: visualize only electric dipole of first transition
        OutF.write("## Electric Transition Dipole \n")
        if j is not 0: OutF.write('#') 
        OutF.write("graphics 0 color %d; " % (2*j+1) )
        OutF.write("vmd_draw_vector 0 { %10.4f %10.4f %10.4f } { %10.4f %10.4f %10.4f }\n"\
        % (tuple(C[k])+tuple(Dipo[k]*Scale)) )
        if MagDipo is not None:
          OutF.write("## Magnetic Transition Dipole \n")
          OutF.write("#graphics 0 color blue; ") # Default: comment this
          OutF.write("vmd_draw_vector 0 { %10.4f %10.4f %10.4f } { %10.4f %10.4f %10.4f }\n"\
          % (tuple(C[k])+tuple(MagDipo[k]*ScaleMag)) )
        
        k += 1
    #end for

    # Write down excitonic (electric) dipoles 
    totcent = np.average(C,axis=0)
    OutF.write("\n## EXCITONIC DIPOLES ...  \n\n")
  
    # All excitonic transitions

    for k in range(sum(NTran)): 
      OutF.write("## EXC STATE %3d \n" % (k+1))
      OutF.write("#graphics 0 color green; ")
      OutF.write("vmd_draw_vector 0 { %10.4f %10.4f %10.4f } { %10.4f %10.4f %10.4f }\n"\
      % (tuple(totcent)+tuple(ExcDipo[k]*Scale)) )


    # End
    OutF.write("display resetview\n")
  # file closed

  # Save dipoles
  OutFN = c.OutFiles['dipo']
  formt = "%10.4f %10.4f %10.4f"
  np.savetxt(OutFN,Dipo*c.PhyCon['ToDeb'],fmt=formt,header="Dipoles in debye")

  # Save excitonic dipoles
  # TODO

  if MagDipo is not None:
    # Save magnetic dipoles
    OutFN = c.OutFiles['magdipo']
    np.savetxt(OutFN,MagDipo,fmt=formt,header="Magnetic Dipoles")


# *****************************************************************************
#
# Save coupling values and distances between chromophore pairs
# 
#
def savecoupdist(Coup,Cent):
  OutCoup = open(c.ExtFiles['out_coup'],'w')
  n = 0
  for i in range(c.NChrom):
    for j in range(i+1,c.NChrom):
      for k in range(c.NTran[i]):
        for l in range(c.NTran[j]):
          ii = sum(c.NTran[0:i])+k
          jj = sum(c.NTran[0:j])+l
          dist = np.linalg.norm(np.array(Cent[jj])-np.array(Cent[ii]))
          OutCoup.write("%5s (%2d) %5s (%2d)  %10.4f %10.4f\n" % (c.ChromList[i],k+1,c.ChromList[j],l+1,dist,Coup[n]))
          n += 1

  OutCoup.close()


def savediag(energy,coeff,coef2):

  TblCoeff = np.column_stack((energy/c.PhyCon['eV2wn'],energy,coeff))
  TblProb  = np.column_stack((energy/c.PhyCon['eV2wn'],energy,coef2))
  outfile = c.OutFiles['diag']
  with open(outfile,'w') as fd:
    np.savetxt(fd,TblCoeff,fmt="%10.4f ",delimiter='',newline='\n')
    fd.write("\n")
    np.savetxt(fd,TblProb,fmt='%10.4f ',delimiter='',newline='\n')

  pass


# ***********************************************

#
# Save a results.out file
#

def resout(energy,dip,LD,R):
  EeV = energy/c.PhyCon['eV2wn']
  n = len(energy)
  outfile = c.OutFiles['results']
  with open(outfile,'w') as out:
    out.write('# Created with %s\n' % (c.PROGVERS))
    out.write('# Initial command: \n')
    out.write('# '+' '.join(sys.argv)+'\n')
    out.write('# State Energy     mu^2      LD       RotStr \n')
    for i in range(n):
      out.write('%3d %10.4f %10.4f %10.4f %10.4f \n'\
       % (i+1,EeV[i],dip[i],LD[i],R[i]))
    if c.v():
      print "   ... file %s containing the excitonic results has been saved!\n"\
       % outfile
  pass

# ***********************************************


def coupforster(Cent,Dipo,NChrom,NTran):
  Coup = [] ; Kappa = []
  n = 0
  for i in range(NChrom):
    for j in range(i+1,NChrom):
      for k in range(NTran[i]):
        for l in range(NTran[j]):
          ii = sum(NTran[0:i])+k
          jj = sum(NTran[0:j])+l
          Ci = np.array(Cent[ii])
          Cj = np.array(Cent[jj])
          Di = np.array(Dipo[ii])
          Dj = np.array(Dipo[jj])
          tcoup,tkappa = forster(Di,Dj,Ci,Cj)
          Coup.append(tcoup)
          Kappa.append(tkappa)
          n += 1
  return Coup,Kappa

#
# Compute FORSTER Couplings
# ----------------------------------------------
# D1 ... transition dipole 1 (a.u.)
# D2 ... transition dipole 2 (a.u.)
# C1 ... center of D1 (Ang)
# C2 ... center of D2 (Ang)
# ECXOpts['refrind']    -> refraction index
# MainOpts['verbosity'] -> verbosity output
# -----------------------------------------------
#
def forster(D1,D2,C1,C2):
  Ang2au      = c.PhyCon['ToAng']
  hartree2cm  = c.PhyCon['Town']
  s = 1.0/(c.OPT['refrind']**2)
  C1 = C1/Ang2au
  C2 = C2/Ang2au
  dist     = C2-C1
  normdist = np.linalg.norm(dist)
  versdist = dist/normdist
  normD1   = np.linalg.norm(D1)
  normD2   = np.linalg.norm(D2)
  versD1   = D1/normD1
  versD2   = D2/normD2
  fact1    = np.dot(versD1,versD2)
  kappa    = np.dot(versD1,versD2)-3*(np.dot(versD1,versdist)*np.dot(versD2,versdist))
  coupvac  = (kappa * normD1 * normD2) / (normdist**3)
  coup     = coupvac*s

  #if c.OPT['verbosity'] > 2:
  #  print "Kappa         : %8.3f"      % kappa
  #  if ( ECXOpts['refrind'] == 1.0 ):
  #    print "Coup     : %8.3f cm-1" % (coupvac*hartree2cm)
  #  else:
  #    print "Coup vacuo    : %8.3f cm-1" % (coupvac*hartree2cm)
  #    print "Coup (s=%4.2f) : %8.3f cm-1" % (s,coup*hartree2cm)

  return(coup*hartree2cm,kappa)



#
# Print out couplings and center distance
#
def prtcoup(system):
  # Print out couplings larger than a threshold
  Cent  = system.Cent
  Kappa = system.Kappa
  Chrom = system.ChromList.Chrom
  NChrom = system.NChrom
  NTran  = system.NTran
  M = system.H

  Dist   = distance.cdist(Cent,Cent, 'euclidean')
  thresh = c.OPT['CouPrtThr']
  if thresh is False: return
  print "\n Couplings larger than %7.2f cm^-1: " % thresh
  print " Chrom Chrom  Tran  Tran   V (cm-1)  Dist (Ang)        Kappa"
  print " -------------------------------------------------------------"

  OutCoupFile = open(c.OutFiles['coup'],'w')

  MKappa = np.zeros(M.shape)

  if Kappa is not None: 
    L = 0
    for igi in range(NChrom):
      idxi = system.get_chrom_tran(igi)
      nci  = len(idxi)
      for igj in range(igi+1,NChrom):
        idxj = system.get_chrom_tran(igj)
        ncj  = len(idxj)
        nblk = len(idxi)*len(idxj)

        MKappa[np.ix_(idxi,idxj)] = Kappa[L:L+nblk].reshape(nci,ncj)
        L += nblk

  
  for i in range(NChrom):
    for j in range(i,NChrom):
      for k in range(NTran[i]):
        for l in range(NTran[j]):
          kk = sum(NTran[:i])+k
          ll = sum(NTran[:j])+l
          if kk == ll: continue
          if abs(M[kk,ll]) > thresh:
            # Print on screen
            print "%6s %6s  %3d %3d   %10.4f   %10.4f   %10.4f" \
            % (Chrom[i],Chrom[j],k+1,l+1,M[kk,ll],Dist[kk,ll],MKappa[kk,ll])
          # Save all on output file
          OutCoupFile.write("%6s %6s  %3d %3d   %10.4f   %10.4f   %10.4f \n" \
          % (Chrom[i],Chrom[j],k+1,l+1,M[kk,ll],Dist[kk,ll],MKappa[kk,ll]))

  OutCoupFile.close()
  return


def prtsite(system):

  NTran = system.NTran
  Site  = system.Site

  OutFile = open(c.OutFiles['site'],'w')
  z = 0
  for i in range(system.NChrom):
    NTran = system.NTran
    fmt   = "%3d   " + ("%10.4f"*NTran[i])+"\n"
    isite = tuple(Site[z:z+NTran[i]])
    prt   = tuple([i+1])+isite
    OutFile.write(fmt % prt)
    z += NTran[i]
  OutFile.close()
  return


def setldaxis():
  if c.OPT['LDAxis'] not in ('x','y','z'):
    axis = c.OPT['LDAxis'].split(',')
    if len(axis) != 3:
      c.error(' Wrong definition of LD Axis: %s' % c.OPT['LDAxis'])
    c.OPT['LDAxis'] = axis 
    


