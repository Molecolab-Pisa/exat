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
# common.py COMMON MODULE
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
import sys, os 
import numpy as np

# VERSION
VERSION  = "2.0.0a"
PROGVERS = "Exat - EXcitonic Analysis Tool - Version %s" % VERSION

# ******************************************************************************
#
# Fill OPT dictionary with default values. If you want to change the default
# options you can modify here.
#
# ******************************************************************************

OPT =    {
            'workdir'   : os.getcwd(),      # Specify the working directory
            'verbosity' : 0 ,               # Verbosity of the output on screen
            'RCalc'     : 'mu',             # Select type of calculation for R
            'Cent'      : 'geom',           # How to compute the center of the chromophore
            'external'  : False,            # Read all data from external files
            'read'      : 'g16',            # Select Gaussian version for EET coupling 
            'logfile'   : None,             # Gaussian LogFile
            'CleanCoup' : 0.0,              # Treshold for couplings
            'ScaleCoup' : 1.0,              # Scaling factor for coupling
            'ScaleDipo' : 1.0,              # Scaling factor for dipoles
            'reorient'  : None,             # Axis to reorient transition dipole moments (length)
            'forcedipo' : False,            # Force dipo to be aligned with --reorient axis
            'anadipo'   : None,             # Request dipole orientation analysis
            'ScaleTran' : None,             # Request dipole and coupling scaling
            'ModSite'   : False,            # Whether to modify site energies
            'ModDipoLen': False,            # Whether to modify transition electric dipoles
            'ModDipoMag': False,            # Whether to modify transition magnetic dipoles
            'ModCoup'   : False,            # Whether to modify couplings
            'ModCent'   : False,            # Whether to modify centers
            'OutPrefix' : False,
            'LDAxis'    : "z",              # LD Axis
            'ctype'     : 'geom' ,          # How define the center of the chromophore
            'coup'      : 'trden',          # Type of coupling calculation
            'refrind'   : 1.0,              # Refraction index for computing screening in point-dipole coupling calculation
            'seltran'   : False,            # Select only particular transitions.
            'CouPrtThr' : False,            # Threshold for coupling printout
            'savetprod' : False,            # Save triple product of (mu) CD
            'savecoeff' : False             # Save hi-acc coefficients in numpy file
          }

# ******************************************************************************
#
# Physical constants and conversion factors
#
# ******************************************************************************

PhyCon =   {
            'Planck'    : 6.62606957E-34  , # Plank constant ( J * sec)
            'Slight'    : 2.99792458E10   , # Speed of ligth (m/s)
            'Avog'      : 6.02214129E23   , # Avogadro number
            'ToAng'     : 0.52917721092   , # Angstroms per Bohr
            'ToKG'      : 1.660538921E-27 , # Kilograms per AMU
            'ToE'       : 1.602176565E-19 , # ESU per electron charge (Coloumb*electron)
            'Town'      : 220000          , # cm-1 per hartree
            'ToDeb'     : 2.54158         , # Debye per electroncharge per Bohr
            'Hartre'    : 4.35974434E-18  , # Joules per Hartre
            'eV2wn'     : 8065.5446811132 , # eV per wavenumber 
            'ToeV'      : 27.211396132    , # eV per Hartree
            'EMKG'      : 4.35974434E-18*10000  # ????
           }


# ******************************************************************************
#
# Atomic Masses
#
# ******************************************************************************

AtMass =  {    1 :    1.0079400,
               2 :    4.0026020,
               3 :    6.9410000,
               4 :    9.0121820,
               5 :   10.8110000,
               6 :   12.0107000,
               7 :   14.0067000,
               8 :   15.9994000,
               9 :   18.9984032,
              10 :   20.1797000,
              11 :   22.9897700,
              12 :   24.3050000,
              13 :   26.9815380,
              14 :   28.0855000,
              15 :   30.9737610,
              16 :   32.0650000,
              17 :   35.4530000,
              18 :   39.9480000,
              19 :   39.0983000,
              20 :   40.0780000,
              21 :   44.9559100,
              22 :   47.8670000,
              23 :   50.9415000,
              24 :   51.9961000,
              25 :   54.9380490,
              26 :   55.8450000,
              27 :   58.9332000,
              28 :   58.6934000,
              29 :   63.5460000,
              30 :   65.4090000,
              31 :   69.7230000,
              32 :   72.6400000,
              33 :   74.9216000,
              34 :   78.9600000,
              35 :   79.9040000,
              36 :   83.7980000,
              37 :   85.4678000,
              38 :   87.6200000,
              39 :   88.9058500,
              40 :   91.2240000,
              41 :   92.9063800,
              42 :   95.9400000,
              43 :   97.9072160,
              44 :  101.0700000,
              45 :  102.9055000,
              46 :  106.4200000,
              47 :  107.8682000,
              48 :  112.4110000,
              49 :  114.8180000,
              50 :  118.7100000,
              51 :  121.7600000,
              52 :  127.6000000,
              53 :  126.9044700,
              54 :  131.2930000,
              55 :  132.9054500,
              56 :  137.3270000,
              57 :  138.9055000,
              58 :  140.1160000,
              59 :  140.9076500,
              60 :  144.2400000,
              61 :  144.9127440,
              62 :  150.3600000,
              63 :  151.9640000,
              64 :  157.2500000,
              65 :  158.9253400,
              66 :  162.5000000,
              67 :  164.9303200,
              68 :  167.2590000,
              69 :  168.9342100,
              70 :  173.0400000,
              71 :  174.9670000,
              72 :  178.4900000,
              73 :  180.9479000,
              74 :  183.8400000,
              75 :  186.2070000,
              76 :  190.2300000,
              77 :  192.2170000,
              78 :  195.0780000,
              79 :  196.9665500,
              80 :  200.5900000,
              81 :  204.3833000,
              82 :  207.2000000,
              83 :  208.9803800,
              84 :  208.9824160,
              85 :  209.9871000,
              86 :  222.0176000,
              87 :  223.0197307,
              88 :  226.0254030,
              89 :  227.0277470,
              90 :  232.0381000,
              91 :  231.0358800,
              92 :  238.0289100,
              93 :  237.0481670,
              94 :  244.0641980,
              95 :  243.0613730,
              96 :  247.0703470,
              97 :  247.0702990,
              98 :  251.0795800,
              99 :  252.0829700,
             100 :  257.0950990,
             101 :  258.0984250,
             102 :  259.1010200,
             103 :  262.1096900,
             104 :  261.1087500,
             105 :  262.1141500,
             106 :  266.1219300,
             107 :  264.1247300,
             108 :  269.1341100,
             109 :  268.1388200
          }

# ******************************************************************************
#
# Standard file name
#
# ******************************************************************************

ExtFiles = {
   # Input Files
   'incoup'       : 'coup.in',        # Coupling values (cm-1)
   'insite'       : 'site.in',        # Site energies (eV)
   'incent'       : 'cent.in',        # Postion of the centers (Ang)
   'dipo'         : 'dipo.in',        # Electric transtion dipoles (lenght formulation) in Debye
   'magdipo'      : 'dipomag.in',     # Magnetic transtion dipoles in a.u.
   'modcoup'      : 'coup.in',        # Coupling modification file
   'crlist'       : 'chromlist.in',   # List of chromophores and selected tranitions
   'refaxis'      : 'reference.in',   # List of chromophores, refaxes and angles (for select tr)
   'ScaleTran'    : 'scaletran.in'    # List of chromophores and scaling factors (for select tr)
   }   

OutFiles = {
   # Output Files
   'results'      : 'results.out',    # Excitonic energies, mu^2,LD,CD
   'matrix'       : 'matrix.dat',     # Excitonic matrix 
   'diag'         : 'diag.dat',       # Excitonic coefficients and squared coefficients
   'coeff'        : 'coeff.npy',      # Excitonic coefficients (numpy file) 
   'tprod'        : 'tprod.dat',      # Triple product matrix 
   'xyz'          : 'geometry.xyz',   # Complete geometry
   'visudipo'     : 'visudipo.vmd',   # VMD script to visualize tr dipoles
   'dipo'         : 'dipo.out',       # Site dipoles
   'magdipo'      : 'magdipo.out',    # Magnetic dipoles
   'rstrength'    : 'component.out',  # Rotational strength components
   'coup'         : 'coup.out',       # List of couplings
   'site'         : 'site.out',       # List of site energies
   'exatdata'     : 'exat.npz',       # Exat data in npz binary format
   # What?
   'dipolen'      : 'dipolen.out',    # List of transition dipole moments (length)
   'dipovel'      : 'dipovel.out',    # List of transition dipole moments (velocity)
   'cent'         : 'cent.out',
   }

# ******************************************************************************
#
# Print the welcome message
#
# ******************************************************************************

def welcome():
   print '                                                                      '
   print '                                        mm                            '
   print '                                     mMMm                             '
   print '                                   mMMMMm         m                   '
   print '                                  mMMMMm          mMm                 '
   print '                                  mMMMMm          mMm                 '
   print '                                  mMMMMMm        mMMm                 '
   print '                                  MMMMMMMMMMMMMMMMMMm                 '
   print '                                 mMMMMMMMMMMMMMMMMMm                  '
   print '        __  ___      __    ____________      __MMMm     __            '
   print '       /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_           '
   print '      / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \          '
   print '     / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /          '
   print '    /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/           '
   print '            /_  __/ __ \/ __ \/ /  / ___/                             '
   print '             / / / / / / / / / /   \__ \                              '
   print '            / / / /_/ / /_/ / /___ __/ /                              '
   print '           /_/  \____/\____/_____/____/                               '
   print '             mMMMMMMMMMMMMMMMm                                        '
   print '           mMMMMMMMMMMMMMMMm                                          '
   print '         mMMMMMMMMMMMMMMMMM   + ------------------------------------ +'
   print '        mMMMMMMMMMMMMMMMMm    |             E  X  A  T               |'
   print '       mMMMMMMMMMMMMMMMMMm    + ------------------------------------ +'
   print '       mMMMMm       mMMMMMm   | S. Jurinovich, L. Cupellini          |'
   print '       mMMMm       mMMMMMMm   | C.A. Guido                           |'
   print '        mMm       mMMMMMMm    |                            ver %-6.5s|' % VERSION
   print '         m       mMMMMMMm     |              molecolab.dcci.unipi.it |'
   print '                mMMMMMm       + ------------------------------------ +'
   print '                                                                      '


# ******************************************************************************
#
# Standard error:
#
# ******************************************************************************

def error(string,where=None):
  # Determine calling function
  if where is None: where = sys._getframe().f_back.f_code.co_name
  msg = "\n------ ERROR in %s ------\n%s\n" % (where,string)
  raise AssertionError(msg)


# ******************************************************************************
#
# Debug message:
#
# ******************************************************************************

def debug():
  # Determine calling function
  where = sys._getframe().f_back.f_code.co_name
  print("\n------ EXIT FOR DEBUGGING -------\nI am in: %s\n\n" % where)
  sys.exit()


# ******************************************************************************
#
# Progression Bar
#
# ******************************************************************************

def ShowProgressBar(istep,ntot):
  percentage = np.round(float(istep)/float(ntot)*100)
  sys.stdout.write('\r')
  sys.stdout.write("[%-50s] %d%%" % ('='*int(percentage/2), percentage))
  sys.stdout.flush()


# ******************************************************************************
#
# Print lists by using a specifc format
#
#   List   : input list to print on screen (list of elements of the same type)
#   NER    : maximum number of element to be print in a row (integer)
#   fmt    : format specification (string)
#   margin : (optional) number of white spaces from left margin (integer)
#   header : (optional) text header (string)
#
# ******************************************************************************

def printout(List,NER,fmt,margin=2,header=''):
  try:
    N     = len(List)          # Total number of elements in List
  except:
    N     = 1
  if N > 1:
    NER   = int(NER)           # Number of element per row
    NFR   = N/NER              # Number of full row
    NLR   = N-NFR*NER          # Number of element in the last row
    if ( N >= NER):
      Fmt = (" "*margin+(fmt)*NER+"\n")*(N/NER)+(" "*margin+(fmt)*NLR+"\n")
    else:
      Fmt = (" "*margin+(fmt)*N+"\n")
    print "\n"+" "*margin+header+"\n"
    print Fmt % tuple(List)
  else:
    print "\n"+" "*margin+header+"\n"
    print fmt % List
  pass


# ******************************************************************************
#
# StringConverter
#
# ******************************************************************************

def stringconverter(In):

  Out = []
  blocks = In.split(',')
  for block in blocks:
    blk = block.split('-')
    if len(blk) > 1:
      sel =  map(int,blk)
      Out += range(sel[0],sel[1]+1,1)
    else:
      Out += [int(block)]
  return Out


# ******************************************************************************
#
# Check if the file exists:
#
# ******************************************************************************

def checkfile(FILENAME):
  if (os.path.isfile(FILENAME) == False):
    print("\n File %s not found!\n" % FILENAME)
    sys.exit()


# ******************************************************************************
#
# Set prefix to output files
#
# ******************************************************************************
def setoutfiles():
  if OPT['OutPrefix'] is not False: 
    for k in OutFiles:
      OutFiles[k] = OPT['OutPrefix']+'.'+OutFiles[k]


# ******************************************************************************
#
# Print out a dictionary or a list (OPT, etc) 
#
# ******************************************************************************

def PrintHeader(title=None):
  length = 70
  sep    = " "+length*'-'
  if title is not None:
    # Find out where to put the title 
    fblank = length/2 - len(title)/2
    fmt     = sep+'\n'+' '*fblank+'%s\n'+sep
    print fmt % (title)
  else:  print sep
  pass

def PrintDict(Stuff,title):
  PrintHeader(title)
  if type(Stuff) is dict:
    for key,val in Stuff.items():  print "   %-14s : %-10s" % (key,val)
  elif type(Stuff) is list:
    for key,val in enumerate(Stuff):  print "   %14d : %-10s" % (key,val)
  pass
  PrintHeader(None)

def PrintOpt():
  PrintDict(OPT,'Selected options:')
  pass

# Return if to print the message!
def v(level=-1):
  v = OPT['verbosity']
  if v is None: v = 0
  return v > level


# ******************************************************************************
#
# VMD Header
#
# ******************************************************************************

def VMDHeader():
  s  = '## Generated with %s\n\n' % PROGVERS 
  s += '## VMD script to view transition dipoles... follow the instructions\n'
  s += '## below the VMD functions to select the dipole(s) to draw         \n\n\n'
  s += '## General Settings                                                \n'
  s += 'menu main on                                                       \n'
  s += 'display projection orthographic                                    \n'
  s += '                                                                   \n'
  s += '## VMD functions to draw a vector                                  \n'
  s += 'proc vmd_draw_arrow {mol start end} {                              \n'
  s += '  set length [veclength [vecsub $end $start]]                      \n'   
  s += '  set conelen [expr max(0.4,0.2*$length) ]                         \n'   
  s += '  set scale [expr max(0.5,(1.0-$conelen/$length))]                 \n'   
  s += '                                                                   \n'   
  s += '  set middle [vecadd $start [vecscale $scale [vecsub $end $start]]]\n'   
  s += '  # You can change the radius of the cylinder/cone                 \n'   
  s += '  graphics $mol cylinder $start $middle radius 0.05 resolution 24  \n'   
  s += '  puts [list cone $middle $end radius 0.15]                        \n'   
  s += '  graphics $mol cone $middle $end radius 0.15 resolution 24        \n'   
  s += '}                                                                  \n'
  s += '                                                                   \n'
  s += 'proc vmd_draw_vector { mol pos val } {                             \n'
  s += '    # Change +1 with any value to scale all vectors!               \n'   
  s += '    set end   [ vecadd $pos [ vecscale +1 $val ] ]                 \n'
  s += '    vmd_draw_arrow $mol $pos $end                                  \n'
  s += '}                                                                  \n\n\n'
  # Instructions
  s += '## Under the comment CHROM X TRAN J there is the command to draw the transition dipole.  \n'
  s += '## The first command prints the electric dipole, the second prints the magnetic dipole.  \n'
  s += '## The electric first transition of each chromophore is displayed by default\n '
  s += '## Below the site transition dipoles, there are also the excitonic (electric) dipoles\n '
  s += '## Uncomment the ones you want to display before running vmd.  \n\n\n'
  return s



# --- END of MODULE ---
