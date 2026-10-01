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
# excsystem.py Exciton objects module
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

from __future__ import print_function
import numpy  as np
from collections import OrderedDict

# Import COMMON Modules
import common as c


class ChromTranList(OrderedDict):
    
    def __init__(self,*args,**kwargs):
        
        super(ChromTranList, self).__init__(*args)
        for k,v in self.items():
            try: 
                iter(v)
                # Copy
                self[k] = [x for x in v]
            except:
                self[k] = [v]

        for k,v in list(self.items()):
            if not str(k) == k:
                self[str(k)] = v
                del self[k]

    
    @classmethod
    def from_crfile(cls,infile,assume_tran=False):
        c.checkfile(infile)
        with open(infile) as f:
            # Strip comments
            lines = [x.split('#')[0] for x in f.readlines()]
        
        ChromList = []; Tran = []

        for line in lines:
            data = line.split()
            if not data: continue
            ChromList.append(data[0])
            ITran = list(map(int,data[1:]))
            if not ITran and assume_tran:
                # Assume one transition per chromophore
                ITran = [1]
            Tran.append(ITran)
        
        return cls(list(zip(ChromList,Tran)))
        
    def set_NTran(self,NTran):
       "Sets the transitions for every chromophore" 

       for i,k in enumerate(self.keys()):
          self[k] = [ j+1 for j in range(NTran[i]) ]

    def TrIdx(self,chrom,itran):
        """Return the ID of the transition for 
           chromophore chrom and transition itran"""

        # Cache NTran once: rebuilt from scratch on every access otherwise
        NTran = self.NTran

        if not str(chrom) == chrom:
            chrom = str(chrom)

        if not chrom in list(self.keys()):
            return None
        ichrom = self.index(chrom)
        if itran > NTran[ichrom]:
            return None
        return sum(NTran[:ichrom]) + itran

    def kChromTran(self,k):
        """ Return the chrom,tran pair of the k-th transition
        """
        if k < 1: 
            return None,None

        NTran = self.NTran

        if k > sum(NTran): 
            return None,None
        
        for i in range(self.NChrom):
            itran = k - sum(NTran[:i])
            chrom = self.Chrom[i]
            if itran <= NTran[i]:
                break
        return chrom,itran
    
    @property
    def NChrom(self):
        return len(self)
    
    @property
    def NTran(self):
        return [len(v) for v in list(self.values()) ]
    
    @property
    def Chrom(self):
        "Returns list with chromophore names/IDs"
        return list(self.keys())
        
    def index(self,x):
        "Returns index of chromophore given name"
        y = str(x)
        return list(self.keys()).index(y)

    def append(self,chrom):
        "Adds a new chromophore w. zero transitions"
 
        ch = str(chrom)
        if ch in self.keys():
            raise ValueError( 'Chromophore %s already in ChromList!' % (ch))

        self[ch] = []

    def add_tran(self,chrom):
        "Adds a transition to the chromophore chrom"

        Ch = self[chrom]
        if len(Ch) > 0:
          JTr = self[chrom][-1] + 1
        else:
          JTr = 1

        self[chrom].append(JTr)
          

    def copy(self):
        "Return deep copy of the object"
        cls = type(self)
        return cls(self.items())

    def __repr__(self):
        cname  = 'excsystem.'+self.__class__.__name__
        params = repr(list(self.items()))
         
        return '%s(%s)' % (cname,params)

    def __str__(self):
        s  = ' ChromList: \n'
        s += ' Chrom     Tran\n'
        for k,v in self.items():
            s += ' %5s     ' % (k) 
            s +=  ' '.join('%4d' % (i) for i in v )
            s += '\n'
        return s
    

class ExcSystem(object):
    """
    Class for an Excitonic System. 
    """

    
    def __init__(self,chromlist,site=None,coup=None,Cent=None,
                 DipoLen=None,DipoVel=None,Mag=None,Kappa=None):
        
        if isinstance(chromlist,ChromTranList):
            self.ChromList = chromlist
        else:
            self.ChromList = ChromTranList(chromlist)

        self.data_init()
            
        if site is not None:
            self.site = site.copy()
        if coup is not None:
            self.coup = coup.copy()
        if Cent is not None:
            self.Cent = Cent.copy()
        if DipoLen is not None:
            self.DipoLen = DipoLen.copy()
        if DipoVel is not None:
            self.DipoVel = DipoVel.copy()
        if Mag is not None:
            self.Mag     = Mag.copy()
        if Kappa is not None:
            self.Kappa   = Kappa.copy()
        else:
            self.Kappa = None

    def data_init(self):

        self.xyz     = None
        self.anum    = None
        self.NAtom   = None
        self.site    = None
        self.coup    = None
        self.Cent    = None
        self.DipoLen = None
        self.DipoVel = None
        self.Mag     = None
    
    @property
    def NTran(self):
        return self.ChromList.NTran
    
    @property
    def NChrom(self):
        return self.ChromList.NChrom
    
    @property
    def Site(self):
        if not self.has_Hamiltonian:
            return self.site
        else:
            return self.H.diagonal()
    
    @property
    def Coup(self):
        if not self.has_Hamiltonian:
            return self.coup

        rows,cols = _pair_order_indices(self.NTran)
        return self.H[rows,cols]

    @property
    def has_Hamiltonian(self):
        try:
            self.H
        except AttributeError:
            return False
        return self.H is not None

    @property
    def has_geom(self):
        try:
            self.xyz,self.anum,self.NAtom
        except AttributeError:
            return False
        return (self.xyz is not None) \
                and (self.anum is not None) and (self.NAtom is not None)

    @property
    def has_eig(self):
        try:
            self.coeff,self.energy
        except AttributeError:
            return False
        return (self.coeff is not None) and (self.energy is not None)

    def add_geom(self,anum,xyz,NAtom):

        nat = len(xyz)

        assert len(anum)  == nat
        assert sum(NAtom) == nat

        self.xyz  = np.array(xyz)
        self.anum = np.array(anum)
        self.NAtom = np.array(NAtom)


    def get_chrom_geom(self,ichrom):

        try: 
            iter(ichrom)
        except:
            ichrom = [ichrom]
        
        idx = []
        for ich in ichrom:
          id1 = sum(self.NAtom[:ich])
          id2 = id1 + self.NAtom[ich]
          idx.extend(list(range(id1,id2)))

        return self.anum[idx],self.xyz[idx]

    def get_chrom_tran(self,ichrom):

        id1 = sum(self.NTran[:ichrom])
        id2 = id1 + self.NTran[ichrom]

        return np.arange(id1,id2)

    def update_sitecoup(self):
        self.coup = self.Coup
        self.site = self.Site

    def buildmatrix(self):
        
        coup = self.coup
        site = self.site

        # Cache NTran/NChrom once: ChromTranList.NTran rebuilds it
        # at every access
        NTran  = self.NTran
        NChrom = self.NChrom

        # Check the dimensions
        dimen = sum(NTran)
        ncoup = 0
        for i in range(NChrom):
          for j in range(i+1,NChrom):
            ncoup += NTran[i]*NTran[j]

        lcoup = len(self.coup)
      
        if lcoup != ncoup :
          print("Couplings Found     : %4d" % lcoup)  
          print("Couplings Requested : %4d" % ncoup)  
          c.error("Confused in the Dimension!","matrixbuilder")
      
        if c.v(1):
          print(" ... Matrix dimension       : %4d" % dimen) 
          print(" ... Number of chromophores : %4d" % NChrom) 
          print(" ... Number of COUPLINGS    : %4d" % ncoup) 
      
      
        # Set all matrix elements to zero:
        self.H = np.zeros((dimen,dimen))
      
        # Write the diagonal part
        self.H[np.diag_indices_from(self.H)] = site

        # Write the off-diagonal blocks:
        # offsets[i] is the starting index of chromophore i's transitions
        offsets = np.concatenate(([0],np.cumsum(NTran)))
        L=0
        for igi in range(NChrom):
          idxi = np.arange(offsets[igi],offsets[igi+1])
          nci  = len(idxi)
          for igj in range(igi+1,NChrom):
            idxj = np.arange(offsets[igj],offsets[igj+1])
            ncj  = len(idxj)
            nblk = len(idxi)*len(idxj)
            
            self.H[np.ix_(idxi,idxj)] = coup[L:L+nblk].reshape(nci,ncj)
            L += nblk

        # Symmetrize Hamiltonian
        idLT = np.tril_indices(dimen,-1)
        self.H[idLT] = self.H.T[idLT]

        self.H[self.H == 0.0] = 0.0 # Fix "real" zeros for printing

        # Just to be sure, delete coeff/energy

        try: 
            del self.coeff,self.energy
        except: 
            pass

        return self.H

    def savematrix(self,outfile):
        if not self.has_Hamiltonian: 
            return
        np.savetxt(outfile,self.H,fmt='%10.1f',delimiter='')

    def copy(self):    
        " Try to return a deep copy of the object "
        cls = type(self)
        crlist_tmp = self.ChromList.copy()

        if self.has_Hamiltonian:
            # H is copied directly below; do not derive site/coup from it
            # here (self.Coup is expensive) only to discard them right away.
            new = cls(crlist_tmp,None,None,
                    self.Cent,self.DipoLen,self.DipoVel,self.Mag,
                    self.Kappa)
            new.H = self.H.copy()
        else:
            # No Hamiltonian yet: site/coup ARE the data, 
            # so carry them over directly.
            new = cls(crlist_tmp,self.site,self.coup,
                    self.Cent,self.DipoLen,self.DipoVel,self.Mag,
                    self.Kappa)

        if self.has_geom:
            new.xyz  = self.xyz.copy()
            new.anum = self.anum.copy()
            new.NAtom = self.NAtom[:]
        return new

    def seltran(self,chromlist):
        """
        Returns a new instance of ExcSystem with transitions selected according to 
        the passed chromlist object
        """


        if isinstance(chromlist,ChromTranList):
            selchromlist = chromlist
        else:
            selchromlist = ChromTranList(chromlist)

        # Cache Chrom/NTran once
        ChromNames = self.ChromList.Chrom
        NTran      = self.NTran
        name_to_idx = {name:i for i,name in enumerate(ChromNames)}
        offsets = np.concatenate(([0],np.cumsum(NTran)))

        IndChrom = [ name_to_idx[str(x)] for x in selchromlist ]    


        mask = np.zeros(sum(NTran),dtype=bool)

        for I in IndChrom:
          Chrom = ChromNames[I]
          idxI  = np.arange(offsets[I],offsets[I+1])
          # beware of off-by-one
          idxT = [j for i,j in enumerate(idxI) if i+1 in selchromlist[Chrom] ]
          mask[idxT] = True

        new = self.copy()

        new.H    = new.H[mask][:,mask]
        new.Mag  = new.Mag[mask]
        new.Cent = new.Cent[mask]
        new.DipoLen = new.DipoLen[mask]
        new.DipoVel = new.DipoVel[mask]

        if new.Kappa is not None:
          new.Kappa   = new.Kappa[get_coupmask(self.ChromList,mask)]

        new.ChromList = selchromlist
        new.update_sitecoup()

        if new.has_geom:
          new.anum,new.xyz = new.get_chrom_geom(IndChrom) 
          new.NAtom = new.NAtom[IndChrom]

        return new


    def diagonalize(self):
        """
        Compute eigenvalues and eigenvectors of the excitonic matrix.
        The Eigenvalues corresponds to the Excitonic Energies and the eigenvectors
        are the contributes of the old state i to the new excitonic state k.
        All the values in the excitonic matrix must be expressed in cm-1.
        
         Excitonic Matrix (M)          eigenval                 eigenvec
        -----------------------      ------------       -------------------------
             E1 V12 V13                   E1'            c1(E1') c1(E2') c1(E2')
            V12  E2 V23                   E2'            c2(E1') c2(E2') c2(E2')
            V13 V23  E3                   E3'            c3(E1') c3(E2') c3(E2')
        
        Eigenvalues are in crescent order of energies, as well as the corresponding
        eigenvectors. An external formatted file is also saved containing the
        eigenvalue vector, eginvector matrix.
        """

        
        if not self.has_Hamiltonian:
            self.buildmatrix()

        eigval,eigvec = eighsort(self.H)

        self.coeff = eigvec.T
        self.coef2 = self.coeff**2

        self.energy = eigval

        return self.energy,self.coeff

    @classmethod
    def from_matrix(cls,M):
        
        dim = M.shape[0]
        crlist = ( (i+1,1) for i in range(dim) )

        new   = cls(crlist)
        new.H = M
        new.update_sitecoup()

        return new


    def __str__(self):
        s  = ' Exciton System with %d chromophores\n' \
                % (self.NChrom)
        s += ' Transitions per chromophore: \n %s \n\n' \
                % (str(self.NTran))

        s += ' Has geometry:    %s \n' % (self.has_geom)
        s += ' Has Hamiltonian: %s \n' % (self.has_Hamiltonian)
        s += ' Is diagonalized: %s \n' % (self.has_eig)

        s += '\n'

        return s

    def save(self,outfile):
        savedict = self.__dict__
        np.savez(outfile,**savedict)

##########################

def _pair_order_indices(NTran):
    """
    Row/column transition indices (i,j) for every cross-chromophore pair I<J,
    in chromophore-pair-major then row-major-within-pair order -- the same
    order used by Coup/Kappa and by buildmatrix()'s off-diagonal blocks
    (which matches how per-pair coupling files are read).
    """
    NTran = np.asarray(NTran)
    M = int(NTran.sum())

    # Chromophore of each transition
    chrom_of = np.repeat(np.arange(len(NTran)),NTran)

    rows,cols = np.triu_indices(M,k=1)

    # Exclude (i,j) when chrom[i] == chrom[j]
    cross = chrom_of[rows] != chrom_of[cols]
    rows,cols = rows[cross],cols[cross]

    # Sort by chromophore pair first (chrom[rows], then chrom[cols]),
    # then by transition within that pair (rows, then cols)
    order = np.lexsort((cols,rows,chrom_of[cols],chrom_of[rows]))
    return rows[order],cols[order]

def get_coupmask(chromlist,mask):

    rows,cols = _pair_order_indices(chromlist.NTran)
    return mask[rows] & mask[cols]

def load_npz(infile):

    data = np.load(infile,allow_pickle=True)

    chromlist = data['ChromList'].item()

    obj  = ExcSystem(chromlist)

    if 'H' in data and data['H'].tolist() is not None:
        obj.H = data['H']
    elif 'site' in data:
        obj.site = data['site']
        if 'coup' in data:
            obj.coup = data['coup']

    obj.DipoLen = data['DipoLen']
    obj.DipoVel = data['DipoVel']
    obj.Mag     = data['Mag']
    obj.Cent    = data['Cent']
    if 'Kappa' in data:
        obj.Kappa   = data['Kappa']
        if obj.Kappa.tolist() is None:
          obj.Kappa = None
    obj.coup    = obj.Coup
    obj.site    = obj.Site
    
    try:
        obj.anum    = data['anum']
        obj.xyz     = data['xyz']
        obj.NAtom   = data['NAtom']
        if obj.anum.tolist() is None:
          obj.anum = None
        if obj.xyz.tolist() is None:
          obj.xyz = None

    except KeyError:
        pass


    return obj


def eighsort(matrix):
    "Diagonalize with ordered eigenvalues "

    e,v = np.linalg.eigh(matrix)
    index = e.argsort()
    e = e[index]
    v = v[:,index]
    return e,v

