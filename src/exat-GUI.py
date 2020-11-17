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
# exat-GUI.py MAIN PROGRAM GUI
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

import os,sys
import gtk,pango

import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt

from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as Canvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar


import numpy as np

import common as c
import util as u
import readdata      # Read data (from Gaussian output)
#import buildmatrix   # Build the excitonic matrix
#import diag          # Diagonalize the excitonic matrix
import trans         # Compute the excitonic dip and R

# Reset program version in common
GUIversion = '0.2.1'
PROGVERS = c.VERSION + '\n GUI version %s ' % (GUIversion)

c.OutFilesBK = c.OutFiles.copy()

def resetcommon():
  # Clear all variables in common 
  c.NChrom = 2
  c.NTran  = []
  c.ChromList = []
  c.CoupList = [] # List of indices k->ij (i,j chromophores)
  c.NAtom = []    # List of atom numbers, for each chromophore
  c.anum  = []    # List of atom numbers
  c.xyz   = []    # Geometry

  resetoutfiles()
  
  pass

def resetoutfiles():
  # Reset out files
  for k in c.OutFiles:
    c.OutFiles[k] = '/dev/stdout'


def guessversion(inlog):
  return 'g16'

def guessexternal(crlist):
  # Guess if there are external files for readexternal
  search = ('site.in','coup.in','dipo.in','cent.in')
  flist  = os.listdir(os.getcwd())
  clname = os.path.split(crlist)[1]
  if 'chromlist.in' in clname:
    base = clname.split('.chromlist.in')[0]
  else:
    base = clname.split('.in')[0]
  # Look for BASE.site.in etc.
  if all([base+'.'+x in flist for x in search]):
    BN = base+'.'
  # Look for site.in etc.
  elif all([x in flist for x in search]):
    BN = ''
  else:
    return False # Found nothing
  # Possibly set external files
  c.ExtFiles['incoup'] = BN+'coup.in'
  c.ExtFiles['incent'] = BN+'cent.in'
  c.ExtFiles['indipo'] = BN+'dipo.in'
  c.ExtFiles['insite'] = BN+'site.in'
  return True

# Spectrum Calculation
def specalc(E,D2,R,broad,shp='gau',xlims=None):

  factOD = 108.8983364615 # w*|mu|^2 to epsilon
  factCD = 4.351846416E-2 # w*R      to Delta(epsilon)

  # xlims = [xmin,xmax,xstep]

  EMin = E.min()
  EMax = E.max()

  if xlims is None:
    xlims = [0,0,0]
    # Compute auto xmin and xmax 
    if shp == 'gau':
      xlims[0] = EMin-4*broad[0]
      xlims[1] = EMax+4*broad[-1]
    else:
      xlims[0] = EMin-8*broad[0]
      xlims[1] = EMax+8*broad[-1]
    xlims[2] = 2.0 #temp

  # Define X-Axis (in wavenumbers)
  x=np.arange(*xlims,dtype=float)
  NStep,NTran = len(x),len(E)
	
  # Initilize the intensity
  OD = np.zeros((NStep,NTran),dtype=float)
  CD = np.zeros((NStep,NTran),dtype=float)
# LD = np.zeros((NStep,NTran),dtype=float)


  if shp == "gau": # Gaussian
    broad = broad/np.sqrt(2*np.log(2)) # HWHH to Sigma
    LineShape = lambda p, x: (x/(p[1]*np.sqrt(2*np.pi))*np.exp(-0.5*((x-p[0])/p[1])**2))
  else :           # Lorentzian
    LineShape = lambda p, x: (x/(np.pi*p[1]))*(p[1]**2/((x-p[0])**2+p[1]**2))
		    
  for j in range(len(E)):
    p = (E[j],broad[j]) 
    OD[:,j] = LineShape(p,x)*D2[j]*factOD
#   LD[:,j] = LineShape(p,x)*dipold[j]*factOD
    CD[:,j] = LineShape(p,x)*R[j]*factCD

  return x,OD,CD


# Custom Nav Toolbar for MatPlotLib
class MyNavToolbar(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]


# Key binding for matplotlib plots
# Restore broken keybindings 
def mpl_keypress(event):
# print 'press', event.key
# sys.stdout.flush()
  if event.key in ('s','ctrl+s'):
    event.canvas.toolbar.save_figure()
  elif event.key in ('h','r','home'):
    event.canvas.toolbar.home()  
  elif event.key in ('o'):
    event.canvas.toolbar.zoom()  
  elif event.key in ('p'):
    event.canvas.toolbar.pan()  
  elif event.key in ('left','c','backspace'):
    event.canvas.toolbar.back()  
  elif event.key in ('v','right'):
    event.canvas.toolbar.forward()  
  elif event.key in ('q','ctrl+w'):
    win = event.canvas.toolbar.get_toplevel()
    win.hide()


# List for List widget
class EasyList:
    def __init__(self, treeview, columns, coltypes):
        """EasyList(gtktreeview, titles, types) -> EasyList object
        
        Creates a new instance of the class EasyList that
        handles the content of a gtktreeview easier than
        the gtk method.
        
        Arguments:
            treeview : The gtk treeview widget to be used.
            columns : A tuple of column titles ("title1", "title2", ...)
            coltypes : Types of column (str, str, int...)
        
        Comments:
            Everything is treated as fundamental types in EasyList."""
                
        if len(columns) != len(coltypes) or len(columns) == 0:
            raise ValueError("columns and coltypes should be equal in size.")
        
        self.treeview = treeview
        self.store = gtk.ListStore(*coltypes)
        self.cellr = gtk.CellRendererText()
        self.cellr.set_property('cell-background', '#efefef')
        self.cellr.set_property('alignment', pango.ALIGN_RIGHT)
        font = pango.FontDescription('courier bold 11')
        self.cellr.set_property('font-desc', font)
        self.cellr.set_property('editable', True)

        for count,column in enumerate(columns):
            current = gtk.TreeViewColumn(column, self.cellr)
            current.add_attribute(self.cellr, "text", count)
            current.set_sort_column_id(count) 
            current.set_resizable(True)
            self.treeview.append_column(current)

        
        self.treeview.set_model(self.store)
    
    def append_row(self, data):
        if type(data) is not list and type(data) is not tuple:
            data = (data,)
        self.store.append(data)

    def clear(self):
      self.store.clear()

    def set_selection_mode(self,mode):
      self.treeview.get_selection().set_mode(mode)

    # Highlight rows where col == data
    def select_data(self,col,data):
      sel = self.treeview.get_selection()
      sel.unselect_all()
      for i,row in enumerate(self.store):
         if data == row[col]:
           sel.select_iter(row.iter)
    # END


class EXATGUI:

  original_data = None
  specopt       = dict()

  def generic_error(self,message=None):
    dl = gtk.MessageDialog(parent=self.window,flags=0,type=gtk.MESSAGE_ERROR,
                       buttons=gtk.BUTTONS_OK,message_format=None)
    if message is not None:
      dl.set_markup('<big>'+message+'</big>')
    else:
      dl.set_markup("Generic error.")
    dl.run()
    dl.destroy()

  def generic_info(self,message):
    dl = gtk.MessageDialog(parent=self.window,flags=0,type=gtk.MESSAGE_INFO,
                       buttons=gtk.BUTTONS_OK,message_format=None)
    dl.set_markup(message)
    dl.run()
    dl.destroy()

  def clearexat(self,setversion=None):
    self.inlog = None
    c.OPT['read'] = 'g16'
    if setversion is not None:
      c.OPT['read'] = setversion
    c.OPT['verbosity'] = 1
    #c.OPT['OutPrefix'] = 'EXAT'
  

    # Possibly delete exat vars
    try:
      del self.Site,self.Coup,self.Mag,self.DipoLen,self.MagInt,\
          self.Cent,self.energy,self.coeff,self.H,self.DipoVel,self.Kappa
    except: pass
    resetcommon()
    self.clearwindow()


  def clearwindow(self):
    # Clear status bar and title
    self.window.set_title('EXAT - EXcitonic Analysis Tool')
    #self.statusbar.remove_all(self.context_id)
    #self.statusbar.remove_all(self.context_id) # Doing twice should do the trick
    # Clear tables
    self.sitelist.clear()
    self.exclist.clear()
    try:
      self.sellist.clear()
      del self.sellist
    except:pass

  def exat_read(self):
    self.clearwindow()
    if self.inlog[-4:] == '.npz':
      c.ExtFiles['load'] = self.inlog
      c.OPT['read'] = 'load'
    elif self.inlog[-3:] != '.in':
      c.OPT['logfile'] = self.inlog
      c.OPT['read'] = guessversion(self.inlog)
    else:
      c.ExtFiles['crlist'] = self.inlog
      os.chdir(os.path.dirname(os.path.abspath(self.inlog)))
      if guessexternal(self.inlog):
        # Ext files are set in function
        c.OPT['read'] = 'external'
        c.OPT['seltran'] = True
      else:
	raise Exception('Could not read ext files')
    # System-independent call (including seltran)
    logname = self.inlog.split('/')[-1]
    self.statusbar.push(self.context_id, "Reading %s ..." % (logname))
    try:
      self.system,self.Sel = readdata.Read()
    except:
      raise Exception('Could not read file %s' % (self.inlog))
    self.statusbar.push(self.context_id, "Loaded %s" % (logname))
    self.window.set_title('EXAT - %s' % (logname))
    self.update_data()
    # Also, save original data in such a way that we can retrieve them
    self.original_data = self.system.copy()
    
  def exat_run(self):
    if c.OPT['RCalc'] == 'mag' and c.OPT['read'] == 'external':
      self.generic_error('No full electric-magnetic treatment with external files!')
      return False

    self.statusbar.push(self.context_id, "Building excitonic Hamiltonian ... ")
    self.H = self.system.buildmatrix()

    self.statusbar.push(self.context_id, "Diagonalize excitonic Hamiltonian ... ")
    self.energy,self.coeff = self.system.diagonalize()

    # Convert some quantities to A.U.
    EEN  =  self.energy/c.PhyCon['Town']  # Excitonic Energies (Hartree)

    self.statusbar.push(self.context_id, "Compute excitonic properties ... ")
    RxDel  = np.cross(self.system.Cent/c.PhyCon['ToAng'],self.system.DipoVel)    
    # Compute internal magnetic moment (DipoMag is gauge-dependent)
    self.MagInt = self.Mag - RxDel

    self.EXCDipoLen = trans.EXCalc(self.coeff,self.DipoLen)
    
    # Compute Linear Absorption Spectrum
    print "\n ... Compute the Linear Absorption Spectrum"
    self.EXCDipo2   = np.sum(self.EXCDipoLen**2,axis=1)
    
    # Compute Linear Dichroism Spectrum
    print "\n ... Compute the Linear Dichroism Spectrum"
    self.LD = trans.LinDichro(self.EXCDipoLen)
    
    # Compute Rotational Strength ...
    print "\n ... Compute the Circular Dichroism Spectrum"
    self.EXCRot = trans.RotStrength(EEN,self.Cent,self.coeff,
              self.DipoLen,self.EXCDipoLen,self.DipoVel,self.MagInt,RxDel,self.Site)
    # Done
    self.statusbar.push(self.context_id, 
       "Ready. Showing %s " % (self.inlog.split('/')[-1]))

    self.update_data()
    self.update_exc_data()



  def seltran(self,selection):
    # Custom seltran that does not read a chromlist file

    if not np.any(selection):
      # No transition selected
      self.generic_error('No transition selected!')
      return False

    # Get selection
    SelChromList = [] ; Sel = []
    for i in range(self.system.NChrom):
      if np.any(selection[i]):
        # Only append Chroms that have selected transitions
        SelChromList.append(i+1)
        ITran  = np.where(selection[i])[0]
        Sel.append(ITran+1)

    self.system = self.system.seltran( zip(SelChromList,Sel) )

    self.update_data()

    # All OK
    return True 


  def savedata(self,outfile):
    c.OPT['OutPrefix'] = outfile
    c.OutFiles = c.OutFilesBK.copy()
    c.setoutfiles()

    # npz
    self.system.save(c.OutFiles['exatdata'])

    # matrix
    np.savetxt(c.OutFiles['matrix'],self.H,fmt='%10.1f',delimiter='',newline='\n')

    #diag
    prob = np.array(self.coeff)**2
    TblCoeff = np.column_stack((self.energy/c.PhyCon['eV2wn'],self.energy,self.coeff))
    TblProb  = np.column_stack((self.energy/c.PhyCon['eV2wn'],self.energy,prob))
    diagfile = c.OutFiles['diag']
    with open(diagfile,'w') as f:
      np.savetxt(f,TblCoeff,fmt="%10.4f ",delimiter='',newline='\n')
      f.write("\n")
      np.savetxt(f,TblProb,fmt='%10.4f ',delimiter='',newline='\n')

    # Visudipo
    u.savegeom(self.system.anum,self.system.xyz.tolist())
    u.savevisudipo(self.system,self.EXCDipoLen,-self.MagInt)

    #results
    u.resout(self.energy,self.EXCDipo2,self.LD,self.EXCRot)

    # Spectrum
    try:     lshape = self.specopt['lshape']
    except:  lshape = 'gau'
    broad = [self.obj('adjustbroad').get_value()]*len(self.energy)
    D2 = self.EXCDipo2*c.PhyCon['ToDeb']**2
    w,OD,CD = specalc(self.energy,D2,self.EXCRot,broad,lshape)
    np.savetxt(outfile+'.OD.dat',np.column_stack((w,OD.sum(axis=1),OD)))
    np.savetxt(outfile+'.CD.dat',np.column_stack((w,CD.sum(axis=1),CD)))

    resetoutfiles()
    return True

  def on_window1_destroy(self,object,data=None):
    gtk.main_quit()

  # Avoid deleting "Levels window"
  def on_windowlevels_delete_event(self,widget,data=None):
    self.winlevels.hide()
    return True

  def on_windowcoeff_delete_event(self,widget,data=None):
    self.wincoeff.hide()
    return True
  
  def on_windowspec_delete_event(self,widget,data=None):
    self.winspec.hide()
    return True

  def on_gtk_quit_activate(self,menuitem,data=None):
    gtk.main_quit()

  def on_gtk_about_activate(self, menuitem, data=None):
    self.response = self.aboutdialog.run()
    self.aboutdialog.hide()

  def on_menu_seltran_activate(self,menuitem,data=None):
    # first of all, retrieve original data
    try:
      self.system = self.original_data.copy()
      self.system.buildmatrix()
      self.update_data()
    except:
      return None
    self.init_winseltran()
    self.response = self.winseltran.run()
    if self.response == 1: 
      select = []
      for x in self.sellist:
        select.append([y for y in x][1:])
      # Select transitions
      self.seltran(select)
      self.exat_run()
    self.winseltran.hide()

  def on_menu_mag_toggled(self,menuitem,data=None):
    # Set RCalc options
    c.OPT['RCalc'] = 'mag' if menuitem.get_active() else 'mu'
    # If possible, re-do EXAT
    try: self.Site
    except: return 
    self.exat_run()
    return

  def on_buttonlevels_clicked(self,object,data=None):
    self.winlevels.show()
    if not self.plot_levels():
      self.winlevels.hide()

  def on_buttoncoeff_clicked(self,object,data=None):
    self.wincoeff.show()
    if not self.plot_coeff_matrix():
      self.wincoeff.hide()

  def on_buttonspec_clicked(self,object,data=None):
    self.winspec.show()
    if not self.plot_spectrum():
      self.winspec.hide()

  def on_gtk_close_activate(self, menuitem, data=None):
    self.clearexat()

  def on_notebook1_switch_page(self,notebook,page,page_num,data=None):
    self.current_tab = notebook.get_nth_page(page_num)
    pass

  def on_units_selector_changed(self,widget,data=None):
    active = widget.get_active()
    item = self.specunilist[active][1]
    self.specopt['units'] = item
    try: self.energy
    except: return
    self.plot_spectrum()
    pass

  def on_lshape_selector_changed(self,widget,data=None):
    active = widget.get_active()
    item = self.lshapelist[active][1]
    self.specopt['lshape'] = item.lower()[:3]
    try: self.energy
    except: return
    self.spectrum_plot_update()
    pass

  def on_gtk_save_as_activate(self, menuitem, data=None):
    self.fcd = gtk.FileChooserDialog("Save as...",None,gtk.FILE_CHOOSER_ACTION_SAVE,
          buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
    if self.current_folder is not None:
      self.fcd.set_current_folder(self.current_folder)

    self.response = self.fcd.run() # Open
    if self.response == gtk.RESPONSE_OK:
      outfile = self.fcd.get_filename()
      # if a file was choosen save the current folder
      self.current_folder = self.fcd.get_current_folder()
      self.fcd.destroy()
      try: 
        self.savedata(outfile)
        self.generic_info('Data saved to %s' % outfile)
      except IOError, e: 
        self.generic_error('Could not save file %s\n %s' % (outfile,e))
      self.curr_outfile = outfile # save outfile info
    else: self.fcd.destroy()
    mpl.rcParams['savefig.directory'] = self.current_folder
    pass

  def on_gtk_save_activate(self, menuitem, data=None):
    try: self.curr_outfile
    except:
      self.on_gtk_save_as_activate(None,None)
      return
    self.savedata(self.curr_outfile)
    pass

  def update_data(self):
    self.sitelist.clear()
    # Update site energy
    self.Site    = self.system.Site/c.PhyCon['eV2wn']
    self.DipoLen = self.system.DipoLen
    self.DipoVel = self.system.DipoVel
    self.Cent    = self.system.Cent
    self.Mag     = self.system.Mag 

    k = 0
    NChrom  = self.system.NChrom
    NTran   = self.system.NTran
    crlist  = self.system.ChromList.Chrom
    for IChr in range(NChrom):
      for ITr in range(NTran[IChr]):
        # Format numbers here
	D2 = np.sum(self.DipoLen[k]**2)*c.PhyCon['ToDeb']**2 #D^2
        self.sitelist.append_row([IChr+1,ITr+1,("%8.4f" % self.Site[k]),("%12.4f"%D2)])
        k += 1
    pass

  def update_exc_data(self):
    self.exclist.clear()
    # Update exciton states
    EeV = self.energy/c.PhyCon['eV2wn']
    n = len(EeV)
    for J in range(n):
      # Format numbers here
      D2 = self.EXCDipo2[J]*c.PhyCon['ToDeb']**2 #D^2
      col2,col3,col4 = ("%8.4f" % EeV[J],"%12.4f" % D2,"%12.4f" % self.EXCRot[J])
      self.exclist.append_row([J+1,col2,col3,col4])

    # Update H 
    HView = self.builder.get_object("textview1")
    HView.modify_font(pango.FontDescription("monospace 11"))
    buff = gtk.TextBuffer()
    HView.set_buffer(buff)

    H_str = ''
    for J in range(n):
      H_str += (' %12.4f'*n + '\n') % tuple(self.H[J].tolist())
    buff.set_text(H_str)

    pass


  def on_file_open_activate(self, menuitem, data=None):
    self.fcd = gtk.FileChooserDialog("Open...",None,gtk.FILE_CHOOSER_ACTION_OPEN,
          buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
    if self.current_folder is not None:
      self.fcd.set_current_folder(self.current_folder)

    # Set file filters
    filtr = gtk.FileFilter()
    filtr.set_name("Gaussian output files (*.log,*.out)")
    filtr.add_pattern("*.log")
    filtr.add_pattern("*.out")
    self.fcd.add_filter(filtr)

    filtr = gtk.FileFilter()
    filtr.set_name("Chrom list files (*.in)")
    filtr.add_pattern("*.in")
    self.fcd.add_filter(filtr)

    filtr = gtk.FileFilter()
    filtr.set_name("All files")
    filtr.add_pattern("*")
    self.fcd.add_filter(filtr)
   
    self.response = self.fcd.run() # Open
    if self.response == gtk.RESPONSE_OK:
      self.inlog = self.fcd.get_filename()
      # if a file was choosen save the current folder
      self.current_folder = self.fcd.get_current_folder()
      self.fcd.destroy()
      self.exat_read()
      self.exat_run()
    else: self.fcd.destroy()


  def init_treeviews(self):
    self.sitelist = EasyList(self.builder.get_object("treeview1"),
          ("Chromophore", "Transition","Site Energy (eV)", "Dipole Strength (D^2)"), (int, int,str,str))

    self.exclist = EasyList(self.builder.get_object("treeview2"),
          ("Exc. State", "Energy (eV)","Dipole Strength (D^2)","Rotatory Strength (10^-40 esu^2 cm^2)"), (int, str, str, str))


    self.sitelist.set_selection_mode(gtk.SELECTION_MULTIPLE)
    self.exclist.set_selection_mode(gtk.SELECTION_MULTIPLE)
    pass 

  # Connect signals to actions
  def connect_more(self):
    self.obj('wlevelsclose').connect('clicked',lambda w: self.winlevels.hide())
    self.obj('adjustlevels').connect("value_changed",lambda w: self.levels_plot_update())
    self.obj('wcoeffclose').connect('clicked',lambda w: self.wincoeff.hide())
    self.obj('wspecclose').connect('clicked',lambda w: self.winspec.hide())
    self.obj('adjustbroad').connect("value_changed",lambda w: self.spectrum_plot_update())
    self.obj('adjustshift').connect("value_changed",lambda w: self.spectrum_plot_update())

    self.obj('seltran_selall').connect('clicked', lambda w: self.on_seltran_alltoggle(True))
    self.obj('seltran_selnone').connect('clicked',lambda w: self.on_seltran_alltoggle(False))

  # Shortcut
  def obj(self,name):
    return self.builder.get_object(name)


  # Initialize seltran window
  def init_winseltran(self):
    # Do not init if no exat data
    try: self.Site
    except: return False
      
    NTrMax = max(self.system.NTran)

    coltyp  = [int]+NTrMax*[bool]    

    try:
      self.sellist
    except:
      # Build custom store/view
      self.sellist = gtk.ListStore(*coltyp)
      tview        = self.obj("treeview_seltran")
      tview.set_model(self.sellist)

      for col in tview.get_columns(): 
        tview.remove_column(col)
  
      # Add cols
      rtext    = gtk.CellRendererText()
      ctext    = gtk.TreeViewColumn("Chrom",rtext)
      tview.append_column(ctext)
      ctext.add_attribute(rtext, "text",0)
      for I in range(NTrMax):
        rtogg  = gtk.CellRendererToggle()
        rtogg.connect('toggled',self.on_seltran_cell_toggled,I+1)
        ctogg  = gtk.TreeViewColumn(str(I+1),rtogg)
        tview.append_column(ctogg)
        ctogg.add_attribute(rtogg, "active",I+1)
      
      # Add rows
      #self.sellist.clear()
      NTran = self.system.NTran
      for J in range(self.system.NChrom):
        self.sellist.append([J+1]+[True]*NTran[J]+[False]*(NTrMax-NTran[J]))


  def on_seltran_cell_toggled(self,widg,path,col):
    self.sellist[path][col] = not self.sellist[path][col]

  def on_seltran_alltoggle(self,value=None):
    value = bool(value)
    N = len(self.sellist)
    M = len(self.sellist[0])
    for I in range(N):
      for J in range(1,M):
        self.sellist[I][J] = value
    pass

  # Initialize stuff needed for plots
  def init_plots(self): 
    # Levels figure 
    self.levelsfigure = plt.figure()
    self.levelsfigure.patch.set_facecolor('#dfdfdf') # Color around plot
    self.levelscanvas = Canvas(self.levelsfigure)
    self.levels_box.pack_start(self.levelscanvas, True, True)
    self.levToolbar = MyNavToolbar(self.levelscanvas, self.winlevels)
    self.levels_box.pack_start(self.levToolbar,False,True,0)
    self.levToolbar.pan() # Activate pan by default
    self.levToolbar.show()
    self.levelscanvas.set_size_request(400, 600)
    self.levelscanvas.set_flags(gtk.HAS_FOCUS | gtk.CAN_FOCUS)
    self.levelscanvas.grab_focus()
    # Connect events and key bindings
    self.levelscanvas.mpl_connect('pick_event', self.onpick_levels)
    self.levelscanvas.mpl_connect('key_press_event',mpl_keypress)

    # Coeff figure
    self.coeffigure = plt.figure()
    self.coeffigure.patch.set_facecolor('#dfdfdf') 
    self.coefcanvas = Canvas(self.coeffigure)
    self.coeff_box.pack_start(self.coefcanvas)
    self.coefcanvas.set_size_request(600, 400)
    self.coefcanvas.set_flags(gtk.HAS_FOCUS | gtk.CAN_FOCUS)
    self.coefcanvas.grab_focus()

    # Spectrum
    self.specfigure = plt.figure() 
    self.specfigure.patch.set_facecolor('#dfdfdf')
    self.speccanvas = Canvas(self.specfigure)
    self.spec_box.pack_start(self.speccanvas)
    self.speccanvas.set_size_request(600, 400)
    self.speccanvas.set_flags(gtk.HAS_FOCUS | gtk.CAN_FOCUS)
    self.speccanvas.grab_focus()
    # toolbar
    self.specToolbar = MyNavToolbar(self.speccanvas, self.winspec)
    self.spec_box.pack_start(self.specToolbar,False,True,0)
    self.specToolbar.pan() # Activate pan by default
    self.specToolbar.show()
    # Connect events and key bindings
    self.speccanvas.mpl_connect('key_press_event',mpl_keypress)
    #self.speccanvas.connect()
    # Unit selector
    self.specunilist = gtk.ListStore(int,str)
    self.specunilist.append([0,'cm^-1'])
    self.specunilist.append([1,'eV'])
    self.specunilist.append([2,'nm'])
    self.spec_unisel.set_model(self.specunilist)
    CRT = gtk.CellRendererText()
    self.spec_unisel.pack_start(CRT,True)
    self.spec_unisel.add_attribute(CRT,'text',1)
    self.spec_unisel.set_active(0)
    # Lineshape selector
    self.lshapelist  = gtk.ListStore(int,str)
    self.lshapelist.append([0,'gaussian'])
    self.lshapelist.append([1,'lorentzian'])
    self.spec_lshsel.set_model(self.lshapelist)
    self.spec_lshsel.pack_start(CRT,True)
    self.spec_lshsel.add_attribute(CRT,'text',1)
    self.spec_lshsel.set_active(0)


    return

  # Plot exciton levels
  def plot_levels(self):

    self.statusbar.push(self.context_id, "Plotting levels ...")

    self.levelsfigure.clf()

    # Hide x-axis and set reasonable range
    ax = self.levelsfigure.gca()
    ax.set_ylabel('Energy (eV)')
    #ax.axes.get_xaxis().set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(axis='x',which='both',bottom='off',top='off',pad=15)
    ax.xaxis.set_ticks([1.15,1.85])
    ax.xaxis.set_ticklabels(['Sites','Excitons'])
    bfont = mpl.font_manager.FontProperties(size=14, weight='bold')
    for lab in ax.xaxis.get_ticklabels():
      lab.set_fontproperties(bfont)

    #ax.axes.autoscale(False)
    ax.set_xlim([0.8,2.2])
    # Margin
    plt.subplots_adjust(left=0.2,right=0.9)

    if not self.levels_plot_update():
      return False

    self.statusbar.push(self.context_id, "Plot Ready")
    return True


  def levels_plot_update(self):

    # Clear subplot
    ax = self.levelsfigure.gca()

    # Possibly clear data
    ax.lines = []

    # See if we have exat results
    try: self.energy
    except:
      self.generic_error('No excitonic levels present')
      return False
    # Plot site energy levels
    for I,E in enumerate(self.Site):
      ax.plot([1,1.3],[E,E],'b-',linewidth=2.0,picker=5)  

    # Plot excitonic levels
    C2     = self.coeff**2
    thresh = self.obj('adjustlevels').get_value()/1e2
    for I,E in enumerate(self.energy/c.PhyCon['eV2wn']):
      ax.plot([1.7,2],[E,E],'b-',linewidth=2.0,picker=5)  
      contribs = np.where(C2[I] > thresh)[0]
      for J in contribs:
        ax.plot([1.3,1.7],[self.Site[J],E],'k--')

    # Draw
    self.levelscanvas.draw()
    self.levelscanvas.show()
    return True 


  def plot_coeff_matrix(self):
    try: self.coeff
    except: 
      self.generic_error('No excitonic data')
      return False
    self.coeffigure.clf()

    ax = self.coeffigure.gca()

    pc = ax.pcolor(self.coeff,cmap=plt.get_cmap('seismic'))
    self.coeffigure.colorbar(pc)

    self.coefcanvas.draw()
    self.coefcanvas.show()
    return True


  def plot_spectrum(self):
    # Plot
    self.specfigure.clf()
    # Set subplots
    ODax = self.specfigure.add_subplot(211)
    CDax = self.specfigure.add_subplot(212,sharex=ODax)
    self.specfigure.tight_layout()
    plt.subplots_adjust(left=0.2,right=0.95,hspace=0.3)
    if not self.spectrum_plot_update():
      return False
    return True



  def spectrum_plot_update(self):
    try:
      E     = self.energy
      mu2   = self.EXCDipo2*c.PhyCon['ToDeb']**2
      rotst = self.EXCRot
    except:
      self.generic_error('No excitonic data')
      return False

    # Compute spectra
    ODsticks  = mu2*E   # intensity = mu^2*energy
    CDsticks  = rotst*E # intensity = R*energy

    broad = [self.obj('adjustbroad').get_value()]*len(E)
    # Spectrum shift in cm^-1
    shift = self.obj('adjustshift').get_value()
    if shift != 0.0: E = E+shift
   
    try:     lshape = self.specopt['lshape']
    except:  lshape = 'gau'
    w,self.OD,self.CD = specalc(E,mu2,rotst,broad,lshape)

    # Set x-axis units
    w      = self.spectrum_get_units(w)
    wstick = self.spectrum_get_units(E)

    # Sum contributions 
    ODtot = self.OD.sum(axis=1)
    CDtot = self.CD.sum(axis=1)

    MaxOD = ODtot.max()
    MaxCD = CDtot.max()
    MinCD = CDtot.min()
    CD_absmax = max(MaxCD,-MinCD)


    # Normalize sticks
    ODsticks *= MaxOD/ODsticks.max()*0.9
    if (CDsticks != 0).any():
      CDsticks *= CD_absmax/np.abs(CDsticks).max()*0.9

    # Plot:
    ODax,CDax = self.specfigure.get_axes()
    # Reset
    ODax.lines = [];  ODax.collections = []; 
    CDax.lines = [];  CDax.collections = []; 
    #OD
    ODax.set_ylabel('Epsilon',fontweight='bold',fontsize=16)
    ODax.set_ylim(0,MaxOD*1.1)
    ODax.set_title('Absorption Spectrum')
    ODax.plot(w,ODtot,linewidth=2.5,linestyle="-",color="red")
    ODax.vlines(wstick,ODsticks,0,linewidth=2.0,color='blue')

    #CD
    CDax.set_ylabel('Delta Epsilon',fontweight='bold',fontsize=16)
    ODax.set_ylim(MinCD*1.1,MaxOD*1.1)
    CDax.set_title('Circular Dichroism Spectrum')
    CDax.axhline(linewidth=1.0,linestyle="-",color="black")
    CDax.plot(w,CDtot,linewidth=2.5,linestyle="-",color="red")
    CDax.vlines(wstick,CDsticks,0,linewidth=2.0,color='blue')

    self.speccanvas.draw()
    self.speccanvas.show()
    return True



  def spectrum_get_units(self,wcm):
    try:     u = self.specopt['units']
    except:  u = 'cm'
    if u == 'eV':      return wcm/c.PhyCon['eV2wn']
    elif u == 'nm':    return 1e7/wcm
    else:              return wcm

  def onpick_levels(self,event):
    A =  event.artist.get_ydata()[0]
    X =  event.artist.get_xdata()[0]
    msg = "Energy: %7.4f eV %7.0f cm^-1" % (A,A*c.PhyCon['eV2wn']) 
    self.statusbar.push(self.context_id,msg)
    print msg
    
    # Select site or exc state based on X-pos
    if X < 1.5:
      self.sitelist.select_data(2,("%8.4f" % A))
    else:
      self.exclist.select_data(1,("%8.4f" % A))

    # Maybe do more interesting things?
    pass

  def __init__(self):
    # Init GUI data from glade
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    self.gladefile = os.path.join(scriptdir,"exat-GUI.glade")
    self.builder = gtk.Builder()
    self.builder.add_from_file(self.gladefile)
    # Connect signals defined in .glade file
    self.builder.connect_signals(self)


    # Get some object shortcuts
    self.window      = self.obj("window1")
    self.aboutdialog = self.obj("aboutdialog1")
    self.aboutdialog.set_version(PROGVERS)
    self.winlevels   = self.obj('windowlevels')
    self.wincoeff    = self.obj('windowcoeff')
    self.winspec     = self.obj('windowspec')
    self.winseltran  = self.obj('dialog_seltran')
    self.statusbar   = self.obj("statusbar1")
    self.context_id  = self.statusbar.get_context_id("status")
    self.notebook    = self.obj('notebook1')
    self.levels_box  = self.obj('levels_box')
    self.coeff_box   = self.obj('coeff_box')
    self.spec_box    = self.obj('spec_box')
    self.spec_unisel = self.obj('units_selector')
    self.spec_lshsel = self.obj('lshape_selector')

    # Connect signals
    self.connect_more()

    # Set current folder 
    self.current_folder = os.path.expanduser(os.getcwd())

    # Init various GUI items
    self.init_treeviews()

    # Init plot items
    self.init_plots()

    # Finally show window
    self.window.show()
    self.clearexat()

    #temp: Possibly open command-line log file
    if (len(sys.argv) > 1):
      self.inlog = sys.argv[1]
      try:
	self.inlog = os.path.abspath(self.inlog)
        self.exat_read()
        self.exat_run()
      except:
        self.generic_error('Could not read %s' % (self.inlog))

if __name__ == "__main__":
  main = EXATGUI()
  gtk.main()

