# LIBTBX_SET_DISPATCHER_NAME phenix.emringer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

"""
Implementation of the EMRinger method, with plots, for use with the 
EMRinger workflow.

References:
  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387

  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS.
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Manuscript in preparation.
"""

# Any software that wants to use the pkl output of this tool 
# should import ringer_residue and ringer_chi from it.
from __future__ import division
import libtbx.phil
import numpy
from libtbx import easy_pickle
from libtbx import easy_mp
from libtbx.str_utils import make_header
from libtbx.utils import Sorry, Usage
from libtbx import adopt_init_args, Auto
from cStringIO import StringIO
import time
import os
import sys
import sqlite3

master_phil = libtbx.phil.parse("""
pdb_file = None
  .type = path
cif_file = None
  .type = path
  .multiple = True
map_coeffs = None
  .type = path
map_label = None
  .type = str
map_file = None
  .type = path
sampling_method = linear *spline direct
  .type = choice(multi=False)
sampling_angle = 5
  .type = int
grid_spacing = 1./5
  .type = float
scaling = *sigma volume
  .type = choice(multi=False)
skip_alt_confs = True
  .type = bool
nproc = 1
  .type = int
show_gui = False
  .type = bool
output_base = None
  .type = str
""")

def run (args, out=None, verbose=True) :
  t0 = time.time()
  if (out is None) : out = sys.stdout
  from iotbx import file_reader
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="pdb_file",
    reflection_file_def="map_coeffs",
    map_file_def="map_file",
    cif_file_def="cif_file",
    usage_string="""\
phenix.emringer model.pdb map.mrc [cif_file ...] [options]

%s
""" % __doc__)
  params = cmdline.work.extract()
  validate_params(params)
  pdb_in = cmdline.get_file(params.pdb_file)
  pdb_in.check_file_type("pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  map_coeffs = ccp4_map = None
  if (params.map_coeffs is not None) :
    mtz_in = cmdline.get_file(params.map_coeffs)
    mtz_in.check_file_type("hkl")
    best_guess = None
    best_labels = []
    all_labels = []
    for array in mtz_in.file_server.miller_arrays :
      if (array.info().label_string() == params.map_label) :
        map_coeffs = array
        break
      elif (params.map_label is None) :
        if (array.is_complex_array()) :
          labels = array.info().label_string()
          all_labels.append(labels)
          if (labels.startswith("2FOFCWT") or labels.startswith("2mFoDFc") or
              labels.startswith("FWT")) :
            best_guess = array
            best_labels.append(labels)
    if (map_coeffs is None) :
      if (len(all_labels) == 0) :
        raise Sorry("No valid (pre-weighted) map coefficients found in file.")
      elif (best_guess is None) :
        raise Sorry("Couldn't automatically determine appropriate map labels. "+
          "Choices:\n  %s" % "  \n".join(all_labels))
      elif (len(best_labels) > 1) :
        raise Sorry("Multiple appropriate map coefficients found in file. "+
          "Choices:\n  %s" % "\n  ".join(best_labels))
      map_coeffs = best_guess
      print >> out, "  Guessing %s for input map coefficients" % best_labels[0]
  else :
    ccp4_map_in = cmdline.get_file(params.map_file)
    ccp4_map_in.check_file_type("ccp4_map")
    ccp4_map = ccp4_map_in.file_object
  make_header("Iterating over residues", out=out)
  t1 = time.time()
  from mmtbx.ringer import iterate_over_residues
  results = iterate_over_residues(
    pdb_hierarchy=hierarchy,
    map_coeffs=map_coeffs,
    ccp4_map=ccp4_map,
    params=params,
    log=out).results
  t2 = time.time()
  if (verbose) :
    print >> out, "Time excluding I/O: %8.1fs" % (t2 - t1)
    print >> out, "Overall runtime:    %8.1fs" % (t2 - t0)
  if (params.output_base is None) :
    pdb_base = os.path.basename(params.pdb_file)
    params.output_base = os.path.splitext(pdb_base)[0] + "_ringer"
  easy_pickle.dump("%s.pkl" % params.output_base, results)
  print >> out, "Wrote %s.pkl" % params.output_base
  csv = "\n".join([ r.format_csv() for r in results ])
  open("%s.csv" % params.output_base, "w").write(csv)
  print >> out, "Wrote %s.csv" % params.output_base
  print >> out, "\nReferences:"
  print >> out, """\
  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387"""
  if (params.show_gui) :
    run_app(results)
  else :
    return results

def validate_params (params) :
  if (params.pdb_file is None) :
    raise Sorry("No PDB file supplied (parameter: pdb_file)")
  if (params.map_coeffs is None) and (params.map_file is None) :
    raise Sorry("No map coefficients supplied (parameter: map_coeffs)")
  return True

########################################################################
# GUI
from wxtbx import plots
import wx

def run_app (results) :
  app = wx.App(0)
  frame = RingerFrame(None, -1, "Ringer results")
  frame.show_results(results)
  frame.Show()
  app.MainLoop()

class RingerFrame (plots.plot_frame) :
  def create_plot_panel (self) :
    plot = RingerPlot(self, figure_size=(6,8))
    plot.canvas.Bind(wx.EVT_CHAR, self.OnChar)
    return plot

  def draw_top_panel (self) :
    self.top_panel = wx.Panel(self, style=wx.SUNKEN_BORDER)
    panel_szr = wx.BoxSizer(wx.VERTICAL)
    self.top_panel.SetSizer(panel_szr)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    panel_szr.Add(szr2)
    txt1 = wx.StaticText(self.top_panel, -1, "Residue to display:")
    szr2.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.chooser = wx.Choice(self.top_panel, -1, size=(200,-1))
    szr2.Add(self.chooser, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnSelect, self.chooser)
    self.Bind(wx.EVT_CHAR, self.OnChar)
    self.chooser.Bind(wx.EVT_CHAR, self.OnChar)
    return self.top_panel

  def OnSelect (self, event) :
    selection = event.GetEventObject().GetSelection()
    self.plot_panel.show_residue(self.results[selection])

  def show_results (self, results) :
    self.results = results
    choices = [ result.format() for result in results ]
    self.chooser.SetItems(choices)
    self.chooser.SetSelection(0)
    self.plot_panel.show_residue(self.results[0])

  def OnChar (self, event) :
    key = event.GetKeyCode()
    if (len(self.results) == 0) : return
    selection = self.chooser.GetSelection()
    if (key in [wx.WXK_TAB, wx.WXK_RETURN, wx.WXK_SPACE]) :
      if (selection < (len(self.results) - 1)) :
        selection += 1
      elif (len(self.results) > 0) :
        selection = 0
    elif (key in [wx.WXK_DELETE, wx.WXK_BACK]) :
      if (selection > 0) :
        selection -= 1
      else :
        selection = len(results) - 1
    self.chooser.SetSelection(selection)
    self.plot_panel.show_residue(self.results[selection])

class RingerPlot (plots.plot_container) :
  def show_residue (self, residue) :
    if (self.disabled) : return
    self.figure.clear()
    subplots = []
    for i in range(1, residue.n_chi + 1) :
      chi = residue.get_angle(i)
      if (chi is None) : continue
      if (len(subplots) > 0) :
        p = self.figure.add_subplot(4, 1, i, sharex=subplots[0])
      else :
        p = self.figure.add_subplot(4, 1, i)
        p.set_title(residue.format())
      p.set_position([0.15, 0.725 - 0.225*(i-1), 0.8, 0.225])
      x = [ k*chi.sampling for k in range(len(chi.densities)) ]
      p.plot(x, chi.densities, 'r-', linewidth=1)
      p.axvline(chi.angle_current, color='b', linewidth=2, linestyle='--')
      p.axvline(chi.peakchi, color='g', linewidth=2, linestyle = '--')
      p.axhline(0, color=(0.4,0.4,0.4), linestyle='--', linewidth=1)
      p.axhspan(0.3,1,facecolor="green",alpha=0.5)
      p.axhspan(-1,0.3,facecolor="grey",alpha=0.5)
      p.set_xlim(0,360)
      ax = p.get_axes()
      ax.set_ylabel("Rho")
      ax.set_xlabel("Chi%d" % i)
      subplots.append(p)
    for p in subplots[:-1] :
      for label in p.get_axes().get_xticklabels() :
        label.set_visible(False)
    p.text(0,-0.5,'Green = Peak, Blue = Modelled',
        transform=ax.transAxes)
    self.canvas.draw()
    self.canvas.Fit()
    self.Layout()
    self.parent.Refresh()

if (__name__ == "__main__") :
  run(sys.argv[1:])
