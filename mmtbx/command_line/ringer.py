# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

"""
Implementation of the Ringer method, with plots

Reference:
  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387
"""

from __future__ import absolute_import, division, print_function
from mmtbx.ringer import * # this is deliberate!
import libtbx.phil
from libtbx import easy_pickle
from libtbx.str_utils import make_header
from libtbx.utils import Sorry
from libtbx import runtime_utils
import mmtbx.model
from iotbx.map_model_manager import map_model_manager
import time
import os
import sys
from six.moves import range

master_phil = libtbx.phil.parse("""
model = None
  .type = path
  .style = file_type:pdb bold input_file
cif_file = None
  .type = path
  .multiple = True
  .short_caption = Restraints
  .style = file_type:cif bold input_file
map_coeffs = None
  .type = path
  .short_caption = Map coefficients
  .style = file_type:hkl bold input_file OnChange:extract_ringer_map_labels
map_label = 2FOFCWT,PH2FOFCWT
  .type = str
  .input_size = 200
  .short_caption = 2Fo-FC map labels
  .style = bold renderer:draw_map_arrays_widget noauto
difference_map_label = FOFCWT,PHFOFCWT
  .type = str
  .input_size = 200
  .short_caption = Fo-Fc map labels
  .style = bold renderer:draw_map_arrays_widget noauto
map_file = None
  .type = path
include scope mmtbx.ringer.ringer_phil_str
sampling_method = linear *spline direct
  .type = choice(multi=False)
grid_spacing = 1./5
  .type = float
gui = False
  .type = bool
  .style = hidden
ignore_symmetry_conflicts = False
  .type = bool
  .help = Allow running with PDB file with symmetry that does not match map
output_base = None
  .type = str
output_dir = None
  .type = path
  .style = output_dir
include scope libtbx.phil.interface.tracking_params
""", process_includes=True)
master_params = master_phil

def run(args, out=None, verbose=True):
  t0 = time.time()
  if (out is None) : out = sys.stdout
  from iotbx import file_reader
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="model",
    reflection_file_def="map_coeffs",
    map_file_def="map_file",
    cif_file_def="cif_file",
    usage_string="""\
mmtbx.ringer model.pdb map_coeffs.mtz [cif_file ...] [options]

%s
""" % __doc__)
  cmdline.work.show()
  params = cmdline.work.extract()
  validate_params(params)
  pdb_in = file_reader.any_file(params.model, force_type="pdb")
  pdb_in.check_file_type("pdb")

  pdb_inp = iotbx.pdb.input(file_name=params.model)
  model = mmtbx.model.manager(
    model_input      = pdb_inp)
  crystal_symmetry_model = model.crystal_symmetry()
  if crystal_symmetry_model is not None:
    crystal_symmetry_model.show_summary()

  hierarchy = model.get_hierarchy()
  map_coeffs = map_inp = difference_map_coeffs = None
  map_data, unit_cell = None, None
  # get miller array if map coefficients are provided
  if (params.map_coeffs is not None):
    mtz_in = file_reader.any_file(params.map_coeffs, force_type="hkl")
    mtz_in.check_file_type("hkl")
    best_guess = None
    best_labels = []
    all_labels = []
    for array in mtz_in.file_server.miller_arrays :
      if (array.is_complex_array()):
        labels = array.info().label_string()
        if (labels == params.map_label):
          map_coeffs = array
        elif (labels == params.difference_map_label):
          difference_map_coeffs = array
        else :
          if (params.map_label is None):
            all_labels.append(labels)
            if (labels.startswith("2FOFCWT") or labels.startswith("2mFoDFc") or
                labels.startswith("FWT")):
              best_guess = array
              best_labels.append(labels)
          if (params.difference_map_label is None):
            if (labels.startswith("FOFCWT") or labels.startswith("DELFWT")):
              difference_map_coeffs = array
    if (map_coeffs is None):
      if (len(all_labels) == 0):
        raise Sorry("No valid (pre-weighted) map coefficients found in file.")
      elif (best_guess is None):
        raise Sorry("Couldn't automatically determine appropriate map labels. "+
          "Choices:\n  %s" % "  \n".join(all_labels))
      elif (len(best_labels) > 1):
        raise Sorry("Multiple appropriate map coefficients found in file. "+
          "Choices:\n  %s" % "\n  ".join(best_labels))
      map_coeffs = best_guess
      print("  Guessing %s for input map coefficients" % best_labels[0], file=out)
  # get map_inp object and do sanity checks if map is provided
  else :
    ccp4_map_in = file_reader.any_file(params.map_file, force_type="ccp4_map")
    ccp4_map_in.check_file_type("ccp4_map")
    map_inp = ccp4_map_in.file_object
    base = map_model_manager(
      map_manager               = map_inp,
      model            = model,
      ignore_symmetry_conflicts = params.ignore_symmetry_conflicts)
    cs_consensus = base.crystal_symmetry()
    hierarchy = base.model().get_hierarchy()
    map_data = base.map_data()
    unit_cell = map_inp.grid_unit_cell()

  hierarchy.atoms().reset_i_seq()

  make_header("Iterating over residues", out=out)
  t1 = time.time()
  results = iterate_over_residues(
    pdb_hierarchy=hierarchy,
    map_coeffs=map_coeffs,
    difference_map_coeffs=difference_map_coeffs,
    map_data  = map_data,
    unit_cell = unit_cell,
    params=params,
    log=out).results
  t2 = time.time()
  if (verbose):
    print("Time excluding I/O: %8.1fs" % (t2 - t1), file=out)
    print("Overall runtime:    %8.1fs" % (t2 - t0), file=out)
  if (params.output_base is None):
    pdb_base = os.path.basename(params.model)
    params.output_base = os.path.splitext(pdb_base)[0] + "_ringer"
  easy_pickle.dump("%s.pkl" % params.output_base, results)
  print("Wrote %s.pkl" % params.output_base, file=out)
  csv = "\n".join([ r.format_csv() for r in results ])
  with open("%s.csv" % params.output_base, "w") as f:
    f.write(csv)
  print("Wrote %s.csv" % params.output_base, file=out)
  print("\nReference:", file=out)
  print("""\
  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387""", file=out)
  if (params.gui):
    run_app(results)
  else :
    return results

def validate_params(params):
  if (params.model is None):
    raise Sorry("No PDB file supplied (parameter: model)")
  if (params.map_coeffs is None) and (params.map_file is None):
    raise Sorry("No map coefficients supplied (parameter: map_coeffs)")
  if (not (1 < params.sampling_angle < 60)):
    raise Sorry("The sampling angle must be an integer between 1 and 60 "+
      "degrees")
  if (not (0 < params.grid_spacing < 0.5)):
    raise Sorry("The grid spacing must be greater than zero but less than 0.5")
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.makedirs(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args), out=sys.stdout)
########################################################################
# GUI
try :
  import wx
except ImportError :
  def run_app(results):
    raise Sorry("wxPython not available.")
else :
  from wxtbx import plots

  def run_app(results):
    app = wx.App(0)
    frame = RingerFrame(None, -1, "Ringer results")
    frame.show_results(results)
    frame.Show()
    app.MainLoop()

  class RingerFrame(plots.plot_frame):
    def create_plot_panel(self):
      plot = RingerPlot(self, figure_size=(6,8))
      plot.canvas.Bind(wx.EVT_CHAR, self.OnChar)
      return plot

    def draw_top_panel(self):
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

    def OnSelect(self, event):
      selection = event.GetEventObject().GetSelection()
      self.plot_panel.show_residue(self.results[selection])
      self.selection_callback(self.results[selection])

    # override in subclasses
    def selection_callback(self, residue):
      pass

    def show_results(self, results):
      self.results = results
      choices = [ result.format() for result in results ]
      self.chooser.SetItems(choices)
      self.chooser.SetSelection(0)
      self.plot_panel.show_residue(self.results[0])

    def OnChar(self, event):
      key = event.GetKeyCode()
      if (len(self.results) == 0) : return
      selection = self.chooser.GetSelection()
      if (key in [wx.WXK_TAB, wx.WXK_RETURN, wx.WXK_SPACE]):
        if (selection < (len(self.results) - 1)):
          selection += 1
        elif (len(self.results) > 0):
          selection = 0
      elif (key in [wx.WXK_DELETE, wx.WXK_BACK]):
        if (selection > 0):
          selection -= 1
        else :
          selection = len(results) - 1
      self.chooser.SetSelection(selection)
      self.plot_panel.show_residue(self.results[selection])

  class RingerPlot(plots.plot_container):
    def show_residue(self, residue, show_background_boxes=True):
      if (self.disabled) : return
      self.figure.clear()
      subplots = []
      for i in range(1, residue.n_chi + 1):
        chi = residue.get_angle(i)
        if (chi is None) : continue
        if (len(subplots) > 0):
          p = self.figure.add_subplot(4, 1, i, sharex=subplots[0])
        else :
          p = self.figure.add_subplot(4, 1, i)
          p.set_title(residue.format())
        p.set_position([0.15, 0.725 - 0.225*(i-1), 0.8, 0.225])
        x = [ k*chi.sampling for k in range(len(chi.densities)) ]
        p.plot(x, chi.densities, 'r-', linewidth=1)
        if (chi.fofc_densities is not None):
          p.plot(x, chi.fofc_densities, linestyle='--', color=[0.5,0.0,1.0])
        p.axvline(chi.angle_current, color='b', linewidth=2, linestyle='--')
        p.axhline(0, color=(0.4,0.4,0.4), linestyle='--', linewidth=1)
        if show_background_boxes:
          p.axhspan(0.3,1,facecolor="green",alpha=0.5)
          p.axhspan(-1,0.3,facecolor="grey",alpha=0.5)
        p.set_xlim(0,360)
        p.set_ylabel("Rho")
        p.set_xlabel("Chi%d" % i)
        subplots.append(p)
      for p in subplots[:-1] :
        for label in p.get_xticklabels():
          label.set_visible(False)
      self.canvas.draw()
      self.canvas.Fit()
      self.Layout()
      self.parent.Refresh()

if (__name__ == "__main__"):
  run(sys.argv[1:])
