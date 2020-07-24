from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import emringer
from six.moves import range

if __name__  ==  '__main__':
  run_program(program_class = emringer.Program)

#  =============================================================================
# old code - maybe necessary until GUI is updated

# LIBTBX_SET_DISPATCHER_NAME phenix.emringer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT = 1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

"""
Implementation of the EMRinger method, with plots, for use with the
EMRinger workflow.

References:
  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS. (2015)
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Nature Methods, in press.

  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387
"""

# Any software that wants to use the pkl output of this tool
# should import ringer_residue and ringer_chi from it.
#from __future__ import absolute_import, division, print_function
import libtbx.phil
from libtbx import easy_pickle
from libtbx.str_utils import make_header
from libtbx import runtime_utils
from libtbx.utils import Sorry
from iotbx.map_model_manager import map_model_manager
import time
import os
import sys

master_phil = libtbx.phil.parse("""
model = None
  .type = path
  .style = file_type:pdb bold input_file
  .short_caption = Model file
map_file = None
  .type = path
  .short_caption = CCP4 or MRC map
  .style = file_type:ccp4_map bold
map_coeffs = None
  .type = path
  .short_caption = Map coefficients
  .style = file_type:hkl input_file OnChange:extract_ringer_map_labels
map_label = 2FOFCWT, PH2FOFCWT
  .type = str
  .input_size = 200
  .short_caption = 2Fo-FC map labels
  .style = renderer:draw_map_arrays_widget noauto
sampling_angle = 5
  .type = int
  .input_size = 64
sampling_method = linear *spline direct
  .type = choice(multi = False)
grid_spacing = 1./5
  .type = float
scaling = *sigma volume
  .type = choice(multi = False)
rolling_window_threshold = 0
  .type = float(value_min = 0)
  .help = Threshold for calculating statistics across rolling windows of residues
skip_alt_confs = True
  .type = bool
nproc = 1
  .type = int
  .short_caption = Number of processors
  .input_size = 64
  .style = renderer:draw_nproc_widget
show_gui = False
  .type = bool
wrapping = None
  .type = bool
  .help = You can specify that the map wraps around the unit cell
ignore_symmetry_conflicts = False
  .type = bool
  .help = You can specify that the model and map symmetryies to not have to \
      match
output_base = None
  .type = str
output_dir = None
  .type = path
  .short_caption = Output directory
  .style = output_dir
include scope libtbx.phil.interface.tracking_params
""", process_includes = True)
master_params = master_phil # XXX Gui hack

def run(args, out = None, verbose = True, plots_dir = None):
  t0 = time.time()
  if (out is None) : out = sys.stdout
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args = args,
    master_phil = master_phil,
    pdb_file_def = "model",
    reflection_file_def = "map_coeffs",
    map_file_def = "map_file",
    usage_string = """\
phenix.emringer model.pdb map.mrc [cif_file ...] [options]

%s
""" % __doc__)
  params = cmdline.work.extract()
  validate_params(params)
  from iotbx.data_manager import DataManager
  dm = DataManager()
  model = dm.get_model(params.model)
  crystal_symmetry_model = model.crystal_symmetry()
  hierarchy = model.get_hierarchy()
  map_coeffs = map_inp = None
  map_data, unit_cell = None, None
  if (params.map_coeffs is not None):
    mtz_in = cmdline.get_file(params.map_coeffs)
    mtz_in.check_file_type("hkl")
    best_guess = None
    best_labels = []
    all_labels = []
    for array in mtz_in.file_server.miller_arrays :
      if (array.info().label_string()  ==  params.map_label):
        map_coeffs = array
        break
      elif (params.map_label is None):
        if (array.is_complex_array()):
          labels = array.info().label_string()
          all_labels.append(labels)
          if (labels.startswith("2FOFCWT") or labels.startswith("2mFoDFc") or
              labels.startswith("FWT")):
            best_guess = array
            best_labels.append(labels)
    if (map_coeffs is None):
      if (len(all_labels)  ==  0):
        raise Sorry("No valid (pre-weighted) map coefficients found in file.")
      elif (best_guess is None):
        raise Sorry("Couldn't automatically determine appropriate map labels. "+
          "Choices:\n  %s" % "  \n".join(all_labels))
      elif (len(best_labels) > 1):
        raise Sorry("Multiple appropriate map coefficients found in file. "+
          "Choices:\n  %s" % "\n  ".join(best_labels))
      map_coeffs = best_guess
      print("  Guessing %s for input map coefficients"% best_labels[0], file = out)
  else :
    ccp4_map_in = cmdline.get_file(params.map_file)
    ccp4_map_in.check_file_type("ccp4_map")
    map_inp = ccp4_map_in.file_object
    base = map_model_manager(
      map_manager               = map_inp,
      model            = model,
      wrapping = params.wrapping,
      ignore_symmetry_conflicts = params.ignore_symmetry_conflicts)
    hierarchy = base.model().get_hierarchy()
    map_data = base.map_data()
    unit_cell = map_inp.grid_unit_cell()

  hierarchy.atoms().reset_i_seq()
  make_header("Iterating over residues", out = out)
  t1 = time.time()
  from mmtbx.ringer import iterate_over_residues
  results = iterate_over_residues(
    pdb_hierarchy = hierarchy,
    map_coeffs = map_coeffs,
    map_data  = map_data,
    unit_cell = unit_cell,
    params = params,
    log = out).results
  t2 = time.time()
  if (verbose):
    print("Time excluding I/O: %8.1fs" % (t2 - t1), file = out)
    print("Overall runtime:    %8.1fs" % (t2 - t0), file = out)
  if (params.output_base is None):
    pdb_base = os.path.basename(params.model)
    params.output_base = os.path.splitext(pdb_base)[0] + "_emringer"
  easy_pickle.dump("%s.pkl" % params.output_base, results)
  print("Wrote %s.pkl" % params.output_base, file = out)
  csv = "\n".join([ r.format_csv() for r in results ])
  open("%s.csv" % params.output_base, "w").write(csv)
  print("Wrote %s.csv" % params.output_base, file = out)
  if (plots_dir is None):
    plots_dir = params.output_base + "_plots"
  if (not os.path.isdir(plots_dir)):
    os.makedirs(plots_dir)
  from mmtbx.ringer import em_rolling
  from mmtbx.ringer import em_scoring
  import matplotlib
  matplotlib.use("Agg")
  make_header("Scoring results", out = out)
  scoring = em_scoring.main(
    file_name = params.output_base,
    ringer_result = results,
    out_dir = plots_dir,
    sampling_angle = params.sampling_angle,
    quiet = False,
    out = out)
  make_header("Inspecting chains", out = out)
  rolling_window_threshold = params.rolling_window_threshold
  rolling = em_rolling.main(
    ringer_results = results,
    dir_name = plots_dir,
    threshold = rolling_window_threshold, #scoring.optimal_threshold,
    graph = False,
    save = True,
    out = out)
  scoring.show_summary(out = out)
  print("\nReferences:", file = out)

  references = """\
  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS. (2015)
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Nature Methods, in press.

  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387"""
  print(references, file = out)
  if (params.show_gui):
    run_app(results)
  else :
    return (results, scoring, rolling)

def validate_params(params):
  if (params.model is None):
    raise Sorry("No PDB file supplied (parameter: model)")
  if (params.map_coeffs is None) and (params.map_file is None):
    raise Sorry("No map coefficients supplied (parameter: map_coeffs)")
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.makedirs(self.output_dir)
    os.chdir(self.output_dir)
    return run(args = list(self.args), out = sys.stdout, plots_dir = "plots")

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
      plot = RingerPlot(self, figure_size = (6, 8))
      plot.canvas.Bind(wx.EVT_CHAR, self.OnChar)
      return plot

    def draw_top_panel(self):
      self.top_panel = wx.Panel(self, style = wx.SUNKEN_BORDER)
      panel_szr = wx.BoxSizer(wx.VERTICAL)
      self.top_panel.SetSizer(panel_szr)
      szr2 = wx.BoxSizer(wx.HORIZONTAL)
      panel_szr.Add(szr2)
      txt1 = wx.StaticText(self.top_panel, -1, "Residue to display:")
      szr2.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.chooser = wx.Choice(self.top_panel, -1, size = (200, -1))
      szr2.Add(self.chooser, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.Bind(wx.EVT_CHOICE, self.OnSelect, self.chooser)
      self.Bind(wx.EVT_CHAR, self.OnChar)
      self.chooser.Bind(wx.EVT_CHAR, self.OnChar)
      return self.top_panel

    def OnSelect(self, event):
      selection = event.GetEventObject().GetSelection()
      self.plot_panel.show_residue(self.results[selection])

    def show_results(self, results):
      self.results = results
      choices = [ result.format() for result in results ]
      self.chooser.SetItems(choices)
      self.chooser.SetSelection(0)
      self.plot_panel.show_residue(self.results[0])

    def OnChar(self, event):
      key = event.GetKeyCode()
      if (len(self.results)  ==  0) : return
      selection = self.chooser.GetSelection()
      if (key in [wx.WXK_TAB, wx.WXK_RETURN, wx.WXK_SPACE]):
        if (selection < (len(self.results) - 1)):
          selection +=  1
        elif (len(self.results) > 0):
          selection = 0
      elif (key in [wx.WXK_DELETE, wx.WXK_BACK]):
        if (selection > 0):
          selection -=  1
        else :
          selection = len(results) - 1
      self.chooser.SetSelection(selection)
      self.plot_panel.show_residue(self.results[selection])

  class RingerPlot(plots.plot_container):
    def show_residue(self, residue, show_background_boxes = False):
      if (self.disabled) : return
      self.figure.clear()
      subplots = []
      for i in range(1, residue.n_chi + 1):
        chi = residue.get_angle(i)
        if (chi is None) : continue
        if (len(subplots) > 0):
          p = self.figure.add_subplot(4, 1, i, sharex = subplots[0])
        else :
          p = self.figure.add_subplot(4, 1, i)
          p.set_title(residue.format())
        p.set_position([0.15, 0.725 - 0.225*(i-1), 0.8, 0.225])
        x = [ k*chi.sampling for k in range(len(chi.densities)) ]
        p.plot(x, chi.densities, 'r-', linewidth = 1)
        p.axvline(chi.angle_current, color = 'b', linewidth = 2, linestyle = '--')
        p.axvline(chi.peak_chi, color = 'g', linewidth = 2, linestyle = '--')
        p.axhline(0, color = (0.4, 0.4, 0.4), linestyle = '--', linewidth = 1)
        if show_background_boxes:
          p.axhspan(0.3, 1, facecolor = "green", alpha = 0.5)
          p.axhspan(-1, 0.3, facecolor = "grey", alpha = 0.5)
        p.set_xlim(0, 360)
        p.set_ylabel("Rho")
        p.set_xlabel("Chi%d" % i)
        subplots.append(p)
      for p in subplots[:-1] :
        for label in p.get_xticklabels():
          label.set_visible(False)
      p.text(0, -0.5, 'Green = Peak, Blue = Modelled',
          transform = ax.transAxes)
      self.canvas.draw()
      self.canvas.Fit()
      self.Layout()
      self.parent.Refresh()

#if (__name__  ==  "__main__"):
#  run(sys.argv[1:])
