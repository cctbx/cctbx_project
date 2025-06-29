"""Find peaks and holes in a difference map"""
# LIBTBX_SET_DISPATCHER_NAME phenix.find_peaks_holes
# LIBTBX_SET_DISPATCHER_NAME mmtbx.find_peaks_holes

# simple frontend to mmtbx.find_peaks, primarily intended for use in quickly
# analyzing structures in the PDB (and storing results)

from __future__ import absolute_import, division, print_function
from mmtbx import utils
from scitbx.array_family import flex
from libtbx.str_utils import make_header, format_value
from libtbx import runtime_utils
from libtbx.utils import Sorry
import libtbx.phil
from libtbx import adopt_init_args, group_args
from iotbx.pdb.hybrid_36 import hy36encode
import operator
import os
import sys
from six.moves import zip

def get_master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    phil_string="""
find_peaks
  .style = auto_align
{
  include scope mmtbx.find_peaks.master_params
}
map_cutoff = 3.0
  .type = float
  .short_caption = mFo-DFc map cutoff (sigma)
anom_map_cutoff = 3.0
  .type = float
  .short_caption = Anomalous map cutoff (sigma)
include_peaks_near_model = False
  .type = bool
  .short_caption = Don't filter peaks by distance to model
  .help = By default, the program will only display peaks that map to points \
    outside of the current model, ignoring those that overlap with atoms.  \
    Setting this option to True is equivalent to specifying a distance \
    cutoff of zero for the filtering step.
  .style = OnChange:toggle_min_model_peak_dist
wavelength = None
  .type = float
  .help = Optional parameter, if defined this will cause all atoms to be \
    treated as anomalous scatterers using the standard Sasaki table to \
    obtain theoretical fp and fpp values.  Only really useful if the Phaser \
    LLG map is being used for the anomalous map.
filter_peaks_by_2fofc = None
  .type = float
  .short_caption = Filter peaks by 2mFo-DFc
  .help = If this is set, peaks outside 2mFo-DFc density at the \
    cutoff will be discarded.  (This does not apply to the analysis of \
    solvent atoms.)  Holes will not be changed.
use_phaser_if_available = True
  .type = bool
  .short_caption = Use Phaser LLG map
  .help = If True, and Phaser is installed and configured, an anomalous LLG \
    map will be used in place of the simple anomalous difference map.  The \
    wavelength should be specified for this to be maximally useful.
write_pdb = True
  .type = bool
  .short_caption = Write peaks to PDB file
write_maps = True
  .type = bool
  .short_caption = Save map coefficients
output_file_prefix = peaks_holes
  .type = str
include scope libtbx.phil.interface.tracking_params
""")

master_phil = get_master_phil()
master_params = master_phil # for phenix GUI

class peaks_holes_container(object):
  def __init__(self, peaks, holes, map_cutoff=3.0, anom_peaks=None,
      anom_map_cutoff=3.0, water_peaks=None, water_anom_peaks=None,
      non_water_anom_peaks=None):
    adopt_init_args(self, locals())
    # XXX pre-sort all lists
    self.peaks.sort(reverse=True)
    self.holes.sort()
    if (self.anom_peaks is not None):
      self.anom_peaks.sort(reverse=True)
    if (self.water_peaks is not None):
      self.water_peaks = sorted(self.water_peaks, key=operator.attrgetter("peak_height"), reverse=True)
    if (self.water_anom_peaks is not None):
      self.water_anom_peaks = sorted(self.water_anom_peaks,
                                     key=operator.attrgetter("peak_height"),
                                     reverse=True)
    self.pdb_file = None
    self.map_file = None

  def show_summary(self, out=sys.stdout):
    self.get_summary().show(out=out)

  def get_summary(self):
    """
    Returns a simple object for harvesting statistics elsewhere.
    """
    n_anom_peaks = None
    if (self.anom_peaks is not None):
      n_anom_peaks = len(self.anom_peaks.heights)
    n_water_peaks = n_water_anom_peaks = None
    if (self.water_peaks is not None):
      n_water_peaks = len(self.water_peaks)
    if (self.water_anom_peaks is not None):
      n_water_anom_peaks = len(self.water_anom_peaks)
    hole_max = peak_max = None
    if (len(self.peaks.heights) > 0):
      peak_max = flex.max(self.peaks.heights)
    if (len(self.holes.heights) > 0):
      hole_max = flex.min(self.holes.heights)
    n_non_water_anom_peaks = None
    if (getattr(self, "non_water_anom_peaks", None) is not None):
      n_non_water_anom_peaks = len(self.non_water_anom_peaks)
    return summary(
      n_peaks_1=(self.peaks.heights > self.map_cutoff).count(True),
      n_peaks_2=(self.peaks.heights > self.map_cutoff + 3).count(True),
      n_peaks_3=(self.peaks.heights > self.map_cutoff + 6).count(True),
      n_holes_1=(self.holes.heights < -self.map_cutoff).count(True),
      n_holes_2=(self.holes.heights < -self.map_cutoff - 3).count(True),
      n_holes_3=(self.holes.heights < -self.map_cutoff - 6).count(True),
      peak_max=peak_max,
      hole_max=hole_max,
      n_anom_peaks=n_anom_peaks,
      n_water_peaks=n_water_peaks,
      n_water_anom_peaks=n_water_anom_peaks,
      map_cutoff=self.map_cutoff,
      anom_map_cutoff=self.anom_map_cutoff,
      n_non_water_anom_peaks=n_non_water_anom_peaks)

  def n_peaks_above_cutoff(self, cutoff):
    assert (cutoff > 0)
    return (self.peaks.heights > cutoff).count(True)

  def n_holes_below_cutoff(self, cutoff):
    assert (cutoff < 0)
    return (self.holes.heights < cutoff).count(True)

  def save_pdb_file(self,
      file_name="peaks.pdb",
      include_holes=True,
      include_anom=True,
      include_water=True,
      log=None):
    """
    Write out a PDB file with up to three chains: A for peaks, B for holes,
    C for anomalous peaks.  Atoms are UNK, with the B-factor set to the height
    or depth of the peak or hole.
    """
    if (log is None) : log = sys.stdout
    import iotbx.pdb.hierarchy
    self.peaks.sort(reverse=True)
    root = iotbx.pdb.hierarchy.root()
    model = iotbx.pdb.hierarchy.model()
    root.append_model(model)
    peaks_chain = iotbx.pdb.hierarchy.chain(id="A")
    model.append_chain(peaks_chain)
    def create_atom(xyz, peak, serial):
      rg = iotbx.pdb.hierarchy.residue_group(resseq=hy36encode(4, serial))
      ag = iotbx.pdb.hierarchy.atom_group(resname="UNK")
      rg.append_atom_group(ag)
      a = iotbx.pdb.hierarchy.atom()
      ag.append_atom(a)
      a.name = " UNK"
      a.element = "X"
      a.xyz = xyz
      a.b = peak
      a.occ = 1.
      a.serial = serial
      return rg
    k = 1
    for peak, xyz in zip(self.peaks.heights, self.peaks.sites):
      rg = create_atom(xyz, peak, k)
      peaks_chain.append_residue_group(rg)
      k += 1
    f = open(file_name, "w")
    f.write("REMARK  Interesting sites from mmtbx.find_peaks_holes\n")
    f.write("REMARK  Chain A is difference map peaks (> %g sigma)\n" %
      self.map_cutoff)
    if (include_holes):
      f.write("REMARK  Chain B is difference map holes (< %g sigma)\n" %
        (- self.map_cutoff))
      holes_chain = iotbx.pdb.hierarchy.chain(id="B")
      model.append_chain(holes_chain)
      k = 1
      for hole, xyz in zip(self.holes.heights, self.holes.sites):
        rg = create_atom(xyz, hole, k)
        holes_chain.append_residue_group(rg)
        k += 1
    if (include_anom) and (self.anom_peaks is not None):
      f.write("REMARK  Chain C is anomalous peaks (> %g sigma)\n" %
        self.anom_map_cutoff)
      anom_chain = iotbx.pdb.hierarchy.chain(id="C")
      model.append_chain(anom_chain)
      k = 1
      for peak, xyz in zip(self.anom_peaks.heights, self.anom_peaks.sites):
        rg = create_atom(xyz, peak, k)
        anom_chain.append_residue_group(rg)
        k += 1
    if (include_water) and (self.water_peaks is not None):
      f.write("REMARK  Chain D is waters with mFo-DFc peaks (> %g sigma)\n" %
        self.map_cutoff)
      waters_chain = iotbx.pdb.hierarchy.chain(id="D")
      model.append_chain(waters_chain)
      for k, peak in enumerate(self.water_peaks):
        rg = create_atom(peak.xyz, peak.peak_height, k+1)
        waters_chain.append_residue_group(rg)
      if (include_anom) and (self.water_anom_peaks is not None):
        f.write("REMARK  Chain E is waters with anom. peaks (> %g sigma)\n" %
          self.anom_map_cutoff)
        waters_chain_2 = iotbx.pdb.hierarchy.chain(id="E")
        model.append_chain(waters_chain_2)
        for k, peak in enumerate(self.water_anom_peaks):
          rg = create_atom(peak.xyz, peak.peak_height, k+1)
          waters_chain_2.append_residue_group(rg)
    if ((include_anom) and
        (getattr(self, "non_water_anom_peaks", None) is not None)):
      f.write("REMARK  Chain F is non-water, non-HD atoms with anom. peaks\n")
      anom_chain_2 = iotbx.pdb.hierarchy.chain(id="F")
      model.append_chain(anom_chain_2)
      for k, peak in enumerate(self.non_water_anom_peaks):
        rg = create_atom(peak.xyz, peak.peak_height, k+1)
        anom_chain_2.append_residue_group(rg)
    f.write(root.as_pdb_string())
    f.close()
    print("Wrote %s" % file_name, file=log)
    self.pdb_file = file_name

  def get_output_file_info(self):
    output_files = []
    if (self.pdb_file is not None):
      output_files.append((self.pdb_file, "Peaks as PDB atoms"))
    if (self.map_file is not None):
      output_files.append((self.map_file, "Map coefficients"))
    return output_files

class summary(group_args):
  def show(self, out=sys.stdout):
    print("", file=out)
    print("SUMMARY OF MAP PEAKS:", file=out)
    cutoffs = [self.map_cutoff, self.map_cutoff + 3.0, self.map_cutoff + 6.0]
    peaks = [ self.n_peaks_1, self.n_peaks_2, self.n_peaks_3 ]
    labels = []
    values = []
    for cutoff, n_peaks in zip(cutoffs, peaks):
      labels.append("mFo-DFc >  %-4g" % cutoff)
      values.append("%6d" % n_peaks)
    labels.append("mFo-DFc max")
    values.append(format_value("%6.2f", self.peak_max))
    holes = [ self.n_holes_1, self.n_holes_2, self.n_holes_3 ]
    for cutoff, n_holes in zip(cutoffs, holes):
      labels.append("mFo-DFc < -%-4g" % cutoff)
      values.append("%6d" % n_holes)
    labels.append("mFo-DFc min")
    values.append(format_value("%6.2f", self.hole_max))
    if (self.n_anom_peaks is not None):
      labels.append("anomalous > %-4g" % self.anom_map_cutoff)
      values.append("%6d" % self.n_anom_peaks)
    if (self.n_water_peaks is not None):
      labels.append("suspicious H2O (mFo-DFC > %g)" % self.map_cutoff)
      values.append("%6d" % self.n_water_peaks)
    if (self.n_water_anom_peaks is not None):
      labels.append("anomalous H2O (anomalous > %g)" % self.map_cutoff)
      values.append("%6d" % self.n_water_anom_peaks)
    if (self.n_non_water_anom_peaks is not None):
      labels.append("anomalous non-water atoms")
      values.append("%6d" % self.n_non_water_anom_peaks)
    labels = [ l.strip() + ":" for l in labels ]
    label_len = max([ len(l) for l in labels ])
    format = "%%-%ds" % label_len
    for label, value in zip(labels, values):
      formatted = format % label
      print("  %s %s" % (formatted, value), file=out)
    print("", file=out)

class water_peak(object):
  def __init__(self, id_str, xyz, peak_height, map_type="mFo-DFc"):
    adopt_init_args(self, locals())

  def show(self, out=sys.stdout):
    print("  %s  map_type=%s  peak=%g" % (self.id_str,
      self.map_type, self.peak_height), file=out)

def find_peaks_holes(
    fmodel,
    pdb_hierarchy,
    params=None,
    map_cutoff=3.0,
    anom_map_cutoff=3.0,
    filter_peaks_by_2fofc=None,
    use_phaser_if_available=True,
    return_llg_map=False,
    include_peaks_near_model=False,
    out=None):
  """
  Find peaks and holes in mFo-DFc map, plus flag solvent atoms with
  suspiciously high mFo-DFc values, plus anomalous peaks if anomalous data are
  present.  Returns a pickle-able object storing all this information (with
  the ability to write out a PDB file with the sites of interest).
  """
  if (out is None) : out = sys.stdout
  if (params is None):
    params = master_phil.fetch().extract().find_peaks
  if (include_peaks_near_model):
    params.map_next_to_model.min_model_peak_dist = 0
  pdb_atoms = pdb_hierarchy.atoms()
  unit_cell = fmodel.xray_structure.unit_cell()
  from mmtbx import find_peaks
  from cctbx import maptbx
  f_map = None
  if (filter_peaks_by_2fofc is not None):
    f_map_ = fmodel.electron_density_map().fft_map(
      resolution_factor=min(0.5, params.grid_step/fmodel.f_obs().d_min()),
      symmetry_flags=maptbx.use_space_group_symmetry,
      map_type="2mFo-DFc",
      use_all_data=True)
    f_map_.apply_sigma_scaling()
    f_map = f_map_.real_map()
  make_header("Positive difference map peaks", out=out)
  coeffs = fmodel.electron_density_map().map_coefficients(
    map_type     = "mFo-DFc",
    fill_missing = False,
    isotropize   = False)
  peaks_result = find_peaks.manager(
    map_coeffs = coeffs,
    xray_structure = fmodel.xray_structure,
    map_cutoff=map_cutoff,
    params=params,
    log=out)
  peaks_result.peaks_mapped()
  peaks_result.show_mapped(pdb_atoms)
  peaks = peaks_result.peaks()
  if (filter_peaks_by_2fofc is not None):
    n_removed = peaks.filter_by_secondary_map(
      map=f_map,
      min_value=filter_peaks_by_2fofc)
    print("", file=out)
    print("%d peaks remaining after 2mFo-DFc filtering" % \
      len(peaks.sites), file=out)
  # very important - sites are initially fractional coordinates!
  peaks.sites = unit_cell.orthogonalize(peaks.sites)
  print("", file=out)
  out.flush()
  make_header("Negative difference map holes", out=out)
  coeffs = fmodel.electron_density_map().map_coefficients(
    map_type     = "mFo-DFc",
    fill_missing = False,
    isotropize   = False)
  holes_result = find_peaks.manager(
    map_coeffs = coeffs,
    xray_structure = fmodel.xray_structure,
    map_cutoff=-map_cutoff,
    params=params,
    log=out)
  holes_result.peaks_mapped()
  holes_result.show_mapped(pdb_atoms)
  holes = holes_result.peaks()
  # XXX is this useful?
  #if (filter_peaks_by_2fofc is not None):
  #  holes.filter_by_secondary_map(
  #    map=f_map,
  #    min_value=filter_peaks_by_2fofc)
  holes.sites = unit_cell.orthogonalize(holes.sites)
  print("", file=out)
  out.flush()
  anom = None
  anom_map_coeffs = None
  if (fmodel.f_obs().anomalous_flag()):
    make_header("Anomalous difference map peaks", out=out)
    anom_map_type = "anom_residual"
    if ((use_phaser_if_available) and (libtbx.env.has_module("phaser")) and
        (not fmodel.twin)):
      import mmtbx.map_tools
      print("Will use Phaser LLG map", file=out)
      anom_map_type = None
      anom_map_coeffs = mmtbx.map_tools.get_phaser_sad_llg_map_coefficients(
        fmodel=fmodel,
        pdb_hierarchy=pdb_hierarchy,
        log=out)
    anom_result = find_peaks.manager(
      fmodel=fmodel,
      map_type=anom_map_type,
      map_coeffs=anom_map_coeffs,
      map_cutoff=anom_map_cutoff,
      params=params,
      log=out)
    anom_result.peaks_mapped()
    anom_result.show_mapped(pdb_atoms)
    anom = anom_result.peaks()
    if (filter_peaks_by_2fofc is not None):
      anom.filter_by_secondary_map(
        map=f_map,
        min_value=filter_peaks_by_2fofc)
      print("", file=out)
      print("%d peaks remaining after 2mFo-DFc filtering" % \
        len(anom.sites), file=out)
    anom.sites = unit_cell.orthogonalize(anom.sites)
    print("", file=out)
    out.flush()
  anom_map = None
  cache = pdb_hierarchy.atom_selection_cache()
  sites_frac = fmodel.xray_structure.sites_frac()
  water_isel = cache.selection(
    "resname HOH and not (element H or element D)").iselection()
  waters_out = [None, None]
  if (len(water_isel) > 0):
    map_types = ["mFo-DFc"]
    map_cutoffs = [ map_cutoff ]
    if (fmodel.f_obs().anomalous_flag()):
      map_types.append("anomalous")
      map_cutoffs.append(anom_map_cutoff)
    for k, map_type in enumerate(map_types):
      fft_map = None
      # re-use Phaser LLG map if it was previously calculated
      if (map_type == "anomalous") and (anom_map_coeffs is not None):
        fft_map = anom_map_coeffs.fft_map(
          resolution_factor=min(0.5, params.grid_step/fmodel.f_obs().d_min()),
          symmetry_flags=maptbx.use_space_group_symmetry)
      else :
        fft_map = fmodel.electron_density_map().fft_map(
          resolution_factor= min(0.5, params.grid_step/fmodel.f_obs().d_min()),
          symmetry_flags=maptbx.use_space_group_symmetry,
          map_type=map_type,
          use_all_data=True)
      real_map = fft_map.apply_sigma_scaling().real_map_unpadded()
      if (map_type == "anomalous") : anom_map = real_map
      suspicious_waters = []
      for i_seq in water_isel :
        atom = pdb_atoms[i_seq]
        rho = real_map.tricubic_interpolation(sites_frac[i_seq])
        if (rho >= map_cutoffs[k]):
          peak = water_peak(
            id_str=atom.id_str(),
            xyz=atom.xyz,
            peak_height=rho,
            map_type=map_type)
          suspicious_waters.append(peak)
      if (len(suspicious_waters) > 0):
        make_header("Water molecules with %s peaks" % map_type, out=out)
        for peak in suspicious_waters :
          peak.show(out=out)
        print("", file=out)
        waters_out[k] = suspicious_waters
  non_water_anom_peaks = None
  if (fmodel.f_obs().anomalous_flag()):
    non_water_anom_peaks = []
    if (anom_map is None):
      fft_map = fmodel.electron_density_map().fft_map(
        resolution_factor=min(0.5, params.grid_step/fmodel.f_obs().d_min()),
        symmetry_flags=maptbx.use_space_group_symmetry,
        map_type="anom",
        use_all_data=True)
      anom_map = fft_map.apply_sigma_scaling().real_map_unpadded()
    non_water_non_H_i_sel = cache.selection(
      "not (resname HOH or element H or element D)").iselection()
    for i_seq in non_water_non_H_i_sel :
      rho = anom_map.tricubic_interpolation(sites_frac[i_seq])
      if (rho >= anom_map_cutoff):
        atom = pdb_atoms[i_seq]
        peak = water_peak(
          id_str=atom.id_str(),
          xyz=atom.xyz,
          peak_height=rho,
          map_type="anomalous")
        non_water_anom_peaks.append(peak)
  all_results = peaks_holes_container(
    peaks=peaks,
    holes=holes,
    anom_peaks=anom,
    map_cutoff=map_cutoff,
    anom_map_cutoff=anom_map_cutoff,
    water_peaks=waters_out[0],
    water_anom_peaks=waters_out[1],
    non_water_anom_peaks=non_water_anom_peaks)
  all_results.show_summary(out=out)
  if (return_llg_map):
    return all_results, anom_map_coeffs
  return all_results

def run(args, out=None):
  if (out is None) : out = sys.stdout
  import mmtbx.command_line
  usage_string = """\
mmtbx.find_peaks_holes - difference map analysis
  Prints a summary of all peaks and holes above the specified cutoff in the
  mFo-DFc map, and flag any water molecules with suspiciously high peaks
  (possible ions).  Will also check the anomalous map if available.
"""
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    create_fmodel=True,
    prefer_anomalous=True,
    usage_string=usage_string)
  fmodel = cmdline.fmodel
  params = cmdline.params
  if (params.wavelength is not None):
    xrs = fmodel.xray_structure
    xrs.set_inelastic_form_factors(
      photon=params.wavelength,
      table="sasaki")
    fmodel.update_xray_structure(xrs, update_f_calc=True)
  out.flush()
  result, llg_map = find_peaks_holes(
    fmodel=fmodel,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    params=params.find_peaks,
    map_cutoff=params.map_cutoff,
    anom_map_cutoff=params.anom_map_cutoff,
    filter_peaks_by_2fofc=params.filter_peaks_by_2fofc,
    use_phaser_if_available=True,
    return_llg_map=True,
    include_peaks_near_model=params.include_peaks_near_model,
    out=out)
  prefix = cmdline.params.output_file_prefix
  if (cmdline.params.write_pdb):
    result.save_pdb_file(file_name="%s.pdb" % prefix, log=out)
  if (cmdline.params.write_maps):
    import mmtbx.maps.utils
    import iotbx.map_tools
    f_map, diff_map = mmtbx.maps.utils.get_maps_from_fmodel(fmodel)
    anom_map = llg_map # use LLG map for anomalous map if available
    if (fmodel.f_obs().anomalous_flag()) and (anom_map is None):
      anom_map = mmtbx.maps.utils.get_anomalous_map(fmodel)
    iotbx.map_tools.write_map_coeffs(
      fwt_coeffs=f_map,
      delfwt_coeffs=diff_map,
      file_name="%s_maps.mtz" % prefix,
      anom_coeffs=anom_map)
    result.map_file = "%s_maps.mtz" % prefix
  return result

def validate_params(params, callback=None):
  from mmtbx.command_line import validate_input_params
  validate_input_params(params)
  if (params.find_peaks.map_next_to_model.min_model_peak_dist < 0):
    raise Sorry("The parameter 'Minimum distance from model' must be at least"+
      " zero.")
  if (params.find_peaks.peak_search.min_cross_distance <= 0):
    raise Sorry("The parameter 'Minimum cross distance' must be greater than "+
      "zero.")

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.chdir(self.output_dir)
    return run(args=list(self.args), out=sys.stdout)

def finish_job(result):
  output_files = []
  stats = []
  if (result is not None):
    output_files = result.get_output_file_info()
  return (output_files, stats)

if (__name__ == "__main__"):
  run(sys.argv[1:])

