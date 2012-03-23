# LIBTBX_SET_DISPATCHER_NAME phenix.find_peaks_holes

# simple frontend to mmtbx.find_peaks, primarily intended for use in quickly
# analyzing structures in the PDB (and storing results)

from mmtbx import utils
from scitbx.array_family import flex
from libtbx.str_utils import make_header, format_value
from libtbx import runtime_utils
from libtbx.utils import Usage, Sorry
import libtbx.phil
from libtbx import adopt_init_args, group_args
from cStringIO import StringIO
import os
import sys

master_phil = libtbx.phil.parse("""
%s
find_peaks {
  include scope mmtbx.find_peaks.master_params
}
map_cutoff = 3.0
  .type = float
  .short_caption = mFo-DFc map cutoff (sigma)
anom_map_cutoff = 3.0
  .type = float
  .short_caption = Anomalous map cutoff (sigma)
filter_peaks_by_2fofc = None
  .type = float
  .short_caption = Filter peaks by 2mFo-DFc
  .help = If this is set, peaks outside 2mFo-DFc density at the \
    cutoff will be discarded.  (This does not apply to the analysis of \
    solvent atoms.)  Holes will not be changed.
write_pdb = True
  .type = bool
  .short_caption = Write peaks to PDB file
write_maps = True
  .type = bool
  .short_caption = Save map coefficients
output_file_prefix = peaks_holes
  .type = str
include scope libtbx.phil.interface.tracking_params
""" % utils.cmdline_input_phil_str,
  process_includes=True)

master_params = master_phil # for phenix GUI

class peaks_holes_container (object) :
  def __init__ (self, peaks, holes, map_cutoff=3.0, anom_peaks=None,
      anom_map_cutoff=3.0, water_peaks=None, water_anom_peaks=None) :
    adopt_init_args(self, locals())
    # XXX pre-sort all lists
    self.peaks.sort(reverse=True)
    self.holes.sort()
    if (self.anom_peaks is not None) :
      self.anom_peaks.sort(reverse=True)
    if (self.water_peaks is not None) :
      self.water_peaks = sorted(self.water_peaks,
        lambda x,y: cmp(y.peak_height, x.peak_height))
    if (self.water_anom_peaks is not None) :
      self.water_anom_peaks = sorted(self.water_anom_peaks,
        lambda x,y: cmp(y.peak_height, x.peak_height))
    self.pdb_file = None
    self.map_file = None

  def show_summary (self, out=sys.stdout) :
    print >> out, ""
    print >> out, "SUMMARY OF MAP PEAKS:"
    cutoffs = [self.map_cutoff, self.map_cutoff + 3.0, self.map_cutoff + 6.0]
    for cutoff in cutoffs :
      n_peaks = (self.peaks.heights > cutoff).count(True)
      print >> out, "  mFo-DFc >  %-4g   : %6d" % (cutoff, n_peaks)
    if (len(self.peaks.heights) > 0) :
      peak_max = flex.max(self.peaks.heights)
    else :
      peak_max = None
    print >> out, "  mFo-DFc max       : %s" % format_value("%6.2f", peak_max)
    for cutoff in cutoffs :
      n_holes = (self.holes.heights < -cutoff).count(True)
      print >> out, "  mFo-DFc < -%-4g   : %6d" % (cutoff, n_holes)
    if (len(self.holes.heights) > 0) :
      hole_max = flex.min(self.holes.heights)
    else :
      hole_max = None
    print >> out, "  mFo-DFc min       : %s" % format_value("%6.2f", hole_max)
    if (self.anom_peaks is not None) :
      print >> out, "  anomalous > %-4g : %6d" % (self.anom_map_cutoff,
        len(self.anom_peaks.heights))
    if (self.water_peaks is not None) :
      print >> out, "  suspicious H2O (mFo-DFC > %g) : %6d" % (self.map_cutoff,
        len(self.water_peaks))
    if (self.water_anom_peaks is not None) :
      print >> out, "  anomalous H2O (anomalous > %g): %6d" % (self.map_cutoff,
        len(self.water_anom_peaks))
    print >> out, ""

  def get_summary (self) :
    """
    Returns a simple object for harvesting statistics elsewhere.
    """
    n_anom_peaks = None
    if (self.anom_peaks is not None) :
      n_anom_peaks = len(self.anom_peaks.heights)
    n_water_peaks = n_water_anom_peaks = None
    if (self.water_peaks is not None) :
      n_water_peaks = len(self.water_peaks)
    if (self.water_anom_peaks is not None) :
      n_water_anom_peaks = len(self.water_anom_peaks)
    return group_args(
      n_peaks_1=(self.peaks.heights > self.map_cutoff).count(True),
      n_peaks_2=(self.peaks.heights > self.map_cutoff + 3).count(True),
      n_peaks_3=(self.peaks.heights > self.map_cutoff + 6).count(True),
      n_holes_1=(self.holes.heights < -self.map_cutoff).count(True),
      n_holes_2=(self.holes.heights < -self.map_cutoff - 3).count(True),
      n_holes_3=(self.holes.heights < -self.map_cutoff - 6).count(True),
      peak_max=flex.max(self.peaks.heights),
      hole_max=flex.min(self.holes.heights),
      n_anom_peaks=n_anom_peaks,
      n_water_peaks=n_water_peaks,
      n_water_anom_peaks=n_water_anom_peaks)

  def n_peaks_above_cutoff (self, cutoff) :
    assert (cutoff > 0)
    return (self.peaks.heights > cutoff).count(True)

  def n_holes_below_cutoff (self, cutoff) :
    assert (cutoff < 0)
    return (self.holes.heights < cutoff).count(True)

  def save_pdb_file (self,
      file_name="peaks.pdb",
      include_holes=True,
      include_anom=True,
      include_water=True,
      log=None) :
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
    def create_atom (xyz, peak, serial) :
      rg = iotbx.pdb.hierarchy.residue_group(resseq=str(serial))
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
    for peak, xyz in zip(self.peaks.heights, self.peaks.sites) :
      rg = create_atom(xyz, peak, k)
      peaks_chain.append_residue_group(rg)
      k += 1
    f = open(file_name, "w")
    f.write("REMARK  Interesting sites from mmtbx.find_peaks_holes\n")
    f.write("REMARK  Chain A is mFo-DFc peaks (> %g sigma)\n" % self.map_cutoff)
    if (include_holes) :
      f.write("REMARK  Chain B is mFo-DFc holes (< -%g sigma)\n" %
        (- self.map_cutoff))
      holes_chain = iotbx.pdb.hierarchy.chain(id="B")
      model.append_chain(holes_chain)
      k = 1
      for hole, xyz in zip(self.holes.heights, self.holes.sites) :
        rg = create_atom(xyz, hole, k)
        holes_chain.append_residue_group(rg)
        k += 1
    if (include_anom) and (self.anom_peaks is not None) :
      f.write("REMARK  Chain C is anomalous peaks (> %g sigma)\n" %
        self.anom_map_cutoff)
      anom_chain = iotbx.pdb.hierarchy.chain(id="C")
      model.append_chain(anom_chain)
      k = 1
      for peak, xyz in zip(self.anom_peaks.heights, self.anom_peaks.sites) :
        rg = create_atom(xyz, peak, k)
        anom_chain.append_residue_group(rg)
        k += 1
    if (include_water) and (self.water_peaks is not None) :
      f.write("REMARK  Chain D is waters with mFo-DFc peaks (> %g sigma)\n" %
        self.map_cutoff)
      waters_chain = iotbx.pdb.hierarchy.chain(id="D")
      model.append_chain(waters_chain)
      for k, peak in enumerate(self.water_peaks) :
        rg = create_atom(peak.xyz, peak.peak_height, k+1)
        waters_chain.append_residue_group(rg)
      if (include_anom) and (self.water_anom_peaks is not None) :
        f.write("REMARK  Chain E is waters with anom. peaks (> %g sigma)\n" %
          self.anom_map_cutoff)
        waters_chain_2 = iotbx.pdb.hierarchy.chain(id="E")
        model.append_chain(waters_chain_2)
        for k, peak in enumerate(self.water_anom_peaks) :
          rg = create_atom(peak.xyz, peak.peak_height, k+1)
          waters_chain_2.append_residue_group(rg)
    f.write(root.as_pdb_string())
    f.close()
    print >> log, "Wrote %s" % file_name
    self.pdb_file = file_name

  def get_output_file_info (self) :
    output_files = []
    if (self.pdb_file is not None) :
      output_files.append((self.pdb_file, "Peaks as PDB atoms"))
    if (self.map_file is not None) :
      output_files.append((self.map_file, "Map coefficients"))
    return output_files

class water_peak (object) :
  def __init__ (self, id_str, xyz, peak_height, map_type="mFo-DFc") :
    adopt_init_args(self, locals())

  def show (self, out=sys.stdout) :
    print >> out, "  %s  map_type=%s  peak=%g" % (self.id_str,
      self.map_type, self.peak_height)

def find_peaks_holes (
    fmodel,
    pdb_hierarchy,
    params=None,
    map_cutoff=3.0,
    anom_map_cutoff=3.0,
    filter_peaks_by_2fofc=None,
    out=None) :
  """
  Find peaks and holes in mFo-DFc map, plus flag solvent atoms with
  suspiciously high mFo-DFc values, plus anomalous peaks if anomalous data are
  present.  Returns a pickle-able object storing all this information (with
  the ability to write out a PDB file with the sites of interest).
  """
  if (out is None) : out = sys.stdout
  if (params is None) :
    params = master_phil.fetch().extract().find_peaks
  pdb_atoms = pdb_hierarchy.atoms()
  unit_cell = fmodel.xray_structure.unit_cell()
  from mmtbx import find_peaks
  from cctbx import maptbx
  f_map = None
  if (filter_peaks_by_2fofc is not None) :
    f_map_ = fmodel.electron_density_map().fft_map(
      resolution_factor=params.resolution_factor,
      symmetry_flags=maptbx.use_space_group_symmetry,
      map_type="2mFo-DFc",
      use_all_data=True)
    f_map_.apply_sigma_scaling()
    f_map = f_map_.real_map()
  make_header("Positive difference map peaks", out=out)
  peaks_result = find_peaks.manager(
    fmodel=fmodel,
    map_type="mFo-DFc",
    map_cutoff=map_cutoff,
    params=params,
    log=out)
  peaks_result.peaks_mapped()
  peaks_result.show_mapped(pdb_atoms)
  peaks = peaks_result.peaks()
  if (filter_peaks_by_2fofc is not None) :
    n_removed = peaks.filter_by_secondary_map(
      map=f_map,
      min_value=filter_peaks_by_2fofc)
    print >> out, ""
    print >> out, "%d peaks remaining after 2mFo-DFc filtering" % \
      len(peaks.sites)
  # very important - sites are initially fractional coordinates!
  peaks.sites = unit_cell.orthogonalize(peaks.sites)
  print >> out, ""
  out.flush()
  make_header("Negative difference map holes", out=out)
  holes_result = find_peaks.manager(
    fmodel=fmodel,
    map_type="mFo-DFc",
    map_cutoff=-map_cutoff,
    params=params,
    log=out)
  holes_result.peaks_mapped()
  holes_result.show_mapped(pdb_atoms)
  holes = holes_result.peaks()
  # XXX is this useful?
  #if (filter_peaks_by_2fofc is not None) :
  #  holes.filter_by_secondary_map(
  #    map=f_map,
  #    min_value=filter_peaks_by_2fofc)
  holes.sites = unit_cell.orthogonalize(holes.sites)
  print >> out, ""
  out.flush()
  anom = None
  if (fmodel.f_obs().anomalous_flag()) :
    make_header("Anomalous difference map peaks", out=out)
    anom_result = find_peaks.manager(
      fmodel=fmodel,
      map_type="anomalous",
      map_cutoff=anom_map_cutoff,
      params=params,
      log=out)
    anom_result.peaks_mapped()
    anom_result.show_mapped(pdb_atoms)
    anom = anom_result.peaks()
    if (filter_peaks_by_2fofc is not None) :
      anom.filter_by_secondary_map(
        map=f_map,
        min_value=filter_peaks_by_2fofc)
      print >> out, ""
      print >> out, "%d peaks remaining after 2mFo-DFc filtering" % \
        len(anom.sites)
    anom.sites = unit_cell.orthogonalize(anom.sites)
    print >> out, ""
    out.flush()
  cache = pdb_hierarchy.atom_selection_cache()
  water_isel = cache.selection("resname HOH").iselection()
  waters_out = [None, None]
  if (len(water_isel) > 0) :
    sites_frac = fmodel.xray_structure.sites_frac()
    map_types = ["mFo-DFc"]
    if (fmodel.f_obs().anomalous_flag()) :
      map_types.append("anomalous")
    for k, map_type in enumerate(map_types) :
      fft_map = fmodel.electron_density_map().fft_map(
        resolution_factor=params.resolution_factor,
        symmetry_flags=maptbx.use_space_group_symmetry,
        map_type=map_type,
        use_all_data=True)
      fft_map.apply_sigma_scaling()
      real_map = fft_map.real_map_unpadded()
      suspicious_waters = []
      for i_seq in water_isel :
        atom = pdb_atoms[i_seq]
        rho = real_map.tricubic_interpolation(sites_frac[i_seq])
        if (rho >= map_cutoff) :
          peak = water_peak(
            id_str=atom.id_str(),
            xyz=atom.xyz,
            peak_height=rho,
            map_type=map_type)
          suspicious_waters.append(peak)
      if (len(suspicious_waters) > 0) :
        make_header("Water molecules with %s peaks" % map_type, out=out)
        for peak in suspicious_waters :
          peak.show(out=out)
        print >> out, ""
        waters_out[k] = suspicious_waters
  all_results = peaks_holes_container(
    peaks=peaks,
    holes=holes,
    anom_peaks=anom,
    map_cutoff=map_cutoff,
    anom_map_cutoff=anom_map_cutoff,
    water_peaks=waters_out[0],
    water_anom_peaks=waters_out[1])
  all_results.show_summary(out=out)
  return all_results

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  if (len(args) == 0) :
    phil_out = StringIO()
    master_phil.show(out=phil_out)
    raise Usage("""
mmtbx.find_peaks_holes - difference map analysis
  Prints a summary of all peaks and holes above the specified cutoff in the
  mFo-DFc map, and flag any water molecules with suspiciously high peaks
  (possible ions).  Will also check the anomalous map if available.

%s""" % phil_out.getvalue())
  cmdline = utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    create_fmodel=True)
  out.flush()
  result = find_peaks_holes(
    fmodel=cmdline.fmodel,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    params=cmdline.params.find_peaks,
    map_cutoff=cmdline.params.map_cutoff,
    anom_map_cutoff=cmdline.params.anom_map_cutoff,
    filter_peaks_by_2fofc=cmdline.params.filter_peaks_by_2fofc,
    out=out)
  prefix = cmdline.params.output_file_prefix
  if (cmdline.params.write_pdb) :
    result.save_pdb_file(file_name="%s.pdb" % prefix, log=out)
  if (cmdline.params.write_maps) :
    import mmtbx.maps.utils
    import iotbx.mtz
    f_map, diff_map = mmtbx.maps.utils.get_maps_from_fmodel(cmdline.fmodel,
      use_filled=False)
    dec = iotbx.mtz.label_decorator(phases_prefix="PH")
    mtz_dat = f_map.as_mtz_dataset(
      column_root_label="2FOFCWT",
      label_decorator=dec)
    mtz_dat.add_miller_array(diff_map,
      column_root_label="FOFCWT",
      label_decorator=dec)
    if (cmdline.fmodel.f_obs().anomalous_flag()) :
      anom_map = mmtbx.maps.utils.get_anomalous_map(cmdline.fmodel)
      mtz_dat.add_miller_array(anom_map,
        column_root_label="ANOM",
        label_decorator=dec)
    mtz_dat.mtz_object().write("%s_maps.mtz" % prefix)
    result.map_file = "%s_maps.mtz" % prefix
  return result

def validate_params (params, callback=None) :
  if params.input.pdb.file_name is None :
    raise Sorry("No PDB file defined.")
  elif params.input.xray_data.file_name is None :
    raise Sorry("No reflection file defined.")
  elif params.input.xray_data.labels is None :
    raise Sorry("No labels chosen for reflection data.")

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.chdir(self.output_dir)
    return run(args=list(self.args), out=sys.stdout)

def finish_job (result) :
  output_files = []
  stats = []
  if (result is not None) :
    output_files = result.get_output_file_info()
  return (output_files, stats)

if (__name__ == "__main__") :
  run(sys.argv[1:])
