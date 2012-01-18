
# simple frontend to mmtbx.find_peaks, primarily intended for use in quickly
# analyzing structures in the PDB (and storing results)
#
# TODO (???) plugin for Coot

from mmtbx import utils
from libtbx.str_utils import make_header
import libtbx.phil
from libtbx import adopt_init_args, group_args
from cStringIO import StringIO
import sys

master_phil = libtbx.phil.parse("""
%s
find_peaks {
  include scope mmtbx.find_peaks.master_params
}
map_cutoff = 3.0
  .type = float
anom_map_cutoff = 3.0
  .type = float
""" % utils.cmdline_input_phil_str,
  process_includes=True)

class peaks_holes_container (object) :
  def __init__ (self, peaks, holes, map_cutoff=3.0, anom_peaks=None,
      anom_map_cutoff=3.0, water_peaks=None, water_anom_peaks=None) :
    adopt_init_args(self, locals())

  def show_summary (self, out=sys.stdout) :
    print >> out, ""
    print >> out, "SUMMARY OF MAP PEAKS:"
    cutoffs = [self.map_cutoff, self.map_cutoff + 2.0, self.map_cutoff + 4.0]
    for cutoff in cutoffs :
      n_peaks = (self.peaks.heights > cutoff).count(True)
      print >> out, "  mFo-DFc >  %g  : %6d" % (cutoff, n_peaks)
    for cutoff in cutoffs :
      n_holes = (self.holes.heights < -cutoff).count(True)
      print >> out, "  mFo-DFc < -%g  : %6d" % (cutoff, n_holes)
    if (self.anom_peaks is not None) :
      print >> out, "  anomalous > %g : %6d" % (self.anom_map_cutoff,
        len(self.anom_peaks.heights))
    if (self.water_peaks is not None) :
      print >> out, "  suspicious H2O (mFo-DFC > %g) : %6d" % (self.map_cutoff,
        len(self.water_peaks))
    if (self.water_anom_peaks is not None) :
      print >> out, "  anomalous H2O (anomalous > %g): %6d" % (self.map_cutoff,
        len(self.water_anom_peaks))
    print >> out, ""

  def get_summary (self) :
    n_anom_peaks = None
    if (self.anom_peaks is not None) :
      n_anom_peaks = len(self.anom_peaks.heights)
    n_water_peaks = None
    if (self.water_peaks is not None) :
      n_water_peaks = len(self.water_peaks)
    if (self.water_anom_peaks is not None) :
      n_water_anom_peaks = len(self.water_anom_peaks)
    return group_args(
      n_peaks_1=(self.peaks.heights > self.map_cutoff).count(True),
      n_peaks_2=(self.peaks.heights > self.map_cutoff + 2).count(True),
      n_peaks_3=(self.peaks.heights > self.map_cutoff + 4).count(True),
      n_holes_1=(self.holes.heights < -self.map_cutoff).count(True),
      n_holes_2=(self.holes.heights < -self.map_cutoff - 2).count(True),
      n_holes_3=(self.holes.heights < -self.map_cutoff - 4).count(True),
      n_anom_peaks=n_anom_peaks,
      n_water_peaks=n_water_peaks,
      n_water_anom_peaks=n_water_anom_peaks)

  def n_peaks_above_cutoff (self, cutoff) :
    assert (cutoff > 0)
    return (self.peaks.heights > cutoff).count(True)

  def n_holes_below_cutoff (self, cutoff) :
    assert (cutoff < 0)
    return (self.holes.heights < cutoff).count(True)

class water_peak (object) :
  def __init__ (self, id_str, xyz, peak_height, map_type="mFo-DFc") :
    adopt_init_args(self, locals())

  def show (self, out=sys.stdout) :
    print >> out, "  %s  map_type=%s  peak=%g" % (self.id_str,
      self.map_type, self.peak_height)

def find_peaks_holes (
    fmodel,
    pdb_hierarchy,
    params,
    map_cutoff=3.0,
    anom_map_cutoff=3.0,
    out=None) :
  if (out is None) : out = sys.stdout
  pdb_atoms = pdb_hierarchy.atoms()
  from mmtbx import find_peaks
  from cctbx import maptbx
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
  print >> out, ""
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
  print >> out, ""
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
    print >> out, ""
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
    master_phil.show(f=phil_out)
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
  return find_peaks_holes(
    fmodel=cmdline.fmodel,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    params=cmdline.params.find_peaks,
    map_cutoff=cmdline.params.map_cutoff,
    anom_map_cutoff=cmdline.params.anom_map_cutoff,
    out=out)

if (__name__ == "__main__") :
  run(sys.argv[1:])
