import sys, math
from cctbx.array_family import flex
from mmtbx import find_peaks
from mmtbx import utils
import iotbx.phil
from scitbx import matrix

master_params_part1 = iotbx.phil.parse("""\
map_type = mFobs-DFmodel
  .type = str
  .help = Map type to be used to find hydrogens
map_cutoff = 2.0
  .type = float
  .help = Map cutoff
""")

master_params_part2 = find_peaks.master_params.fetch(iotbx.phil.parse("""\
use_sigma_scaled_maps = True
resolution_factor = 1./4.
map_next_to_model
{
  min_model_peak_dist = 0.7
  max_model_peak_dist = 1.2
  min_peak_peak_dist = 1.0
  use_hydrogens = False
}
peak_search
{
  peak_search_level = 1
  min_cross_distance = 1.0
}
"""))

def all_master_params():
  return iotbx.phil.parse("""\
    include scope mmtbx.hydrogens.find.master_params_part1
    include scope mmtbx.hydrogens.find.master_params_part2
""", process_includes=True)

class h_peak(object):
  def __init__(self, site_frac,
                     height,
                     dist,
                     scatterer_o,
                     atom_attribute_o,
                     i_seq_o):
    self.site_frac        = site_frac
    self.height           = height
    self.dist             = dist
    self.scatterer_o      = scatterer_o
    self.atom_attribute_o = atom_attribute_o
    self.i_seq_o          = i_seq_o

def water_bond_angle(o,h1,h2):
  result = None
  a = h1[0]-o[0], h1[1]-o[1], h1[2]-o[2]
  b = h2[0]-o[0], h2[1]-o[1], h2[2]-o[2]
  a = matrix.col(a)
  b = matrix.col(b)
  return a.angle(b, deg=True)

def find_hydrogen_peaks(fmodel,
                        atom_attributes_list,
                        params,
                        log):
  fp_manager = find_peaks.manager(fmodel     = fmodel,
                                  map_type   = params.map_type,
                                  map_cutoff = params.map_cutoff,
                                  params     = params,
                                  log        = log)
  result = fp_manager.peaks_mapped()
  fp_manager.show_mapped(atom_attributes_list = atom_attributes_list)
  return result

def extract_hoh_peaks(peaks, atom_attributes_list, xray_structure, log):
  scatterers = xray_structure.scatterers()
  assert scatterers.size() == len(atom_attributes_list)
  assert peaks.sites.size() == peaks.heights.size()
  assert peaks.heights.size() == peaks.iseqs_of_closest_atoms.size()
  perm = flex.sort_permutation(peaks.iseqs_of_closest_atoms)
  sites = peaks.sites.select(perm)
  heights = peaks.heights.select(perm)
  iseqs_of_closest_atoms = peaks.iseqs_of_closest_atoms.select(perm)
  result = {}
  unit_cell = xray_structure.unit_cell()
  get_class = iotbx.pdb.common_residue_names_get_class
  for s, h, i_seq in zip(sites, heights, iseqs_of_closest_atoms):
    aa = atom_attributes_list[i_seq]
    if(get_class(name = aa.resName) == "common_water"):
      scatterer_o = scatterers[i_seq]
      hp = h_peak(
        site_frac        = s,
        height           = h,
        dist             = unit_cell.distance(s, scatterer_o.site),
        scatterer_o      = scatterer_o,
        atom_attribute_o = aa,
        i_seq_o          = i_seq)
      result.setdefault(i_seq,[]).append(hp)
  for key in result.keys():
    print >> log, atom_attributes_list[int(key)].pdb_format()
    sz = len(result[key])
    for j, jhp in enumerate(result[key]):
      if(sz == 1):
        print >> log, "  ", jhp.height, jhp.dist
      else:
        angles = []
        for k, khp in enumerate(result[key]):
          if(j != k):
            assert jhp.scatterer_o.site == khp.scatterer_o.site
            angl = water_bond_angle(
              o=unit_cell.orthogonalize(jhp.scatterer_o.site),
              h1=unit_cell.orthogonalize(jhp.site_frac),
              h2=unit_cell.orthogonalize(khp.site_frac))
            angles.append(angl)
            print >> log, "  ", jhp.height, jhp.dist, angl
  return result

#def rebuild_hydrogens(h_peaks, model):
#  items = h_peaks.items()
#  items.sort()
#  items.reverse()
#  h_peaks_as_tuple = [(key, value) for key, value in items]
#  # add H ato xray_structure and aal
#  # important: asume atoms are aded from tail to head of the list
#  #            asume only one bonded atom is added
#  sca = model.xray_structure.scatterers()
#  aal = model.atom_attributres_list
#  xh_ct = model
#  for hp in h_peaks_as_tuple:
#    i = hp[0]
#    label = "H1"
#    sc = xray.scatterer(label           = label,
#                        scattering_type = "H",
#                        site            = hp[1].site_frac,
#                        u               = hp.scatterer_o.u_iso, # XXX if aniso, then convert to iso
#                        occupancy       = hp.scatterer_o.occupancy)
#    sca = sca[:i+1]+sc+sca[i+1:]
#    #HETATM   35  H1  HOH S   4       7.424  12.171   7.250  1.00 15.00           H
#    fmt = "HETATM    0  %2s  HOH%2s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f%12s"
#
#  # invalidate relevant memebrs of model

def run(fmodel, model, log):
  params = all_master_params().extract()
  peaks = find_hydrogen_peaks(
    fmodel               = fmodel,
    atom_attributes_list = model.atom_attributes_list,
    params               = params,
    log                  = log)
  h_peaks = extract_hoh_peaks(
    peaks                = peaks,
    atom_attributes_list = model.atom_attributes_list,
    xray_structure       = model.xray_structure,
    log                  = log)

