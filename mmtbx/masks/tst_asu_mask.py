import sys, os, time
import iotbx.pdb
from scitbx.array_family import flex
import iotbx
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.shelx import fcf
from cStringIO import StringIO
import mmtbx.f_model
from mmtbx import utils
from mmtbx import masks
from cctbx import adptbx
from cctbx import miller
from cctbx.sgtbx import space_group_info
from cctbx.development import random_structure
from cctbx import crystal
from cctbx import xray


# from
#    Phenix/cctbx/miller.py
def r1_factor(obs, other, scale=1.0, assume_index_matching=False):
    """ sum ||F| - |F'|| / sum |F|
    where F is obs.data() and F' is other.data() """
    assert (obs.observation_type() is None
            or obs.is_complex_array() or obs.is_xray_amplitude_array())
    assert (obs.observation_type() is None
            or other.is_complex_array() or other.is_xray_amplitude_array())
    if assume_index_matching:
      data0, data = obs.data(), other.data()
    else:
      matching = miller.match_indices(obs.indices(), other.indices())
      # assert not matching.have_singles()
      data0 = obs.select(matching.pairs().column(0)).data()
      data = other.select(matching.pairs().column(1)).data()
      assert data0.size() == data.size() & data.size() > 1
    data  = scale*flex.abs(data)
    data0 = flex.abs(data0)
    return flex.sum(flex.abs(data - data0)) / flex.sum(data0)

def scale_factor(obs, other, assume_index_matching=False):
    """ |F| ~ k*|F'|   k = sum (|F| * |F'|) / sum (|F'|*|F'|)
    where F is obs.data() and F' is other.data() """
    assert (obs.observation_type() is None
            or obs.is_complex_array() or obs.is_xray_amplitude_array())
    assert (obs.observation_type() is None
            or other.is_complex_array() or other.is_xray_amplitude_array())
    if assume_index_matching:
      data0, data = obs.data(), other.data()
    else:
      matching = miller.match_indices(obs.indices(), other.indices())
      # assert not matching.have_singles()
      data0 = obs.select(matching.pairs().column(0)).data()
      data = other.select(matching.pairs().column(1)).data()
      assert data0.size() == data.size() & data.size() > 1
    data  = flex.abs(data)
    data0 = flex.abs(data0)
    scale = flex.sum(data * data0) / flex.sum(data*data)
    return  scale


def get_radii(structure):
  from cctbx.eltbx import van_der_waals_radii
  unknown = []
  atom_radii = []
  for i_seq, scatterer in enumerate(structure.scatterers()):
    try:
      atom_radii.append(
           van_der_waals_radii.vdw.table[scatterer.element_symbol()])
    except:
      unknown.append(scatterer.element_symbol())
  return atom_radii



Elements = ("N", "C", "O", "H", "Ca", "C", "C", "O", "O", "N", "H", "H", "Mg", "Se")

def build_struc():
  cell = (4.57428, 5.94656, 7.77627, 90, 109, 90)
  site = (0.266883, 0.568905, 0.291467)
  sg = "Cc"
  symmetry = crystal.symmetry(unit_cell=cell,
                              space_group_symbol=sg)
  structure = xray.structure(crystal_symmetry=symmetry)
  scatterer = xray.scatterer(
                 site = site,
                 u = 0.1,
                 occupancy = 1.0,
                 scattering_type = "C")
  structure.add_scatterer(scatterer)
  return structure

def exercize_1(sg, atoms, molvol):
  params = mmtbx.masks.mask_master_params.extract()
  params.ignore_hydrogens = False
  params.ignore_zero_occupancy_atoms = False
  # params.shrink_truncation_radius = 0.0
  # params.solvent_radius = 0.0
  group = space_group_info(sg)
  group.show_summary()
  struc = random_structure.xray_structure(
      space_group_info = group,
      # unit_cell = (14.6, 26.1, 29.2, 90, 90, 90),
      volume_per_atom = molvol,
      general_positions_only = True,
      elements = atoms
      )
  struc.show_summary()
  fc = struc.structure_factors(d_min = 2.10145471291).f_calc()
  struc_p1 = struc.expand_to_p1()
  fc_p1 = struc_p1.structure_factors(d_min = 2.10145471291).f_calc()
  asu_mask = masks.atom_mask(
      unit_cell = struc.unit_cell(),
      group = struc.space_group(),
      resolution = fc.d_min(),
      grid_method = 2,
      grid_step_factor = params.grid_step_factor,
      solvent_radius = params.solvent_radius,
      shrink_truncation_radius = params.shrink_truncation_radius )
  grid =  asu_mask.data.focus()
  print "asu mask grid = ", grid
  blk_p1 = masks.bulk_solvent(
    xray_structure = struc_p1,
    gridding_n_real = grid,
    ignore_zero_occupancy_atoms = params.ignore_zero_occupancy_atoms,
    ignore_hydrogen_atoms = params.ignore_hydrogens,
    solvent_radius = params.solvent_radius,
    shrink_truncation_radius = params.shrink_truncation_radius)
  blk_p1.show_summary()
  fm_p1 = blk_p1.structure_factors( miller_set = fc_p1 )
  radii = get_radii(struc)
  assert len(radii) == len(struc.sites_frac())
  asu_mask.compute( struc.sites_frac(), radii )
  fm_asu = asu_mask.structure_factors( fc.indices() )
  fm_asu = fc.set().array( data =  fm_asu )
  print " solvent content = ", asu_mask.contact_surface_fraction*100.0
  k_asu = scale_factor( fm_asu, fm_p1 )
  print "Asu-P1  Scale-factor = ", k_asu
  r1_asu = r1_factor( fm_asu, fm_p1 )
  print "Asu-P1  R-factor = ", r1_asu
  r1_asu_k = r1_factor(fm_asu, fm_p1, scale=k_asu)
  print "Asu-P1 Scaled R-factor = ", r1_asu_k
  print len(asu_mask.data)
  print asu_mask.data.size()
  assert r1_asu_k < 1.0E-6


def run():
  sg = "P 21 21 21"
  natoms = 0
  molvol = 50.0
  atoms = ("C", "H", "O")
  if len(sys.argv)>1 :
    sg = sys.argv[1]
  if len(sys.argv)>2 :
    natoms = int(sys.argv[2])
  if len(sys.argv)>3 :
    molvol = float(sys.argv[3])
  if natoms>0 :
    atoms = []
    for i in xrange(natoms):
      atoms.append("C")
  if sg == "std":
    exercize_0()
  if sg == "all" :
    for isg in xrange(1,231):
      exercize_1(str(isg), atoms, molvol)
  else:
    exercize_1(sg, atoms, molvol)

if (__name__ == "__main__"):
  run()

