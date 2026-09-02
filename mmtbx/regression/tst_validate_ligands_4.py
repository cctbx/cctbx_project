from __future__ import absolute_import, division, print_function
'''
Tests for the cryo-EM (real-space map) branch of validate_ligands.

The maps are simulated: Fcalc from a model, FFT'd and sigma-scaled, wrapped in a
map_manager. Building the map from a *perturbed* model (ligand deleted) while
the tool correlates against the *correct* model is what makes the RSCC
discriminating, so these tests check both ends of the scale, not just that a
number comes out.
'''
import os
import time
import iotbx.pdb
import mmtbx.model
from cctbx import maptbx, miller
from iotbx.map_manager import map_manager as map_manager_class
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from mmtbx.regression.tst_validate_ligands import find_lr
from mmtbx.validation import validate_ligands as vlmod

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)  # Only show critical errors

D_MIN = 3.0
LIG_SEL = 'chain A and resseq 1 and resname GOL'

# ------------------------------------------------------------------------------

def run():
  run_test28()
  run_test29()
  run_test30()
  run_test31()
  run_test32()
  run_test33()

# ------------------------------------------------------------------------------

# A glycerol with a short poly-ALA alongside it: the peptide gives the ligand an
# environment, so the 'sites' RSCC has something to be computed over. ALA 100 CB
# sits ~3.1 A from GOL O2, inside the default within_radius of 3.0 A.
_map_pdb_str = '''
CRYST1   30.000   40.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  C1  GOL A   1       5.578   9.079   8.959  1.00 20.00           C
HETATM    2  C2  GOL A   1       5.404  10.193   9.989  1.00 20.00           C
HETATM    3  C3  GOL A   1       4.003  10.183  10.608  1.00 20.00           C
HETATM    4  O1  GOL A   1       5.482   7.794   9.563  1.00 20.00           O
HETATM    5  O2  GOL A   1       5.628  11.463   9.370  1.00 20.00           O
HETATM    6  O3  GOL A   1       3.905  11.289  11.512  1.00 20.00           O
ATOM      7  N   ALA A 100       7.700  15.050   8.400  1.00 20.00           N
ATOM      8  CA  ALA A 100       6.400  15.450   8.900  1.00 20.00           C
ATOM      9  C   ALA A 100       5.700  16.250   7.900  1.00 20.00           C
ATOM     10  O   ALA A 100       4.900  17.150   8.200  1.00 20.00           O
ATOM     11  CB  ALA A 100       5.628  14.350   9.370  1.00 20.00           C
ATOM     12  N   ALA A 101       6.100  15.950   6.650  1.00 20.00           N
ATOM     13  CA  ALA A 101       5.500  16.750   5.600  1.00 20.00           C
ATOM     14  C   ALA A 101       6.400  17.950   5.300  1.00 20.00           C
ATOM     15  O   ALA A 101       7.500  17.950   5.850  1.00 20.00           O
ATOM     16  CB  ALA A 101       5.300  15.900   4.350  1.00 20.00           C
ATOM     17  N   ALA A 102       6.000  18.950   4.500  1.00 20.00           N
ATOM     18  CA  ALA A 102       6.800  20.150   4.200  1.00 20.00           C
ATOM     19  C   ALA A 102       6.000  21.150   3.400  1.00 20.00           C
ATOM     20  O   ALA A 102       4.800  21.050   3.300  1.00 20.00           O
ATOM     21  CB  ALA A 102       7.300  20.800   5.480  1.00 20.00           C
'''

# ------------------------------------------------------------------------------

def _model(pdb_str=_map_pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  model.setup_scattering_dictionaries(scattering_table='electron')
  return model


def _simulated_map_manager(model, omit_selection=None, d_min=D_MIN):
  '''Fcalc map from the model, optionally with `omit_selection` deleted.

  The gridding is the same either way, so the resulting map_manager is directly
  comparable with one built from the full model.
  '''
  xrs = model.get_xray_structure()
  cg = maptbx.crystal_gridding(
    unit_cell        = xrs.unit_cell(),
    space_group_info = xrs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.6)
  if omit_selection is not None:
    xrs = xrs.select(~model.selection(omit_selection))
  fft_map = miller.fft_map(
    crystal_gridding     = cg,
    fourier_coefficients = xrs.structure_factors(d_min=d_min).f_calc())
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  return map_manager_class(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = True)


def _map_manager_and_vl(omit_selection=None, resolution=D_MIN):
  model = _model()
  mm = _simulated_map_manager(model, omit_selection=omit_selection)
  params = vlmod.master_params().extract().validate_ligands
  params.ligand_code = []
  params.resolution = resolution
  vl = vlmod.manager(model=model, fmodel=None, map_manager=mm, params=params,
                     log=null_out())
  vl.run()
  return model, mm, vl

# ------------------------------------------------------------------------------

def run_test28():
  '''A ligand that is genuinely in the map scores high, and the map path fills
  in the site and per-fragment CCs rather than leaving them None.'''
  print('test28')
  model, mm, vl = _map_manager_and_vl()
  lr = find_lr(vl, LIG_SEL)
  ccs = lr.get_ccs()
  assert ccs is not None
  assert ccs.rscc > 0.9, ccs.rscc
  assert ccs.rscc_sites is not None
  assert ccs.rscc_sites > 0.9, ccs.rscc_sites
  assert ccs.frag_ccs, ccs.frag_ccs
  assert len(ccs.frag_ccs) == len(lr.ligand_rigid_components_isels)

# ------------------------------------------------------------------------------

def run_test29():
  '''A ligand that is not in the map at all must score below the highlight
  cutoff; the mask must not be so wide that neighbouring density carries it.'''
  print('test29')
  model, mm, vl = _map_manager_and_vl(omit_selection=LIG_SEL)
  lr = find_lr(vl, LIG_SEL)
  ccs = lr.get_ccs()
  assert ccs is not None
  assert ccs.rscc < 0.6, ccs.rscc
  # the environment is present in that map, so its RSCC stays high:
  assert ccs.rscc_sites > 0.9, ccs.rscc_sites

# ------------------------------------------------------------------------------

def run_test30():
  '''Reflections and a map together should raise a Sorry.'''
  print('test30')
  model = _model()
  mm = _simulated_map_manager(model)
  params = vlmod.master_params().extract().validate_ligands
  params.resolution = D_MIN
  try:
    vlmod.manager(model=model, fmodel='not-an-fmodel', map_manager=mm,
                  params=params, log=null_out())
  except Sorry as e:
    assert 'map' in str(e).lower(), str(e)
  else:
    raise AssertionError('expected Sorry for fmodel + map_manager')

# ------------------------------------------------------------------------------

def run_test31():
  '''The two mask radii are deliberately different; guard against silent drift.'''
  print('test31')
  assert approx_equal(vlmod.CC_MASK_RADIUS_MAP, 2.0)
  assert approx_equal(vlmod.CC_MASK_RADIUS_XRAY, 1.5)

# ------------------------------------------------------------------------------

def run_test32():
  '''End to end through the Program: a map on the command line works, and the
  resolution is taken from the map when it is not supplied.'''
  print('test32')
  from iotbx.cli_parser import run_program
  from mmtbx.programs import validate_ligands as val_lig

  model = _model()
  mm = _simulated_map_manager(model)
  pdb_fn = 'tst_validate_ligands_4.pdb'
  map_fn = 'tst_validate_ligands_4.map'
  with open(pdb_fn, 'w') as fh:
    fh.write(_map_pdb_str)
  mm.write_map(map_fn)
  try:
    result = run_program(
      program_class = val_lig.Program,
      args          = [pdb_fn, map_fn, 'run_reduce2=False'],
      logger        = null_out())
    lr = find_lr(result.ligand_manager, LIG_SEL)
    ccs = lr.get_ccs()
    assert ccs is not None and ccs.rscc is not None
    assert ccs.rscc > 0.9, ccs.rscc
  finally:
    for fn in (pdb_fn, map_fn):
      if os.path.isfile(fn):
        os.remove(fn)

# ------------------------------------------------------------------------------

def run_test33():
  '''A map plus reflection data is rejected up front.'''
  print('test33')
  from iotbx.cli_parser import run_program
  from mmtbx.programs import validate_ligands as val_lig

  model = _model()
  mm = _simulated_map_manager(model)
  xrs = model.get_xray_structure()
  f_obs = abs(xrs.structure_factors(d_min=D_MIN).f_calc())
  f_obs = f_obs.set_observation_type_xray_amplitude()

  pdb_fn = 'tst_validate_ligands_4b.pdb'
  map_fn = 'tst_validate_ligands_4b.map'
  mtz_fn = 'tst_validate_ligands_4b.mtz'
  with open(pdb_fn, 'w') as fh:
    fh.write(_map_pdb_str)
  mm.write_map(map_fn)
  f_obs.as_mtz_dataset(column_root_label='FOBS').mtz_object().write(mtz_fn)
  try:
    try:
      run_program(
        program_class = val_lig.Program,
        args          = [pdb_fn, map_fn, mtz_fn, 'run_reduce2=False'],
        logger        = null_out())
    except Sorry as e:
      assert 'map' in str(e).lower(), str(e)
    else:
      raise AssertionError('expected Sorry for map + reflection data')
  finally:
    for fn in (pdb_fn, map_fn, mtz_fn):
      if os.path.isfile(fn):
        os.remove(fn)

# ------------------------------------------------------------------------------

if (__name__ == '__main__'):
  t0 = time.time()
  run()
  print('OK. Time: %8.3f' % (time.time() - t0))
