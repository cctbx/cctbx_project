from __future__ import division, print_function
import time, os
import iotbx.pdb
import mmtbx.model
from cctbx import crystal
from cctbx import xray
from cctbx.array_family import flex
from iotbx.data_manager import DataManager
from libtbx.utils import null_out

pdb_str = """
CRYST1    4.870   10.060   30.660  94.85  90.26  99.99 P 1           1
ATOM      1  N   GLN A   1       1.153   1.465  -2.445  1.00  9.04           N
ATOM      2  CA  GLN A   1       1.681   2.444  -1.508  1.00  4.57           C
ATOM      3  C   GLN A   1       1.010   2.351  -0.141  1.00  3.21           C
ATOM      4  O   GLN A   1      -0.214   2.379  -0.045  1.00  5.48           O
ATOM      5  CB  GLN A   1       1.514   3.851  -2.090  1.00  6.45           C
ATOM      6  CG  GLN A   1       2.245   4.927  -1.325  1.00  8.01           C
ATOM      7  CD  GLN A   1       1.858   6.297  -1.796  1.00  5.89           C
ATOM      8  OE1 GLN A   1       0.717   6.720  -1.624  1.00  7.78           O
ATOM      9  NE2 GLN A   1       2.793   6.997  -2.405  1.00  6.95           N
"""

def tst000():
  """
  Test scattering table consistency between model object and xray_structure.

  This function sets up an atomic model from a PDB string and assigns
  different scattering tables ("n_gaussian", "wk1995", "it1992", "electron").
  It then verifies that the scattering table assigned via
  `model.setup_scattering_dictionaries` matches the table retrieved from
  the `xray_structure` object.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If the scattering table retrieved from `xray_structure` does not match
      the expected table.
  """
  for table in ["n_gaussian", "wk1995", "it1992", "electron"]:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
    model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
    # set the scattering table via model object
    model.setup_scattering_dictionaries(scattering_table=table)
    # check that this table is consistent with what can be accessed with xrs
    xrs = model.get_xray_structure()
    assert (xrs.scattering_type_registry_params.table == table)
    assert (xrs.get_scattering_table() == table)

# ------------------------------------------------------------------------------

def tst001():
  """
  Test the scattering table consistency between what is supplied in the data
  manager and the xray_structure retrieved from the fmodel object.

  This function writes atomic model from a PDB string to a file,
  and then uses different scattering tables to create an fmodel object via
  the data manager. It tests if the xray structure object correctly retains
  the assigned scattering table.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If the scattering table in xray_structure obtained from
      fmodel via the data manager does not match the assigned table.

  Notes
  -----
  - The generated model is saved as 'model.pdb'.
  - Structure factors are computed at d_min = 1.0 Å.
  """
  phil_str = '''
  data_manager {
      model {
        file = tst_scattering_table_tst_001_model.pdb
        type = %s
      }
      miller_array {
        file = tst_scattering_table_tst_001_data.mtz
        labels {
          name = FOBS
          type = %s
        }
        labels {
          name = R-free-flags
          type = %s
        }
    }
  }'''
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  pdb_fn = 'tst_scattering_table_tst_001_model.pdb'
  mtz_fn = 'tst_scattering_table_tst_001_data.mtz'
  label_fobs = 'FOBS'
  label_flag = 'R-free-flags'
  with open(pdb_fn,"w") as fo:
    fo.write(model.model_as_pdb())
  for table in ["n_gaussian", "wk1995", "it1992", "electron"]:
    if table in ["n_gaussian", "wk1995", "it1992"]:
      _type = '*x_ray neutron electron'
    else:
      _type = 'x_ray neutron *electron'
    model.setup_scattering_dictionaries(scattering_table=table)
    fc = model.get_xray_structure().structure_factors(d_min=1,
      algorithm="direct").f_calc()
    mtz_dataset = abs(fc).as_mtz_dataset(column_root_label = 'FOBS')
    mtz_dataset.add_miller_array(
      miller_array      = fc.generate_r_free_flags(),
      column_root_label = label_flag)
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = mtz_fn)
    dm_phil = iotbx.phil.parse(phil_str%(_type + ' reference', _type, _type))
    dm = DataManager(['model', 'miller_array'])
    dm.load_phil_scope(dm_phil)
    fmodel = dm.get_fmodel(scattering_table=table)
    #print(fmodel.r_work())
    xrs = fmodel.xray_structure
    assert (xrs.get_scattering_table() == table)
    assert (xrs.scattering_type_registry_params.table == table)
  os.remove(pdb_fn)
  os.remove(mtz_fn)

# ------------------------------------------------------------------------------

def tst002():
  """
  Test scattering table consistency after selections.

  This function sets up an atomic model from a PDB string and assigns
  different scattering tables ("n_gaussian", "wk1995", "it1992", "electron").
  It then verifies that the scattering tables are consistent after xrs.select().

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If the scattering table retrieved after applying selections do not match
      the expected table.
  """
  for table in ["n_gaussian", "wk1995", "it1992", "electron"]:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
    model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
    # get a selection
    sel = model.selection('name N or name CA or name C or name O')
    # copy model
    m1 = model.deep_copy()
    # -----
    # Test 1: check that table is consistent after applying xrs.select()
    # -----
    # set the scattering table via model object
    m1.setup_scattering_dictionaries(scattering_table=table)
    # check that this table is consistent with what can be accessed with xrs
    xrs1 = m1.get_xray_structure()
    assert (xrs1.scattering_type_registry_params.table == table)
    assert (xrs1.get_scattering_table() == table)
    # Check that xrs remembers scattering table after selection
    xrs_selected = xrs1.select(sel)
    assert (xrs_selected.get_scattering_table() == table)
    # -----
    # Test 2: check that table is consistent after deep_copy()
    # -----
    xrs_copy = xrs1.deep_copy_scatterers()
    assert (xrs_copy.get_scattering_table() == table)

# ------------------------------------------------------------------------------

def tst003():
  cs = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer("o", (0.5, 0, 0)),
    xray.scatterer("c", (0, 0, 0))))
  xrs = xray.structure(sp, scatterers)
  # we have not set any scattering table yet, so it should be None
  assert (xrs.get_scattering_table() is None)
  #fc = xrs.structure_factors(d_min=1).f_calc()
  # here we create a default, which is n_gaussian
  xrs.scattering_type_registry()
  assert (xrs.get_scattering_table() == 'n_gaussian')
  # Now we set it explicitly
  xrs.scattering_type_registry(table = 'electron')
  assert (xrs.get_scattering_table() == 'electron')

# ------------------------------------------------------------------------------

if(__name__ == "__main__"):
  """
  Run the scattering table tests.

  This script executes the following tests:
  - `tst000()`: Ensures cconsistency between model object and xray_structure.
  - `tst001()`: consistency between what is supplied in the data manager and
     the xray_structure retrieved from the fmodel object
  - `tst002()`: Ensure consistency after xrs.select()
  - `tst003()`: illustrate toy example

  Returns
  -------
  None
  """
  t0 = time.time()
  tst000()
  tst001()
  tst002()
  tst003() # toy example
  print("OK. Time: %8.3f"%(time.time()-t0))
