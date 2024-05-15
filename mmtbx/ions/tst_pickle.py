from __future__ import absolute_import, division, print_function
import os
from pickle import loads, dumps
from types import MethodType
import libtbx
from mmtbx import ions
from mmtbx.ions import environment
from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
import mmtbx.utils
import iotbx.phil
from iotbx.data_manager import DataManager

def exercise(prefix="mmtbx_ions_tst_pickle"):
  wavelength = 1.025
  mtz_file, pdb_file = generate_zinc_inputs(
      file_base=prefix,
      anonymize = False)
  null_out = libtbx.utils.null_out()

  dm = DataManager()
  m = dm.get_model(pdb_file)
  xrs = m.get_xray_structure()
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")

  os.remove(pdb_file)
  os.remove(mtz_file)

  xrs.set_inelastic_form_factors(photon = wavelength, table = "sasaki")
  fmo.update_xray_structure(xrs, update_f_calc = True)
  m.process(make_restraints=True)

  phil_str = '''
include scope mmtbx.ions.identify.ion_master_phil
include scope mmtbx.ions.svm.svm_phil_str
'''
  params = iotbx.phil.parse(input_string = phil_str,
    process_includes=True).extract()

  manager = ions.identify.create_manager(
    pdb_hierarchy = m.get_hierarchy(),
    fmodel = fmo,
    geometry_restraints_manager = m.get_restraints_manager().geometry,
    wavelength = wavelength,
    params = params,
    nproc = 1,
    log = null_out
    )

  manager.validate_ions(out = null_out)

  fo_map = manager.get_map("mFo")
  fofc_map = manager.get_map("mFo-DFc")

  for atom_props in manager.atoms_to_props.values():
    chem_env = environment.ChemicalEnvironment(
      atom_props.i_seq,
      manager.find_nearby_atoms(atom_props.i_seq, far_distance_cutoff = 3.5),
      manager
      )
    new_chem_env = loads(dumps(chem_env))
    for attr in dir(chem_env):
      if attr == "atom":
        # The two won't be directly comparable, but we will trust atom_labels is
        # tested fully in its own module
        assert chem_env.atom.id_str() == new_chem_env.atom.id_str()
      elif not attr.startswith("_") and \
        not isinstance(getattr(chem_env, attr), MethodType):
        assert getattr(chem_env, attr) == getattr(new_chem_env, attr)

    scatter_env = environment.ScatteringEnvironment(
      atom_props.i_seq, manager, fo_map, fofc_map
      )
    new_scatter_env = loads(dumps(scatter_env))
    for attr in dir(scatter_env):
      if not attr.startswith("_") and \
        not isinstance(getattr(scatter_env, attr), MethodType):
        assert getattr(scatter_env, attr) == getattr(new_scatter_env, attr)

  del fo_map
  del fofc_map

  print("OK")

if __name__ == "__main__":
  exercise()
