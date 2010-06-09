import os, sys
import libtbx.load_env
from libtbx.utils import Sorry
from libtbx.test_utils import Exception_expected
from iotbx.file_reader import any_file
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
import cPickle

if (libtbx.env.has_module("ccp4io")):
  from iotbx import mtz
else:
  mtz = None

def exercise_others () :
  #--- CIF
  cif_data = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
HOH      .   'water                               ' solvent             3   1 .
#
data_comp_HOH
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 HOH           O      O    OH2      -0.408
 HOH           H1     H    HOH2      0.204
 HOH           H2     H    HOH2      0.204
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 HOH      O      n/a    .      END
 HOH      H1     O      .      .
 HOH      H2     O      .      .
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 HOH      O      H1        coval       0.980    0.020
 HOH      O      H2        coval       0.980    0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 HOH      H1     O      H2      106.800    3.000
"""

  f = open("tmp1.cif", "w")
  f.write(cif_data)
  f.close()
  cif = any_file("tmp1.cif")
  cif.assert_file_type("cif")
  os.remove("tmp1.cif")

  #--- PDB
  pdb_data = """
HEADER    HYDROLASE(METALLOPROTEINASE)            17-NOV-93   1THL
ATOM      1  N   ILE     1       9.581  51.813  -0.720  1.00 31.90      1THL 158
ATOM      2  CA  ILE     1       8.335  52.235  -0.041  1.00 52.95      1THL 159
ATOM      3  C   ILE     1       7.959  53.741   0.036  1.00 26.88      1THL 160
END
"""

  f = open("tmp1.pdb", "w")
  f.write(pdb_data)
  f.close()
  pdb = any_file("tmp1.pdb")
  try :
    pdb.assert_file_type("txt")
  except Sorry :
    pass
  else :
    raise Exception_expected(
      "Expected exception for 'pdb.assert_file_type(\"txt\")'.")
  pdb.assert_file_type("pdb")
  os.remove("tmp1.pdb")

  #--- PHIL
  phil_data = """\
refinement {
  input {
    include scope mmtbx.utils.neutron_data_str
    include scope mmtbx.utils.xray_data_str
  }
  refine {
    strategy = *individual_adp *individual_sites *tls occupancies group_adp
      .type = choice(multi=True)
    adp {
      tls = None
        .type = str
        .multiple = True
    }
  }
  main {
    number_of_macro_cycles = 3
      .type = int
    simulated_annealing = True
      .type = bool
  }
}"""

  f = open("tmp1.phil", "w")
  f.write(phil_data)
  f.close()
  phil = any_file("tmp1.phil")
  assert phil.file_type == "phil"
  params = phil.file_object.extract()
  assert params.refinement.main.number_of_macro_cycles == 3
  assert params.refinement.input.xray_data.r_free_flags.test_flag_value == None
  os.remove("tmp1.phil")

  #--- Pickle
  f = open("tmp1.pkl", "wb")
  cPickle.dump(params, f)
  f.close()
  pkl = any_file("tmp1.pkl")
  assert pkl.file_type == "pkl"
  params = pkl.file_object
  assert params.refinement.main.number_of_macro_cycles == 3
  os.remove("tmp1.pkl")

def exercise_maps () :
  xplor_map = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/cns.map",
    test=os.path.isfile)
  if xplor_map is not None :
    f = any_file(xplor_map)
    assert f.file_type == "xplor_map"
  ccp4_map = libtbx.env.under_dist(
    module_name="iotbx",
    path="ccp4_map/tst_input.map")
  f = any_file(ccp4_map)
  assert f.file_type == "ccp4_map"

def exercise_hkl () :
  #--- HKL
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  input_arrays = [miller_set.array(
    data=flex.random_double(size=miller_set.indices().size()))
      .set_observation_type_xray_amplitude()
        for i in [0,1]]
  mtz_dataset = input_arrays[0].as_mtz_dataset(column_root_label="F0")
  mtz_dataset.mtz_object().write("tmp1.mtz")
  hkl = any_file("tmp1.mtz")
  assert hkl.file_type == "hkl"
  #assert hkl.file_server
  assert hkl.file_server.miller_arrays[0].info().labels == ["F0"]
  os.remove("tmp1.mtz")

def exercise () :
  exercise_others()
  exercise_maps()
  if mtz is None :
    print "Skipping mtz file tests"
  else :
    exercise_hkl()
  print "OK"
if __name__ == "__main__" :
  exercise()
#---end
