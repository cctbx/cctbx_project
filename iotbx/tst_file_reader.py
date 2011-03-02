import os
import libtbx.load_env
from libtbx.utils import Sorry
from libtbx.test_utils import Exception_expected
from iotbx.file_reader import any_file, sort_by_file_type, group_files
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
  pdb.check_file_type("pdb")
  try :
    pdb.check_file_type("hkl")
  except Sorry :
    pass
  else :
    raise Exception_expected
  try :
    pdb.check_file_type(multiple_formats=["hkl","seq"])
  except Sorry :
    pass
  else :
    raise Exception_expected
  pdb.check_file_type(multiple_formats=["hkl","pdb"])
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

  #--- sequence
  seqs = ["AAKDVKFGVNVLADAV",
          "AAKDVKFGNDARVKML",
          "AVKDVRYGNEARVKIL"]
  headers = ["anb.pdb chain A",
             "anb.pdb chain B",
             "anb.pdb chain C"]
  f = open("sequence.pir", "w")
  f.write("""\
> %s

%s
*""" % (headers[0], seqs[0]))
  f.close()
  f = any_file("sequence.pir")
  assert (f.file_type == "seq")
  assert (f.file_object[0].sequence == (seqs[0] + "*"))
  os.remove("sequence.pir")
  f = open("sequence.dat", "w")
  f.write(seqs[0])
  f.close()
  f = any_file("sequence.dat")
  assert (f.file_type == "seq")
  assert (f.file_object[0].sequence == seqs[0])
  os.remove("sequence.dat")
  f = open("sequence.fa", "w")
  f.write("""\
> %s
%s""" % (headers[0], seqs[0]))
  f.close()
  f = any_file("sequence.fa", "w")
  assert (f.file_type == "seq")
  os.remove("sequence.fa")
  f = open("sequences.fa", "w")
  for header, seq in zip(headers, seqs) :
    f.write("""\
> %s
%s
""" % (header, seq))
  f.close()
  f = any_file("sequences.fa")
  assert (f.file_type == "seq")
  for i, seq_object in enumerate(f.file_object) :
    assert (seq_object.sequence == seqs[i])
  os.remove("sequences.fa")

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

def exercise_misc () :
  file_names = ["foo.pdb", "foo.mtz", "bar.pdb", "bar.mtz", "seq.dat"]
  file_names = sort_by_file_type(file_names, sort_order=["pdb","hkl","seq"])
  assert (file_names == ['foo.pdb','bar.pdb','foo.mtz','bar.mtz','seq.dat'])

def exercise_groups () :
  files = [
    "/data/user/project1/p1.sca",
    "/data/user/project2/p2.pdb",
    "/data/user/project1/p1.pdb",
    "/data/user/project1/process/mosflm.log",
    "/data/user/project4/refine/current.pdb",
    "/data/user/project3/data.sca",
    "/data/user/project3/model.pdb",
    "/data/user/project2/native.mtz",
    "/data/user/project3/ligands.cif",
    "/data/user/project4/data.mtz",
    "/data/user/project4/refine/ligands.cif",
    #"/data/user/project2/p2_reindex.pdb",
    #"/data/user/project2/p2_reindex.mtz",
    "/data/user/new_data.hkl",
    "/data/user/project5/model.pdb",
    "/data/user/project5/model2.pdb",
    "/data/user/project5/anom.mtz",
  ]
  g = group_files(files)
  assert (g.ambiguous_files == ['/data/user/new_data.hkl',
                                '/data/user/project5/anom.mtz'])
  assert (g.ungrouped_files == [])
  assert (g.grouped_files == [
    ['/data/user/project2/p2.pdb',
      '/data/user/project2/native.mtz'],
    ['/data/user/project1/p1.pdb',
      '/data/user/project1/p1.sca',
      '/data/user/project1/process/mosflm.log'],
    ['/data/user/project4/refine/current.pdb',
      '/data/user/project4/data.mtz',
      '/data/user/project4/refine/ligands.cif'],
    ['/data/user/project3/model.pdb',
      '/data/user/project3/data.sca',
      '/data/user/project3/ligands.cif'],
    ['/data/user/project5/model.pdb'],
    ['/data/user/project5/model2.pdb'],
  ])
  g = group_files(files, template_format="hkl")
  assert (g.ambiguous_files == [])
  assert (g.ungrouped_files == [])
  assert (g.grouped_files == [
    ['/data/user/project1/p1.sca',
      '/data/user/project1/p1.pdb',
      '/data/user/project1/process/mosflm.log'],
    ['/data/user/project3/data.sca',
      '/data/user/project3/model.pdb',
      '/data/user/project3/ligands.cif'],
    ['/data/user/project2/native.mtz',
      '/data/user/project2/p2.pdb'],
    ['/data/user/project4/data.mtz',
      '/data/user/project4/refine/current.pdb',
      '/data/user/project4/refine/ligands.cif'],
    ['/data/user/new_data.hkl'],
    ['/data/user/project5/anom.mtz',
      '/data/user/project5/model.pdb',
      '/data/user/project5/model2.pdb']
  ])

def exercise () :
  exercise_groups()
  exercise_misc()
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
