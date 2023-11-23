from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
from mmtbx.geometry_restraints.torsion_restraints.reference_model import \
    reference_model, reference_model_params
from mmtbx.geometry_restraints.torsion_restraints import utils
from mmtbx.validation.rotalyze import rotalyze
import mmtbx.model
from cctbx.array_family import flex
import iotbx.phil
import iotbx.pdb
from libtbx.test_utils import show_diff
import libtbx.load_env
from six.moves import cStringIO as StringIO
import sys, os, time
from six.moves import range

model_raw_records = """\
CRYST1   41.566   72.307   92.870 108.51  93.02  90.06 P 1           4
ATOM   5466  N   ASN C 236      17.899  72.943  29.028  1.00 60.13           N
ATOM   5467  CA  ASN C 236      16.519  72.435  29.114  1.00 60.52           C
ATOM   5468  C   ASN C 236      16.377  70.925  29.327  1.00 60.49           C
ATOM   5469  O   ASN C 236      15.429  70.294  28.863  1.00 60.60           O
ATOM   5470  CB  ASN C 236      15.689  72.896  27.916  1.00 60.55           C
ATOM   5471  CG  ASN C 236      14.357  73.447  28.338  1.00 61.75           C
ATOM   5472  OD1 ASN C 236      14.256  74.609  28.768  1.00 62.86           O
ATOM   5473  ND2 ASN C 236      13.319  72.616  28.247  1.00 61.22           N
ATOM   5474  N   LEU C 237      17.316  70.364  30.068  1.00 60.55           N
ATOM   5475  CA  LEU C 237      17.444  68.931  30.166  1.00 60.48           C
ATOM   5476  C   LEU C 237      17.815  68.555  31.581  1.00 60.06           C
ATOM   5477  O   LEU C 237      17.335  67.547  32.097  1.00 60.41           O
ATOM   5478  CB  LEU C 237      18.518  68.464  29.178  1.00 60.91           C
ATOM   5479  CG  LEU C 237      18.542  67.095  28.491  1.00 62.25           C
ATOM   5480  CD1 LEU C 237      17.407  66.153  28.923  1.00 63.18           C
ATOM   5481  CD2 LEU C 237      18.563  67.309  26.965  1.00 62.89           C
"""

reference_raw_records = """\
CRYST1   40.688   71.918   93.213 108.16  93.25  90.40 P 1           4
ATOM   5485  N   ASN C 236      16.417  72.834  29.095  1.00  7.17           N
ATOM   5486  CA  ASN C 236      15.051  72.312  29.173  1.00  7.74           C
ATOM   5487  C   ASN C 236      15.000  70.818  29.431  1.00  7.38           C
ATOM   5488  O   ASN C 236      14.047  70.141  29.024  1.00  7.80           O
ATOM   5489  CB  ASN C 236      14.281  72.645  27.887  1.00  8.78           C
ATOM   5490  CG  ASN C 236      12.769  72.657  28.088  1.00 13.44           C
ATOM   5491  OD1 ASN C 236      12.265  73.196  29.082  1.00 20.19           O
ATOM   5492  ND2 ASN C 236      12.032  72.114  27.109  1.00 16.07           N
ATOM   5493  N   LEU C 237      16.010  70.282  30.134  1.00  6.60           N
ATOM   5494  CA  LEU C 237      16.122  68.825  30.270  1.00  7.41           C
ATOM   5495  C   LEU C 237      16.481  68.430  31.697  1.00  6.01           C
ATOM   5496  O   LEU C 237      15.944  67.448  32.224  1.00  6.47           O
ATOM   5497  CB  LEU C 237      17.151  68.239  29.297  1.00  8.10           C
ATOM   5498  CG  LEU C 237      17.384  66.726  29.347  1.00 10.94           C
ATOM   5499  CD1 LEU C 237      16.055  65.956  29.107  1.00 13.10           C
ATOM   5500  CD2 LEU C 237      18.455  66.271  28.343  1.00 11.63           C
"""

reference_raw_records_alt_seq = """\
CRYST1   40.688   71.918   93.213 108.16  93.25  90.40 P 1           4
ATOM   5485  N   ASN B 246      16.417  72.834  29.095  1.00  7.17           N
ATOM   5486  CA  ASN B 246      15.051  72.312  29.173  1.00  7.74           C
ATOM   5487  C   ASN B 246      15.000  70.818  29.431  1.00  7.38           C
ATOM   5488  O   ASN B 246      14.047  70.141  29.024  1.00  7.80           O
ATOM   5489  CB  ASN B 246      14.281  72.645  27.887  1.00  8.78           C
ATOM   5490  CG  ASN B 246      12.769  72.657  28.088  1.00 13.44           C
ATOM   5491  OD1 ASN B 246      12.265  73.196  29.082  1.00 20.19           O
ATOM   5492  ND2 ASN B 246      12.032  72.114  27.109  1.00 16.07           N
ATOM   5493  N   LEU B 247      16.010  70.282  30.134  1.00  6.60           N
ATOM   5494  CA  LEU B 247      16.122  68.825  30.270  1.00  7.41           C
ATOM   5495  C   LEU B 247      16.481  68.430  31.697  1.00  6.01           C
ATOM   5496  O   LEU B 247      15.944  67.448  32.224  1.00  6.47           O
ATOM   5497  CB  LEU B 247      17.151  68.239  29.297  1.00  8.10           C
ATOM   5498  CG  LEU B 247      17.384  66.726  29.347  1.00 10.94           C
ATOM   5499  CD1 LEU B 247      16.055  65.956  29.107  1.00 13.10           C
ATOM   5500  CD2 LEU B 247      18.455  66.271  28.343  1.00 11.63           C
"""

reference_raw_records_match = """\
CRYST1   40.688   71.918   93.213 108.16  93.25  90.40 P 1           4
ATOM   5485  N   ASN C 270      16.417  72.834  29.095  1.00  7.17           N
ATOM   5486  CA  ASN C 270      15.051  72.312  29.173  1.00  7.74           C
ATOM   5487  C   ASN C 270      15.000  70.818  29.431  1.00  7.38           C
ATOM   5488  O   ASN C 270      14.047  70.141  29.024  1.00  7.80           O
ATOM   5489  CB  ASN C 270      14.281  72.645  27.887  1.00  8.78           C
ATOM   5490  CG  ASN C 270      12.769  72.657  28.088  1.00 13.44           C
ATOM   5491  OD1 ASN C 270      12.265  73.196  29.082  1.00 20.19           O
ATOM   5492  ND2 ASN C 270      12.032  72.114  27.109  1.00 16.07           N
ATOM   5493  N   ALA C 271      16.010  70.282  30.134  1.00  6.60           N
ATOM   5494  CA  ALA C 271      16.122  68.825  30.270  1.00  7.41           C
ATOM   5495  C   ALA C 271      16.481  68.430  31.697  1.00  6.01           C
ATOM   5496  O   ALA C 271      15.944  67.448  32.224  1.00  6.47           O
ATOM   5497  CB  ALA C 271      17.151  68.239  29.297  1.00  8.10           C
"""


def exercise_reference_model(args, mon_lib_srv, ener_lib):
  log = StringIO()
  work_params = reference_model_params.extract()
  work_params.reference_model.enabled = True
  work_params.reference_model.fix_outliers = False
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=flex.split_lines(model_raw_records),
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()
  reference_hierarchy_list = []
  tmp_hierarchy = iotbx.pdb.input(
    source_info=None,
    lines=flex.split_lines(reference_raw_records)).construct_hierarchy()

  reference_hierarchy_list.append(tmp_hierarchy)
  rm = reference_model(
         model=model,
         reference_hierarchy_list=reference_hierarchy_list,
         params=work_params.reference_model,
         log=log)
  assert rm.get_n_proxies() == 5, "Got %d, expected 5" % rm.get_n_proxies()

  reference_hierarchy_list_alt_seq = []
  tmp_hierarchy = iotbx.pdb.input(
    source_info=None,
    lines=flex.split_lines(reference_raw_records_alt_seq)).\
            construct_hierarchy()
  reference_hierarchy_list_alt_seq.append(tmp_hierarchy)
  reference_hierarchy_list_ref_match = []
  tmp_hierarchy = iotbx.pdb.input(
    source_info=None,
    lines=flex.split_lines(reference_raw_records_match)).\
            construct_hierarchy()
  reference_hierarchy_list_ref_match.append(tmp_hierarchy)

  i_seq_name_hash = utils.build_name_hash(
    pdb_hierarchy=pdb_h)
  assert i_seq_name_hash == \
    {0: ' N   ASN C 236     ', 1: ' CA  ASN C 236     ',
     2: ' C   ASN C 236     ', 3: ' O   ASN C 236     ',
     4: ' CB  ASN C 236     ', 5: ' CG  ASN C 236     ',
     6: ' OD1 ASN C 236     ', 7: ' ND2 ASN C 236     ',
     8: ' N   LEU C 237     ', 9: ' CA  LEU C 237     ',
     10: ' C   LEU C 237     ', 11: ' O   LEU C 237     ',
     12: ' CB  LEU C 237     ', 13: ' CG  LEU C 237     ',
     14: ' CD1 LEU C 237     ', 15: ' CD2 LEU C 237     '}
  i_seq_element_hash = utils.build_element_hash(
    pdb_hierarchy=pdb_h)
  assert i_seq_element_hash == \
       {0: 'N', 1: 'C', 2: 'C', 3: 'O', 4: 'C', 5: 'C', 6: 'O', 7: 'N', 8: 'N',
        9: 'C', 10: 'C', 11: 'O', 12: 'C', 13: 'C', 14: 'C', 15: 'C'}

  ref_pdb_hierarchy = reference_hierarchy_list[0]
  dihedral_proxies = \
    utils.get_complete_dihedral_proxies(pdb_hierarchy=ref_pdb_hierarchy)
  sites_cart_ref = ref_pdb_hierarchy.atoms().extract_xyz()
  dihedral_hash = rm.build_dihedral_hash(
    dihedral_proxies=dihedral_proxies,
    sites_cart=sites_cart_ref,
    pdb_hierarchy=ref_pdb_hierarchy,
    include_hydrogens=False,
    include_main_chain=True,
    include_side_chain=True)
  assert len(dihedral_hash) == 5
  reference_dihedral_proxies = rm.reference_dihedral_proxies.deep_copy()
  assert reference_dihedral_proxies is not None
  assert len(reference_dihedral_proxies) == len(dihedral_hash)
  for rdp in reference_dihedral_proxies:
    assert rdp.limit == work_params.reference_model.limit

  r1 = rotalyze(pdb_hierarchy=pdb_h, outliers_only=False)
  out1 = StringIO()
  r1.show_old_output(out=out1)
  r2 = rotalyze(pdb_hierarchy=ref_pdb_hierarchy, outliers_only=False)
  out2 = StringIO()
  r2.show_old_output(out=out2)

  assert not show_diff(out1.getvalue(), """\
 C 236  ASN:1.00:0.2:227.3:80.2:::OUTLIER:OUTLIER
 C 237  LEU:1.00:0.0:209.6:357.2:::OUTLIER:OUTLIER
""")

  assert not show_diff(out2.getvalue(), """\
 C 236  ASN:1.00:39.1:203.2:43.6:::Favored:t0
 C 237  LEU:1.00:60.8:179.1:57.3:::Favored:tp
""")

  xray_structure = pdb_h.extract_xray_structure()
  rm.set_rotamer_to_reference(
    xray_structure=xray_structure,
    mon_lib_srv=mon_lib_srv,
    quiet=True)
  pdb_h.adopt_xray_structure(xray_structure)
  r2 = rotalyze(pdb_hierarchy=pdb_h, outliers_only=False)
  out3 = StringIO()
  r2.show_old_output(out=out3)
  assert not show_diff(out3.getvalue(), """\
 C 236  ASN:1.00:39.1:203.2:43.6:::Favored:t0
 C 237  LEU:1.00:60.8:179.1:57.3:::Favored:tp
""")

  match_map = rm.match_map['ref0']
  assert match_map == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11,
   12: 12, 13: 13, 14: 14, 15: 15}

  master_phil_str_overrides = """
  reference_model.reference_group {
    reference= chain B and resseq 246:247
    selection= chain C and resid 236:237
  }
  """
  def_pars = reference_model_params
  pars = iotbx.phil.parse(master_phil_str_overrides)
  all_pars = def_pars.fetch(pars).extract()
  all_pars.reference_model.enabled = True
  rm = reference_model(
         model = model,
         reference_hierarchy_list=reference_hierarchy_list_alt_seq,
         params=all_pars.reference_model,
         log=log)
  match_map = rm.match_map
  assert match_map['ref0'] == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11,
   12: 12, 13: 13, 14: 14, 15: 15}

  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(file_name=pdb_file,
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()
  # pdb_hierarchy = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
  reference_file_list = []
  reference_file_list.append(pdb_file)

  work_pars = reference_model_params.extract()
  work_pars.reference_model.fix_outliers = False
  work_pars.reference_model.enabled = True

  rm = reference_model(
         model=model,
         reference_file_list=reference_file_list,
         params=work_pars.reference_model,
         log=log)
  reference_dihedral_proxies = rm.reference_dihedral_proxies
  standard_weight = 0
  for dp in reference_dihedral_proxies:
    if dp.weight == 1.0:
      standard_weight += 1
  assert standard_weight == 1181, "Expecting 1181, got %d" % standard_weight
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipping KSDSSP tests: ksdssp module not available.")
  else:
    work_pars = reference_model_params.extract()
    work_pars.reference_model.secondary_structure_only = True
    work_pars.reference_model.enabled = True
    rm.params = work_pars.reference_model
    rm.get_reference_dihedral_proxies(model=model)
    reference_dihedral_proxies = rm.reference_dihedral_proxies
    ss_weight = 0
    for dp in reference_dihedral_proxies:
      if dp.weight == 1.0:
        ss_weight += 1
    assert ss_weight == 694, "expecting 694 proxies, got %d" % ss_weight

  #test SSM alignment
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/ncs/rnase-s.pdb",
    test=os.path.isfile)
  pdb_hierarchy = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
  reference_file_list = []
  reference_file_list.append(pdb_file)
  pdb_hierarchy.reset_i_seq_if_necessary()

  import ccp4io_adaptbx
  ssm = ccp4io_adaptbx.SecondaryStructureMatching(
    reference=pdb_hierarchy.models()[0].chains()[0],
    moving=pdb_hierarchy.models()[0].chains()[1])
  alignment = ccp4io_adaptbx.SSMAlignment.residue_groups(match=ssm)
  assert ssm.GetQvalues()[0] > 0.98

def exercise_multiple_to_one(args, mon_lib_srv, ener_lib):

  pdb_str_original = """\
CRYST1   69.211   49.956   52.557  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   3      51.193  44.956  23.993  1.00 80.52           N
ATOM      2  CA  THR A   3      50.812  43.732  23.211  1.00 80.52           C
ATOM      4  CB  THR A   3      50.446  42.559  24.181  1.00 79.62           C
ATOM      6  OG1 THR A   3      50.206  41.358  23.433  1.00 79.62           O
ATOM      8  CG2 THR A   3      49.239  42.888  25.066  1.00 79.62           C
ATOM     12  C   THR A   3      49.657  44.014  22.221  1.00 80.52           C
ATOM     13  O   THR A   3      48.520  44.223  22.631  1.00 80.52           O
ATOM     17  N   GLY A   4      49.963  44.013  20.917  1.00 79.31           N
ATOM     18  CA  GLY A   4      49.030  44.458  19.892  1.00 79.31           C
ATOM     21  C   GLY A   4      48.761  43.480  18.761  1.00 79.31           C
ATOM     22  O   GLY A   4      47.790  42.725  18.808  1.00 79.31           O
ATOM     24  N   ALA A   5      49.581  43.499  17.715  1.00 78.81           N
ATOM     25  CA  ALA A   5      49.395  42.604  16.581  1.00 78.81           C
ATOM     27  CB  ALA A   5      49.774  43.314  15.283  1.00 77.40           C
ATOM     31  C   ALA A   5      50.195  41.315  16.714  1.00 78.81           C
ATOM     32  O   ALA A   5      50.258  40.537  15.757  1.00 78.81           O
ATOM     34  N   GLN A   6      50.816  41.073  17.872  1.00 80.55           N
ATOM     35  CA  GLN A   6      51.642  39.880  18.018  1.00 80.55           C
ATOM     37  CB  GLN A   6      52.383  39.879  19.354  1.00 79.84           C
ATOM     40  CG  GLN A   6      53.264  41.072  19.596  1.00 79.84           C
ATOM     43  CD  GLN A   6      52.490  42.211  20.225  1.00 79.84           C
ATOM     44  OE1 GLN A   6      51.290  42.091  20.489  1.00 79.84           O
ATOM     45  NE2 GLN A   6      53.167  43.325  20.468  1.00 79.84           N
ATOM     48  C   GLN A   6      50.788  38.631  17.945  1.00 80.55           C
ATOM     49  O   GLN A   6      51.148  37.659  17.273  1.00 80.55           O
ATOM     51  N   VAL A   7      49.643  38.651  18.631  1.00 79.06           N
ATOM     52  CA  VAL A   7      48.822  37.460  18.795  1.00 79.06           C
ATOM     54  CB  VAL A   7      47.610  37.794  19.688  1.00 78.99           C
ATOM     56  CG1 VAL A   7      46.649  36.606  19.794  1.00 78.99           C
ATOM     60  CG2 VAL A   7      48.075  38.245  21.063  1.00 78.99           C
ATOM     64  C   VAL A   7      48.399  36.907  17.450  1.00 79.06           C
ATOM     65  O   VAL A   7      47.962  35.755  17.360  1.00 79.06           O
ATOM     67  N   TYR A   8      48.538  37.700  16.390  1.00 79.78           N
ATOM     68  CA  TYR A   8      48.445  37.147  15.051  1.00 79.78           C
ATOM     70  CB  TYR A   8      48.732  38.228  14.014  1.00 77.69           C
ATOM     73  CG  TYR A   8      48.634  37.736  12.583  1.00 77.69           C
ATOM     74  CD1 TYR A   8      47.404  37.638  11.944  1.00 77.69           C
ATOM     76  CE1 TYR A   8      47.308  37.187  10.640  1.00 77.69           C
ATOM     78  CZ  TYR A   8      48.444  36.802   9.966  1.00 77.69           C
ATOM     79  OH  TYR A   8      48.355  36.348   8.672  1.00 77.69           O
ATOM     81  CE2 TYR A   8      49.672  36.872  10.580  1.00 77.69           C
ATOM     83  CD2 TYR A   8      49.763  37.333  11.883  1.00 77.69           C
ATOM     85  C   TYR A   8      49.416  35.991  14.857  1.00 79.78           C
ATOM     86  O   TYR A   8      49.202  35.164  13.967  1.00 79.78           O
ATOM     88  N   ALA A   9      50.475  35.912  15.671  1.00 79.03           N
ATOM     89  CA  ALA A   9      51.463  34.844  15.546  1.00 79.02           C
ATOM     91  CB  ALA A   9      52.444  34.896  16.719  1.00 79.18           C
ATOM     95  C   ALA A   9      50.833  33.459  15.484  1.00 79.02           C
ATOM     96  O   ALA A   9      51.470  32.524  14.982  1.00 79.02           O
ATOM     98  N   ASN A  10      49.611  33.298  16.002  1.00 79.63           N
ATOM     99  CA  ASN A  10      48.890  32.036  15.896  1.00 79.63           C
ATOM    101  CB  ASN A  10      47.838  31.938  17.002  1.00 78.91           C
ATOM    104  CG  ASN A  10      48.455  31.885  18.387  1.00 78.91           C
ATOM    105  OD1 ASN A  10      49.636  31.603  18.527  1.00 78.91           O
ATOM    106  ND2 ASN A  10      47.648  32.113  19.418  1.00 78.91           N
ATOM    109  C   ASN A  10      48.213  31.859  14.543  1.00 79.63           C
ATOM    110  O   ASN A  10      47.724  30.767  14.246  1.00 79.63           O
TER      58      ASN A  10
ATOM   1990  N   THR B   3      21.107   5.000  45.226  1.00 82.71           N
ATOM   1991  CA  THR B   3      21.187   6.487  45.312  1.00 82.71           C
ATOM   1993  CB  THR B   3      20.105   7.035  46.286  1.00 80.11           C
ATOM   1995  OG1 THR B   3      20.201   6.377  47.557  1.00 80.11           O
ATOM   1997  CG2 THR B   3      18.701   6.831  45.702  1.00 80.11           C
ATOM   2001  C   THR B   3      22.604   6.951  45.721  1.00 82.71           C
ATOM   2002  O   THR B   3      23.561   6.189  45.599  1.00 82.71           O
ATOM   2006  N   GLY B   4      22.752   8.203  46.153  1.00 80.69           N
ATOM   2007  CA  GLY B   4      24.064   8.716  46.532  1.00 80.69           C
ATOM   2010  C   GLY B   4      25.028   8.902  45.376  1.00 80.69           C
ATOM   2011  O   GLY B   4      26.250   8.861  45.572  1.00 80.69           O
ATOM   2013  N   ALA B   5      24.503   9.142  44.177  1.00 80.08           N
ATOM   2014  CA  ALA B   5      25.268   9.118  42.937  1.00 80.09           C
ATOM   2016  CB  ALA B   5      26.031   7.798  42.787  1.00 77.84           C
ATOM   2020  C   ALA B   5      24.301   9.316  41.777  1.00 80.09           C
ATOM   2021  O   ALA B   5      24.660   9.874  40.734  1.00 80.09           O
ATOM   2023  N   GLN B   6      23.035   8.849  42.004  1.00 81.52           N
ATOM   2024  CA  GLN B   6      21.978   8.970  41.003  1.00 81.53           C
ATOM   2026  CB  GLN B   6      20.722   8.250  41.506  1.00 84.25           C
ATOM   2029  CG  GLN B   6      19.920   9.000  42.596  1.00 84.24           C
ATOM   2032  CD  GLN B   6      20.032  10.516  42.500  1.00 84.25           C
ATOM   2033  OE1 GLN B   6      19.770  11.098  41.447  1.00 84.24           O
ATOM   2034  NE2 GLN B   6      20.441  11.159  43.593  1.00 84.25           N
ATOM   2037  C   GLN B   6      21.660  10.426  40.679  1.00 81.52           C
ATOM   2038  O   GLN B   6      21.344  10.750  39.530  1.00 81.52           O
ATOM   2040  N   VAL B   7      21.740  11.307  41.646  1.00 80.27           N
ATOM   2041  CA  VAL B   7      21.376  12.702  41.416  1.00 80.28           C
ATOM   2043  CB  VAL B   7      21.371  13.503  42.738  1.00 79.22           C
ATOM   2045  CG1 VAL B   7      21.092  15.002  42.494  1.00 79.23           C
ATOM   2049  CG2 VAL B   7      20.346  12.946  43.687  1.00 79.22           C
ATOM   2053  C   VAL B   7      22.311  13.348  40.415  1.00 80.27           C
ATOM   2054  O   VAL B   7      21.937  14.328  39.759  1.00 80.27           O
ATOM   2056  N   TYR B   8      23.517  12.809  40.259  1.00 79.95           N
ATOM   2057  CA  TYR B   8      24.474  13.363  39.313  1.00 79.95           C
ATOM   2059  CB  TYR B   8      25.847  12.697  39.486  1.00 79.66           C
ATOM   2062  CG  TYR B   8      26.909  13.218  38.529  1.00 79.66           C
ATOM   2063  CD1 TYR B   8      27.478  14.478  38.703  1.00 79.66           C
ATOM   2065  CE1 TYR B   8      28.444  14.958  37.831  1.00 79.66           C
ATOM   2067  CZ  TYR B   8      28.865  14.173  36.779  1.00 79.66           C
ATOM   2068  OH  TYR B   8      29.825  14.640  35.913  1.00 79.66           O
ATOM   2070  CE2 TYR B   8      28.325  12.919  36.585  1.00 79.66           C
ATOM   2072  CD2 TYR B   8      27.353  12.445  37.459  1.00 79.66           C
ATOM   2074  C   TYR B   8      23.951  13.207  37.884  1.00 79.95           C
ATOM   2075  O   TYR B   8      24.569  13.705  36.937  1.00 79.95           O
ATOM   2077  N   ALA B   9      22.809  12.526  37.712  1.00 80.47           N
ATOM   2078  CA  ALA B   9      22.221  12.378  36.382  1.00 80.47           C
ATOM   2080  CB  ALA B   9      21.051  11.395  36.420  1.00 78.95           C
ATOM   2084  C   ALA B   9      21.758  13.717  35.823  1.00 80.47           C
ATOM   2085  O   ALA B   9      21.827  13.949  34.609  1.00 80.47           O
ATOM   2087  N   ASN B  10      21.261  14.606  36.684  1.00 78.19           N
ATOM   2088  CA  ASN B  10      20.912  15.948  36.235  1.00 78.18           C
ATOM   2090  CB  ASN B  10      20.105  16.644  37.329  1.00 78.39           C
ATOM   2093  CG  ASN B  10      18.743  16.000  37.542  1.00 78.39           C
ATOM   2094  OD1 ASN B  10      18.177  15.401  36.628  1.00 78.39           O
ATOM   2095  ND2 ASN B  10      18.229  16.094  38.762  1.00 78.39           N
ATOM   2098  C   ASN B  10      22.147  16.764  35.859  1.00 78.19           C
ATOM   2099  O   ASN B  10      22.037  17.714  35.076  1.00 78.18           O
TER     116      ASN B  10
ATOM   3968  N   THR C   3      12.127   9.313  24.749  1.00 79.35           N
ATOM   3969  CA  THR C   3      10.942   8.737  24.046  1.00 79.35           C
ATOM   3971  CB  THR C   3      11.262   7.332  23.448  1.00 79.78           C
ATOM   3973  OG1 THR C   3      11.663   6.434  24.490  1.00 79.78           O
ATOM   3975  CG2 THR C   3      12.389   7.415  22.416  1.00 79.78           C
ATOM   3979  C   THR C   3       9.763   8.654  25.028  1.00 79.35           C
ATOM   3980  O   THR C   3       9.889   8.068  26.102  1.00 79.35           O
ATOM   3984  N   GLY C   4       8.622   9.230  24.652  1.00 79.79           N
ATOM   3985  CA  GLY C   4       7.500   9.425  25.558  1.00 79.79           C
ATOM   3988  C   GLY C   4       7.491  10.798  26.210  1.00 79.79           C
ATOM   3989  O   GLY C   4       8.275  11.700  25.895  1.00 79.79           O
ATOM   3991  N   ALA C   5       6.558  10.952  27.145  1.00 80.31           N
ATOM   3992  CA  ALA C   5       6.415  12.204  27.871  1.00 80.31           C
ATOM   3994  CB  ALA C   5       5.000  12.293  28.444  1.00 76.75           C
ATOM   3998  C   ALA C   5       7.437  12.387  28.995  1.00 80.31           C
ATOM   3999  O   ALA C   5       7.578  13.512  29.487  1.00 80.31           O
ATOM   4001  N   GLN C   6       8.160  11.330  29.402  1.00 79.57           N
ATOM   4002  CA  GLN C   6       8.950  11.380  30.637  1.00 79.57           C
ATOM   4004  CB  GLN C   6       9.511   9.998  31.005  1.00 81.75           C
ATOM   4007  CG  GLN C   6      10.700   9.510  30.181  1.00 81.75           C
ATOM   4010  CD  GLN C   6      10.314   9.110  28.788  1.00 81.75           C
ATOM   4011  OE1 GLN C   6       9.147   9.195  28.407  1.00 81.75           O
ATOM   4012  NE2 GLN C   6      11.293   8.668  28.008  1.00 81.75           N
ATOM   4015  C   GLN C   6      10.108  12.363  30.557  1.00 79.57           C
ATOM   4016  O   GLN C   6      10.641  12.760  31.599  1.00 79.57           O
ATOM   4018  N   VAL C   7      10.531  12.731  29.349  1.00 79.61           N
ATOM   4019  CA  VAL C   7      11.538  13.775  29.192  1.00 79.61           C
ATOM   4021  CB  VAL C   7      11.695  14.094  27.694  1.00 78.73           C
ATOM   4023  CG1 VAL C   7      10.350  14.552  27.106  1.00 79.29           C
ATOM   4027  CG2 VAL C   7      12.788  15.133  27.480  1.00 80.51           C
ATOM   4031  C   VAL C   7      11.178  15.010  30.006  1.00 79.62           C
ATOM   4032  O   VAL C   7      12.062  15.759  30.443  1.00 79.61           O
ATOM   4034  N   TYR C   8       9.882  15.234  30.243  1.00 78.80           N
ATOM   4035  CA  TYR C   8       9.422  16.321  31.101  1.00 78.80           C
ATOM   4037  CB  TYR C   8       7.887  16.311  31.116  1.00 79.05           C
ATOM   4040  CG  TYR C   8       7.242  17.382  31.967  1.00 79.05           C
ATOM   4041  CD1 TYR C   8       7.143  18.691  31.510  1.00 79.05           C
ATOM   4043  CE1 TYR C   8       6.548  19.676  32.279  1.00 79.05           C
ATOM   4045  CZ  TYR C   8       6.045  19.358  33.521  1.00 79.05           C
ATOM   4046  OH  TYR C   8       5.457  20.342  34.283  1.00 79.05           O
ATOM   4048  CE2 TYR C   8       6.125  18.064  33.998  1.00 79.05           C
ATOM   4050  CD2 TYR C   8       6.720  17.084  33.219  1.00 79.05           C
ATOM   4052  C   TYR C   8       9.983  16.231  32.521  1.00 78.80           C
ATOM   4053  O   TYR C   8       9.801  17.170  33.302  1.00 78.80           O
ATOM   4055  N   ALA C   9      10.675  15.139  32.866  1.00 79.52           N
ATOM   4056  CA  ALA C   9      11.171  14.948  34.228  1.00 79.52           C
ATOM   4058  CB  ALA C   9      12.014  13.674  34.293  1.00 78.34           C
ATOM   4062  C   ALA C   9      11.983  16.145  34.702  1.00 79.52           C
ATOM   4063  O   ALA C   9      11.793  16.641  35.818  1.00 79.52           O
ATOM   4065  N   ASN C  10      12.896  16.627  33.865  1.00 80.25           N
ATOM   4066  CA  ASN C  10      13.672  17.797  34.239  1.00 80.25           C
ATOM   4068  CB  ASN C  10      14.712  18.063  33.172  1.00 78.17           C
ATOM   4071  CG  ASN C  10      15.782  17.007  33.161  1.00 78.17           C
ATOM   4072  OD1 ASN C  10      16.004  16.325  34.166  1.00 78.17           O
ATOM   4073  ND2 ASN C  10      16.442  16.845  32.028  1.00 78.17           N
ATOM   4076  C   ASN C  10      12.798  19.015  34.457  1.00 80.25           C
ATOM   4077  O   ASN C  10      13.290  20.040  34.941  1.00 80.25           O
TER     174      ASN C  10
ATOM   5959  N   THR D   3      60.805  23.774   6.731  1.00 77.43           N
ATOM   5960  CA  THR D   3      61.763  22.725   7.191  1.00 77.43           C
ATOM   5962  CB  THR D   3      62.603  22.175   6.010  1.00 78.92           C
ATOM   5964  OG1 THR D   3      63.305  23.243   5.360  1.00 78.92           O
ATOM   5966  CG2 THR D   3      61.703  21.469   5.000  1.00 78.92           C
ATOM   5970  C   THR D   3      62.675  23.293   8.284  1.00 77.43           C
ATOM   5971  O   THR D   3      62.761  24.506   8.443  1.00 77.43           O
ATOM   5975  N   GLY D   4      63.363  22.412   9.022  1.00 79.20           N
ATOM   5976  CA  GLY D   4      64.130  22.797  10.196  1.00 79.20           C
ATOM   5979  C   GLY D   4      63.309  22.788  11.472  1.00 79.20           C
ATOM   5980  O   GLY D   4      62.145  22.393  11.509  1.00 79.20           O
ATOM   5982  N   ALA D   5      63.950  23.233  12.557  1.00 80.19           N
ATOM   5983  CA  ALA D   5      63.257  23.361  13.836  1.00 80.19           C
ATOM   5985  CB  ALA D   5      64.211  23.993  14.857  1.00 75.84           C
ATOM   5989  C   ALA D   5      61.970  24.181  13.714  1.00 80.19           C
ATOM   5990  O   ALA D   5      60.999  23.931  14.438  1.00 80.19           O
ATOM   5992  N   GLN D   6      61.942  25.142  12.784  1.00 78.97           N
ATOM   5993  CA  GLN D   6      60.843  26.092  12.591  1.00 78.97           C
ATOM   5995  CB  GLN D   6      61.204  27.062  11.469  1.00 80.12           C
ATOM   5998  CG  GLN D   6      61.464  26.355  10.144  1.00 80.12           C
ATOM   6001  CD  GLN D   6      61.853  27.306   9.032  1.00 80.12           C
ATOM   6002  OE1 GLN D   6      62.179  28.464   9.288  1.00 80.12           O
ATOM   6003  NE2 GLN D   6      61.851  26.812   7.790  1.00 80.12           N
ATOM   6006  C   GLN D   6      59.510  25.447  12.245  1.00 78.96           C
ATOM   6007  O   GLN D   6      58.509  26.166  12.139  1.00 78.96           O
ATOM   6009  N   VAL D   7      59.474  24.140  11.995  1.00 78.86           N
ATOM   6010  CA  VAL D   7      58.194  23.449  11.865  1.00 78.86           C
ATOM   6012  CB  VAL D   7      58.425  21.993  11.421  1.00 81.21           C
ATOM   6014  CG1 VAL D   7      58.877  21.975   9.986  1.00 81.21           C
ATOM   6018  CG2 VAL D   7      59.474  21.288  12.321  1.00 81.21           C
ATOM   6022  C   VAL D   7      57.423  23.523  13.168  1.00 78.86           C
ATOM   6023  O   VAL D   7      56.190  23.411  13.186  1.00 78.86           O
ATOM   6025  N   TYR D   8      58.138  23.697  14.277  1.00 79.34           N
ATOM   6026  CA  TYR D   8      57.515  23.918  15.568  1.00 79.34           C
ATOM   6028  CB  TYR D   8      58.584  23.823  16.649  1.00 79.01           C
ATOM   6031  CG  TYR D   8      58.096  24.160  18.024  1.00 79.01           C
ATOM   6032  CD1 TYR D   8      57.220  23.317  18.688  1.00 79.01           C
ATOM   6034  CE1 TYR D   8      56.778  23.608  19.947  1.00 79.01           C
ATOM   6036  CZ  TYR D   8      57.227  24.739  20.578  1.00 79.01           C
ATOM   6037  OH  TYR D   8      56.779  25.015  21.845  1.00 79.01           O
ATOM   6039  CE2 TYR D   8      58.111  25.590  19.948  1.00 79.01           C
ATOM   6041  CD2 TYR D   8      58.544  25.294  18.680  1.00 79.01           C
ATOM   6043  C   TYR D   8      56.807  25.258  15.636  1.00 79.34           C
ATOM   6044  O   TYR D   8      55.950  25.447  16.505  1.00 79.34           O
ATOM   6046  N   ALA D   9      57.137  26.174  14.730  1.00 78.81           N
ATOM   6047  CA  ALA D   9      56.591  27.522  14.741  1.00 78.81           C
ATOM   6049  CB  ALA D   9      56.758  28.183  13.374  1.00 79.37           C
ATOM   6053  C   ALA D   9      55.127  27.498  15.121  1.00 78.81           C
ATOM   6054  O   ALA D   9      54.764  27.896  16.226  1.00 78.81           O
ATOM   6056  N   ASN D  10      54.284  26.983  14.233  1.00 80.12           N
ATOM   6057  CA  ASN D  10      52.848  27.017  14.467  1.00 80.13           C
ATOM   6059  CB  ASN D  10      52.140  26.402  13.274  1.00 80.26           C
ATOM   6062  CG  ASN D  10      52.645  25.031  12.969  1.00 80.26           C
ATOM   6063  OD1 ASN D  10      53.101  24.311  13.860  1.00 80.26           O
ATOM   6064  ND2 ASN D  10      52.586  24.655  11.705  1.00 80.26           N
ATOM   6067  C   ASN D  10      52.420  26.314  15.753  1.00 80.14           C
ATOM   6068  O   ASN D  10      51.225  26.333  16.068  1.00 80.16           O
TER     232      ASN D  10
END
"""

  pdb_str_ref_minimized = """\
CRYST1   69.211   49.956   52.557  90.00  90.00  90.00 P 1
SCALE1      0.014449  0.000000  0.000000        0.00000
SCALE2      0.000000  0.020018  0.000000        0.00000
SCALE3      0.000000  0.000000  0.019027        0.00000
ATOM      1  N   THR A   3      50.767  43.905  24.734  1.00 80.52           N
ATOM      2  CA  THR A   3      50.582  43.115  23.523  1.00 80.52           C
ATOM      4  CB  THR A   3      49.583  41.964  23.746  1.00 79.62           C
ATOM      6  OG1 THR A   3      49.442  41.209  22.536  1.00 79.62           O
ATOM      8  CG2 THR A   3      48.225  42.510  24.160  1.00 79.62           C
ATOM     12  C   THR A   3      50.093  43.985  22.370  1.00 80.52           C
ATOM     13  O   THR A   3      49.756  45.154  22.562  1.00 80.52           O
ATOM     17  N   GLY A   4      50.055  43.408  21.174  1.00 79.31           N
ATOM     18  CA  GLY A   4      49.609  44.126  19.994  1.00 79.31           C
ATOM     21  C   GLY A   4      49.582  43.269  18.744  1.00 79.31           C
ATOM     22  O   GLY A   4      48.530  43.075  18.136  1.00 79.31           O
ATOM     24  N   ALA A   5      50.746  42.754  18.361  1.00 78.81           N
ATOM     25  CA  ALA A   5      50.858  41.913  17.175  1.00 78.81           C
ATOM     27  CB  ALA A   5      51.904  42.473  16.224  1.00 77.40           C
ATOM     31  C   ALA A   5      51.199  40.476  17.556  1.00 78.81           C
ATOM     32  O   ALA A   5      51.983  39.814  16.877  1.00 78.81           O
ATOM     34  N   GLN A   6      50.604  40.001  18.645  1.00 80.55           N
ATOM     35  CA  GLN A   6      50.843  38.643  19.118  1.00 80.55           C
ATOM     37  CB  GLN A   6      51.379  38.655  20.554  1.00 79.84           C
ATOM     40  CG  GLN A   6      52.763  39.273  20.711  1.00 79.84           C
ATOM     43  CD  GLN A   6      52.740  40.791  20.698  1.00 79.84           C
ATOM     44  OE1 GLN A   6      51.676  41.408  20.641  1.00 79.84           O
ATOM     45  NE2 GLN A   6      53.919  41.400  20.750  1.00 79.84           N
ATOM     48  C   GLN A   6      49.570  37.807  19.041  1.00 80.55           C
ATOM     49  O   GLN A   6      49.417  36.823  19.765  1.00 80.55           O
ATOM     51  N   VAL A   7      48.660  38.205  18.158  1.00 79.06           N
ATOM     52  CA  VAL A   7      47.399  37.495  17.985  1.00 79.06           C
ATOM     54  CB  VAL A   7      46.201  38.453  18.098  1.00 78.99           C
ATOM     56  CG1 VAL A   7      44.896  37.694  17.913  1.00 78.99           C
ATOM     60  CG2 VAL A   7      46.222  39.173  19.437  1.00 78.99           C
ATOM     64  C   VAL A   7      47.382  36.766  16.647  1.00 79.06           C
ATOM     65  O   VAL A   7      47.123  35.564  16.586  1.00 79.06           O
ATOM     67  N   TYR A   8      47.661  37.501  15.575  1.00 79.78           N
ATOM     68  CA  TYR A   8      47.677  36.928  14.235  1.00 79.78           C
ATOM     70  CB  TYR A   8      47.601  38.032  13.178  1.00 77.69           C
ATOM     73  CG  TYR A   8      47.544  37.522  11.755  1.00 77.69           C
ATOM     74  CD1 TYR A   8      46.344  37.106  11.194  1.00 77.69           C
ATOM     76  CE1 TYR A   8      46.286  36.641   9.894  1.00 77.69           C
ATOM     78  CZ  TYR A   8      47.437  36.589   9.136  1.00 77.69           C
ATOM     79  OH  TYR A   8      47.384  36.126   7.842  1.00 77.69           O
ATOM     81  CE2 TYR A   8      48.641  36.998   9.670  1.00 77.69           C
ATOM     83  CD2 TYR A   8      48.689  37.462  10.971  1.00 77.69           C
ATOM     85  C   TYR A   8      48.925  36.077  14.022  1.00 79.78           C
ATOM     86  O   TYR A   8      48.903  35.104  13.267  1.00 79.78           O
ATOM     88  N   ALA A   9      50.010  36.450  14.691  1.00 79.03           N
ATOM     89  CA  ALA A   9      51.269  35.722  14.577  1.00 79.02           C
ATOM     91  CB  ALA A   9      52.423  36.578  15.075  1.00 79.18           C
ATOM     95  C   ALA A   9      51.210  34.407  15.346  1.00 79.02           C
ATOM     96  O   ALA A   9      51.871  33.434  14.982  1.00 79.02           O
ATOM     98  N   ASN A  10      50.416  34.384  16.411  1.00 79.63           N
ATOM     99  CA  ASN A  10      50.270  33.188  17.233  1.00 79.63           C
ATOM    101  CB  ASN A  10      50.147  33.564  18.711  1.00 78.91           C
ATOM    104  CG  ASN A  10      51.365  34.304  19.227  1.00 78.91           C
ATOM    105  OD1 ASN A  10      52.472  34.134  18.716  1.00 78.91           O
ATOM    106  ND2 ASN A  10      51.167  35.132  20.246  1.00 78.91           N
ATOM    109  C   ASN A  10      49.059  32.369  16.797  1.00 79.63           C
ATOM    110  O   ASN A  10      49.101  31.139  16.787  1.00 79.63           O
TER
  """

  ref_file = open("ref.pdb", 'w')
  ref_file.write(pdb_str_ref_minimized)
  ref_file.close()
  log = StringIO()
  # log = sys.stdout
  # orig_file = open("start.pdb", "w")
  # orig_file.write(pdb_str_original)
  # orig_file.close()
  def_pars = reference_model_params
  params_text = """\
 reference_model {
    reference_group {
      reference = chain 'A'
      selection = chain 'A'
      file_name = "ref.pdb"
    }
    reference_group {
      reference = chain 'A'
      selection = chain 'B'
      file_name = "ref.pdb"
    }
    reference_group {
      reference = chain 'A'
      selection = chain 'C'
      file_name = "ref.pdb"
    }
    reference_group {
      reference = chain 'A'
      selection = chain 'D'
      file_name = "ref.pdb"
    }
  }  """
  pars = iotbx.phil.parse(params_text)
  all_pars = def_pars.fetch(pars).extract()
  all_pars.reference_model.enabled = True

  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=flex.split_lines(pdb_str_original),
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()
  rm = reference_model(
         model = model,
         reference_file_list=['ref.pdb'],
         params=all_pars.reference_model,
         log=log)
  # rm.show_reference_summary(log=log)
  assert rm.get_n_proxies() == 124, \
      "Expecting 124 proxies, got %d" % rm.get_n_proxies()
  # STOP()
  new_h = pdb_h.deep_copy()
  xray_structure = new_h.extract_xray_structure()
  rm.set_rotamer_to_reference(
    xray_structure=xray_structure)
  new_h.adopt_xray_structure(xray_structure)
  r1 = rotalyze(pdb_hierarchy=new_h, outliers_only=False)
  assert r1.n_outliers == 0
  # new_h.write_pdb_file(file_name="final.pdb")

  #
  # The same, but from multiple files
  for i in range(4):
    ref_file = open("ref_%d.pdb" % i, 'w')
    ref_file.write(pdb_str_ref_minimized)
    ref_file.close()
  def_pars = reference_model_params
  params_text = """\
  reference_model {
    file = ref_0.pdb
    file = ref_1.pdb
    file = ref_2.pdb
    file = ref_3.pdb
    reference_group {
      reference = chain 'A'
      selection = chain 'A'
      file_name = "ref_0.pdb"
    }
    reference_group {
      reference = chain 'A'
      selection = chain 'B'
      file_name = "ref_1.pdb"
    }
    reference_group {
      reference = chain 'A'
      selection = chain 'C'
      file_name = "ref_2.pdb"
    }
    reference_group {
      reference = chain 'A'
      selection = chain 'D'
      file_name = "ref_3.pdb"
    }
  }  """
  pars = iotbx.phil.parse(params_text)
  all_pars = def_pars.fetch(pars).extract()
  all_pars.reference_model.enabled = True
  rm = reference_model(
         model=model,
         reference_file_list=['ref_0.pdb', 'ref_1.pdb', 'ref_2.pdb', 'ref_3.pdb'],
         params=all_pars.reference_model,
         log=log)
  assert rm.get_n_proxies() == 124, \
      "Expecting 124 proxies, got %d" % rm.get_n_proxies()
  for i in range(4):
    os.remove("ref_%d.pdb" % i)
  #
  # The same, 1 group, should be 116/4=29 proxies
  ref_file = open("ref_0.pdb", 'w')
  ref_file.write(pdb_str_ref_minimized)
  ref_file.close()
  def_pars = reference_model_params
  params_text = """\
  reference_model {
    file = ref_0.pdb
    reference_group {
      reference = chain 'A'
      selection = chain 'A'
      file_name = "ref_0.pdb"
    }
  }  """
  pars = iotbx.phil.parse(params_text)
  all_pars = def_pars.fetch(pars).extract()
  all_pars.reference_model.enabled = True
  rm = reference_model(
         model=model,
         reference_file_list=['ref_0.pdb'],
         params=all_pars.reference_model,
         log=log)
  assert rm.get_n_proxies() == 31, \
      "Expecting 31 proxies, got %d" % rm.get_n_proxies()
  all_pars.reference_model.side_chain=False
  rm = reference_model(
         model=model,
         reference_file_list=['ref_0.pdb'],
         params=all_pars.reference_model,
         log=log)
  assert rm.get_n_proxies() == 21, \
      "Expecting 21 proxies, got %d" % rm.get_n_proxies()
  all_pars.reference_model.side_chain=True
  all_pars.reference_model.main_chain=False
  rm = reference_model(
         model=model,
         reference_file_list=['ref_0.pdb'],
         params=all_pars.reference_model,
         log=log)
  assert rm.get_n_proxies() == 10, \
      "Expecting 10 proxies, got %d" % rm.get_n_proxies()

  # just throw all in without specifying:
  all_pars = def_pars.fetch().extract()
  all_pars.reference_model.enabled = True
  all_pars.reference_model.file = 'ref_0.pdb'
  rm = reference_model(
         model=model,
         reference_file_list=['ref_0.pdb'],
         params=all_pars.reference_model,
         log=log)
  assert rm.get_n_proxies() == 124, \
      "Expecting 124 proxies, got %d" % rm.get_n_proxies()
  os.remove("ref_0.pdb")

  # reference on self and make sure it is chains A<->A, B<->B etc
  log = StringIO()
  def_pars = reference_model_params
  all_pars = def_pars.fetch().extract()
  all_pars.reference_model.enabled = True
  all_pars.reference_model.use_starting_model_as_reference = True
  rm = reference_model(
      model=model,
      reference_hierarchy_list=\
          [model.get_hierarchy()],
      params=all_pars.reference_model,
      log=log)
  rm.show_reference_summary(log=log)
  log_strings = log.getvalue().split("\n")
  # print rm.get_n_proxies()
  # print "=========="
  # print "\n".join(log_strings)
  # print "=========="
  assert rm.get_n_proxies() == 124, \
      "Expecting 124 proxies, got %d" % rm.get_n_proxies()
  for needed_string in [
      "GLN A   6  <=====>  GLN A   6",
      "ALA A   9  <=====>  ALA A   9",
      "ASN A  10  <=====>  ASN A  10",
      "THR B   3  <=====>  THR B   3",
      "GLN B   6  <=====>  GLN B   6",
      "ALA B   9  <=====>  ALA B   9",
      "ASN B  10  <=====>  ASN B  10",
      "THR C   3  <=====>  THR C   3",
      "GLN C   6  <=====>  GLN C   6",
      "ALA D   5  <=====>  ALA D   5",
      "GLN D   6  <=====>  GLN D   6",
      ]:
    assert needed_string in log_strings, "'%s' not in log!" % needed_string



def exercise_multiple_ncs_groups_found(mon_lib_srv, ener_lib):
  pdb_str_original = """\
CRYST1   49.945   53.842   33.425  90.00  90.00  90.00 P 1
ATOM   5466  N   ASN C 236       9.580  47.176  25.356  1.00 60.13           N
ATOM   5467  CA  ASN C 236       8.200  46.668  25.442  1.00 60.52           C
ATOM   5468  C   ASN C 236       8.058  45.158  25.655  1.00 60.49           C
ATOM   5469  O   ASN C 236       7.110  44.527  25.191  1.00 60.60           O
ATOM   5470  CB  ASN C 236       7.370  47.129  24.244  1.00 60.55           C
ATOM   5471  CG  ASN C 236       6.038  47.680  24.666  1.00 61.75           C
ATOM   5472  OD1 ASN C 236       5.937  48.842  25.096  1.00 62.86           O
ATOM   5473  ND2 ASN C 236       5.000  46.849  24.575  1.00 61.22           N
ATOM   5474  N   LEU C 237       8.997  44.597  26.396  1.00 60.55           N
ATOM   5475  CA  LEU C 237       9.125  43.164  26.494  1.00 60.48           C
ATOM   5476  C   LEU C 237       9.496  42.788  27.909  1.00 60.06           C
ATOM   5477  O   LEU C 237       9.016  41.780  28.425  1.00 60.41           O
ATOM   5478  CB  LEU C 237      10.199  42.697  25.506  1.00 60.91           C
ATOM   5479  CG  LEU C 237      10.223  41.328  24.819  1.00 62.25           C
ATOM   5480  CD1 LEU C 237       9.088  40.386  25.251  1.00 63.18           C
ATOM   5481  CD2 LEU C 237      10.244  41.542  23.293  1.00 62.89           C
TER
ATOM      1  N   THR A   3      42.874  19.189  20.321  1.00 80.52           N
ATOM      2  CA  THR A   3      42.493  17.965  19.539  1.00 80.52           C
ATOM      4  CB  THR A   3      42.127  16.792  20.509  1.00 79.62           C
ATOM      6  OG1 THR A   3      41.887  15.591  19.761  1.00 79.62           O
ATOM      8  CG2 THR A   3      40.920  17.121  21.394  1.00 79.62           C
ATOM     12  C   THR A   3      41.338  18.247  18.549  1.00 80.52           C
ATOM     13  O   THR A   3      40.201  18.456  18.959  1.00 80.52           O
ATOM     17  N   GLY A   4      41.644  18.246  17.245  1.00 79.31           N
ATOM     18  CA  GLY A   4      40.711  18.691  16.220  1.00 79.31           C
ATOM     21  C   GLY A   4      40.442  17.713  15.089  1.00 79.31           C
ATOM     22  O   GLY A   4      39.471  16.958  15.136  1.00 79.31           O
ATOM     24  N   ALA A   5      41.262  17.732  14.043  1.00 78.81           N
ATOM     25  CA  ALA A   5      41.076  16.837  12.909  1.00 78.81           C
ATOM     27  CB  ALA A   5      41.455  17.547  11.611  1.00 77.40           C
ATOM     31  C   ALA A   5      41.876  15.548  13.042  1.00 78.81           C
ATOM     32  O   ALA A   5      41.939  14.770  12.085  1.00 78.81           O
ATOM     34  N   GLN A   6      42.497  15.306  14.200  1.00 80.55           N
ATOM     35  CA  GLN A   6      43.323  14.113  14.346  1.00 80.55           C
ATOM     37  CB  GLN A   6      44.064  14.112  15.682  1.00 79.84           C
ATOM     40  CG  GLN A   6      44.945  15.305  15.924  1.00 79.84           C
ATOM     43  CD  GLN A   6      44.171  16.444  16.553  1.00 79.84           C
ATOM     44  OE1 GLN A   6      42.971  16.324  16.817  1.00 79.84           O
ATOM     45  NE2 GLN A   6      44.848  17.558  16.796  1.00 79.84           N
ATOM     48  C   GLN A   6      42.469  12.864  14.273  1.00 80.55           C
ATOM     49  O   GLN A   6      42.829  11.892  13.601  1.00 80.55           O
ATOM     51  N   VAL A   7      41.324  12.884  14.959  1.00 79.06           N
ATOM     52  CA  VAL A   7      40.503  11.693  15.123  1.00 79.06           C
ATOM     54  CB  VAL A   7      39.291  12.027  16.016  1.00 78.99           C
ATOM     56  CG1 VAL A   7      38.330  10.839  16.122  1.00 78.99           C
ATOM     60  CG2 VAL A   7      39.756  12.478  17.391  1.00 78.99           C
ATOM     64  C   VAL A   7      40.080  11.140  13.778  1.00 79.06           C
ATOM     65  O   VAL A   7      39.643   9.988  13.688  1.00 79.06           O
ATOM     67  N   TYR A   8      40.219  11.933  12.718  1.00 79.78           N
ATOM     68  CA  TYR A   8      40.126  11.380  11.379  1.00 79.78           C
ATOM     70  CB  TYR A   8      40.413  12.461  10.342  1.00 77.69           C
ATOM     73  CG  TYR A   8      40.315  11.969   8.911  1.00 77.69           C
ATOM     74  CD1 TYR A   8      39.085  11.871   8.272  1.00 77.69           C
ATOM     76  CE1 TYR A   8      38.989  11.420   6.968  1.00 77.69           C
ATOM     78  CZ  TYR A   8      40.125  11.035   6.294  1.00 77.69           C
ATOM     79  OH  TYR A   8      40.036  10.581   5.000  1.00 77.69           O
ATOM     81  CE2 TYR A   8      41.353  11.105   6.908  1.00 77.69           C
ATOM     83  CD2 TYR A   8      41.444  11.566   8.211  1.00 77.69           C
ATOM     85  C   TYR A   8      41.097  10.224  11.185  1.00 79.78           C
ATOM     86  O   TYR A   8      40.883   9.397  10.295  1.00 79.78           O
ATOM     88  N   ALA A   9      42.156  10.145  11.999  1.00 79.03           N
ATOM     89  CA  ALA A   9      43.144   9.077  11.874  1.00 79.02           C
ATOM     91  CB  ALA A   9      44.125   9.129  13.047  1.00 79.18           C
ATOM     95  C   ALA A   9      42.514   7.692  11.812  1.00 79.02           C
ATOM     96  O   ALA A   9      43.151   6.757  11.310  1.00 79.02           O
ATOM     98  N   ASN A  10      41.292   7.531  12.330  1.00 79.63           N
ATOM     99  CA  ASN A  10      40.571   6.269  12.224  1.00 79.63           C
ATOM    101  CB  ASN A  10      39.519   6.171  13.330  1.00 78.91           C
ATOM    104  CG  ASN A  10      40.136   6.118  14.715  1.00 78.91           C
ATOM    105  OD1 ASN A  10      41.317   5.836  14.855  1.00 78.91           O
ATOM    106  ND2 ASN A  10      39.329   6.346  15.746  1.00 78.91           N
ATOM    109  C   ASN A  10      39.894   6.092  10.871  1.00 79.63           C
ATOM    110  O   ASN A  10      39.405   5.000  10.574  1.00 79.63           O
TER
END
  """

  pdb_str_ref = """\
CRYST1   49.945   53.842   33.425  90.00  90.00  90.00 P 1
ATOM   5466  N   ASN C 236      10.328  45.698  25.449  1.00 60.13           N
ATOM   5467  CA  ASN C 236       8.971  45.973  25.787  1.00 60.52           C
ATOM   5468  C   ASN C 236       8.271  44.664  25.724  1.00 60.49           C
ATOM   5469  O   ASN C 236       7.276  44.532  25.017  1.00 60.60           O
ATOM   5470  CB  ASN C 236       8.337  46.962  24.776  1.00 60.55           C
ATOM   5471  CG  ASN C 236       7.235  47.762  25.415  1.00 61.75           C
ATOM   5472  OD1 ASN C 236       6.331  47.222  26.063  1.00 62.86           O
ATOM   5473  ND2 ASN C 236       7.315  49.079  25.302  1.00 61.22           N
ATOM   5474  N   LEU C 237       8.820  43.663  26.441  1.00 60.55           N
ATOM   5475  CA  LEU C 237       8.420  42.305  26.286  1.00 60.48           C
ATOM   5476  C   LEU C 237       8.713  41.508  27.558  1.00 60.06           C
ATOM   5477  O   LEU C 237       7.907  41.421  28.503  1.00 60.41           O
ATOM   5478  CB  LEU C 237       9.159  41.598  25.114  1.00 60.91           C
ATOM   5479  CG  LEU C 237       9.365  42.136  23.662  1.00 62.25           C
ATOM   5480  CD1 LEU C 237      10.605  42.996  23.496  1.00 63.18           C
ATOM   5481  CD2 LEU C 237       9.419  40.966  22.765  1.00 62.89           C
TER
ATOM      1  N   THR A   3      40.527  19.363  20.612  1.00 80.52           N
ATOM      2  CA  THR A   3      41.278  18.625  19.636  1.00 80.52           C
ATOM      4  CB  THR A   3      40.971  17.090  19.710  1.00 79.62           C
ATOM      6  OG1 THR A   3      40.039  16.849  20.760  1.00 79.62           O
ATOM      8  CG2 THR A   3      42.308  16.246  19.999  1.00 79.62           C
ATOM     12  C   THR A   3      40.899  19.134  18.229  1.00 80.52           C
ATOM     13  O   THR A   3      39.780  19.542  17.983  1.00 80.52           O
ATOM     17  N   GLY A   4      41.890  19.246  17.384  1.00 79.31           N
ATOM     18  CA  GLY A   4      41.732  19.850  16.092  1.00 79.31           C
ATOM     21  C   GLY A   4      41.306  18.930  14.985  1.00 79.31           C
ATOM     22  O   GLY A   4      40.121  18.885  14.657  1.00 79.31           O
ATOM     24  N   ALA A   5      42.279  18.233  14.402  1.00 78.81           N
ATOM     25  CA  ALA A   5      41.969  17.264  13.392  1.00 78.81           C
ATOM     27  CB  ALA A   5      42.474  17.741  12.001  1.00 77.40           C
ATOM     31  C   ALA A   5      42.643  15.914  13.751  1.00 78.81           C
ATOM     32  O   ALA A   5      43.503  15.474  12.983  1.00 78.81           O
ATOM     34  N   GLN A   6      42.216  15.310  14.835  1.00 80.55           N
ATOM     35  CA  GLN A   6      42.871  14.115  15.363  1.00 80.55           C
ATOM     37  CB  GLN A   6      43.590  14.383  16.698  1.00 79.84           C
ATOM     40  CG  GLN A   6      44.888  15.121  16.536  1.00 79.84           C
ATOM     43  CD  GLN A   6      44.671  16.613  16.295  1.00 79.84           C
ATOM     44  OE1 GLN A   6      44.164  17.330  17.155  1.00 79.84           O
ATOM     45  NE2 GLN A   6      45.100  17.105  15.149  1.00 79.84           N
ATOM     48  C   GLN A   6      41.888  12.972  15.564  1.00 80.55           C
ATOM     49  O   GLN A   6      42.024  12.228  16.514  1.00 80.55           O
ATOM     51  N   VAL A   7      40.933  12.858  14.656  1.00 79.06           N
ATOM     52  CA  VAL A   7      40.101  11.677  14.619  1.00 79.06           C
ATOM     54  CB  VAL A   7      38.947  11.709  15.573  1.00 78.99           C
ATOM     56  CG1 VAL A   7      39.330  11.128  16.941  1.00 78.99           C
ATOM     60  CG2 VAL A   7      38.334  13.128  15.699  1.00 78.99           C
ATOM     64  C   VAL A   7      39.594  11.421  13.214  1.00 79.06           C
ATOM     65  O   VAL A   7      38.407  11.279  12.954  1.00 79.06           O
ATOM     67  N   TYR A   8      40.568  11.433  12.304  1.00 79.78           N
ATOM     68  CA  TYR A   8      40.360  10.983  10.905  1.00 79.78           C
ATOM     70  CB  TYR A   8      40.783  12.069   9.904  1.00 77.69           C
ATOM     73  CG  TYR A   8      40.349  11.670   8.527  1.00 77.69           C
ATOM     74  CD1 TYR A   8      39.008  11.604   8.192  1.00 77.69           C
ATOM     76  CE1 TYR A   8      38.600  11.184   6.931  1.00 77.69           C
ATOM     78  CZ  TYR A   8      39.528  10.864   5.979  1.00 77.69           C
ATOM     79  OH  TYR A   8      39.195  10.466   4.696  1.00 77.69           O
ATOM     81  CE2 TYR A   8      40.880  10.918   6.304  1.00 77.69           C
ATOM     83  CD2 TYR A   8      41.286  11.303   7.563  1.00 77.69           C
ATOM     85  C   TYR A   8      41.107   9.702  10.624  1.00 79.78           C
ATOM     86  O   TYR A   8      40.892   9.064   9.584  1.00 79.78           O
ATOM     88  N   ALA A   9      42.003   9.312  11.526  1.00 79.03           N
ATOM     89  CA  ALA A   9      42.888   8.166  11.317  1.00 79.02           C
ATOM     91  CB  ALA A   9      44.052   8.246  12.231  1.00 79.18           C
ATOM     95  C   ALA A   9      42.102   6.856  11.504  1.00 79.02           C
ATOM     96  O   ALA A   9      42.154   5.981  10.647  1.00 79.02           O
ATOM     98  N   ASN A  10      41.404   6.751  12.642  1.00 79.63           N
ATOM     99  CA  ASN A  10      40.465   5.684  12.913  1.00 79.63           C
ATOM    101  CB  ASN A  10      39.947   5.766  14.373  1.00 78.91           C
ATOM    104  CG  ASN A  10      41.037   5.501  15.391  1.00 78.91           C
ATOM    105  OD1 ASN A  10      42.073   4.895  15.058  1.00 78.91           O
ATOM    106  ND2 ASN A  10      40.820   5.957  16.635  1.00 78.91           N
ATOM    109  C   ASN A  10      39.283   5.748  11.958  1.00 79.63           C
ATOM    110  O   ASN A  10      39.365   5.382  10.797  1.00 79.63           O
TER
  """
  ref_file = open("ref.pdb", 'w')
  ref_file.write(pdb_str_ref)
  ref_file.close()
  log = StringIO()
  # log = sys.stdout

  def_pars = reference_model_params
  all_pars = def_pars.fetch().extract()
  all_pars.reference_model.file = 'ref.pdb'
  all_pars.reference_model.enabled = True
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=flex.split_lines(pdb_str_original),
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()
  rm = reference_model(
         model=model,
         reference_file_list=['ref.pdb'],
         params=all_pars.reference_model,
         log=log)
  assert rm.get_n_proxies() == 36, \
      "Expecting 36 proxies, got %d" % rm.get_n_proxies()
  os.remove("ref.pdb")


def exercise_cutted_residue(mon_lib_srv, ener_lib):
  pdb_str_original = """\
CRYST1  117.739  195.224  119.094  90.00 101.60  90.00 P 1 21 1
ATOM   6368  N   THR K 332       4.163  72.088  52.141  1.00171.28           N
ATOM   6369  CA  THR K 332       2.830  71.741  52.608  1.00153.71           C
ATOM   6370  C   THR K 332       1.990  70.958  51.609  1.00132.45           C
ATOM   6371  O   THR K 332       2.224  71.000  50.405  1.00130.38           O
ATOM   6372  CB  THR K 332       2.047  72.996  53.035  1.00155.45           C
ATOM   6373  N   VAL K 333       1.006  70.246  52.144  1.00121.58           N
ATOM   6374  CA  VAL K 333       0.085  69.440  51.360  1.00129.11           C
ATOM   6375  C   VAL K 333      -1.326  69.771  51.818  1.00146.57           C
ATOM   6376  O   VAL K 333      -1.517  70.242  52.935  1.00151.92           O
ATOM   6377  CB  VAL K 333       0.342  67.942  51.562  1.00126.37           C
ATOM   6378  N   SER K 334      -2.318  69.535  50.968  1.00156.08           N
ATOM   6379  CA  SER K 334      -3.687  69.866  51.335  1.00158.16           C
ATOM   6380  C   SER K 334      -4.197  69.116  52.555  1.00157.55           C
ATOM   6381  O   SER K 334      -4.066  67.905  52.664  1.00161.93           O
ATOM   6382  CB  SER K 334      -4.630  69.614  50.166  1.00162.09           C
ATOM   6383  OG  SER K 334      -5.836  69.041  50.632  1.00170.98           O
END
  """

  pdb_str_ref = """\
CRYST1  117.739  195.224  119.094  90.00 101.60  90.00 P 1 21 1
ATOM      1  N   THR G 332       4.195  72.012  51.895  1.00171.28           N
ATOM      2  CA  THR G 332       2.946  71.699  52.580  1.00153.71           C
ATOM      3  C   THR G 332       1.980  70.971  51.651  1.00132.45           C
ATOM      4  O   THR G 332       2.092  71.062  50.429  1.00130.38           O
ATOM      5  CB  THR G 332       2.291  72.982  53.125  1.00 20.00           C
ATOM      6  OG1 THR G 332       2.036  73.887  52.046  1.00 20.00           O
ATOM      7  CG2 THR G 332       3.269  73.749  54.003  1.00 20.00           C
ATOM      8  N   VAL G 333       1.033  70.248  52.240  1.00121.58           N
ATOM      9  CA  VAL G 333       0.047  69.503  51.468  1.00129.11           C
ATOM     10  C   VAL G 333      -1.363  69.905  51.883  1.00146.57           C
ATOM     11  O   VAL G 333      -1.552  70.599  52.882  1.00151.92           O
ATOM     12  CB  VAL G 333       0.216  67.983  51.643  1.00 20.00           C
ATOM     13  CG1 VAL G 333      -0.905  67.237  50.935  1.00 20.00           C
ATOM     14  CG2 VAL G 333       1.574  67.534  51.125  1.00 20.00           C
ATOM     15  N   SER G 334      -2.351  69.465  51.111  1.00156.08           N
ATOM     16  CA  SER G 334      -3.745  69.778  51.397  1.00158.16           C
ATOM     17  C   SER G 334      -4.297  68.870  52.492  1.00157.55           C
ATOM     18  O   SER G 334      -3.964  67.686  52.556  1.00161.93           O
ATOM     19  CB  SER G 334      -4.595  69.652  50.131  1.00162.09           C
ATOM     20  OG  SER G 334      -5.954  69.950  50.396  1.00170.98           O
  """
  params_text = """\
  reference_model {
    reference_group {
      reference = chain 'G'
      selection = chain 'K'
      file_name = "ref.pdb"
    }
  }
  """
  ref_file = open("ref.pdb", 'w')
  ref_file.write(pdb_str_ref)
  ref_file.close()
  log = StringIO()
  # log = sys.stdout
  # orig_file = open("start.pdb", "w")
  # orig_file.write(pdb_str_original)
  # orig_file.close()
  def_pars = reference_model_params
  pars = iotbx.phil.parse(params_text)
  all_pars = def_pars.fetch(pars).extract()
  all_pars.reference_model.enabled = True
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=flex.split_lines(pdb_str_original),
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()

  rm = reference_model(
         model=model,
         reference_file_list=['ref.pdb'],
         params=all_pars.reference_model,
         log=log)
  rm.show_reference_summary(log=log)
  new_h = pdb_h.deep_copy()
  xray_structure = new_h.extract_xray_structure()
  rm.set_rotamer_to_reference(
    xray_structure=xray_structure)
  new_h.adopt_xray_structure(xray_structure)
  r1 = rotalyze(pdb_hierarchy=new_h, outliers_only=False)
  assert r1.n_outliers == 0

def exercise_dna(mon_lib_srv, ener_lib):
  pdb_str_original = """\
CRYST1   25.287   40.217   65.471  90.00  90.00  90.00 P 21 21 21    8
SCALE1      0.039546  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024865  0.000000        0.00000
SCALE3      0.000000  0.000000  0.015274        0.00000
ATOM     80  P    DA A   5      -8.062  -5.965 -15.755  1.00 42.17           P
ATOM     81  OP1  DA A   5      -8.426  -7.228 -16.405  1.00 50.61           O
ATOM     82  OP2  DA A   5      -8.689  -5.557 -14.457  1.00 51.75           O
ATOM     83  O5'  DA A   5      -6.496  -5.961 -15.638  1.00 34.89           O
ATOM     84  C5'  DA A   5      -5.791  -6.321 -16.790  1.00 30.71           C
ATOM     85  C4'  DA A   5      -4.355  -5.917 -16.600  1.00 34.43           C
ATOM     86  O4'  DA A   5      -4.303  -4.509 -16.239  1.00 33.96           O
ATOM     87  C3'  DA A   5      -3.630  -6.687 -15.491  1.00 35.56           C
ATOM     88  O3'  DA A   5      -2.407  -7.257 -16.020  1.00 33.08           O
ATOM     89  C2'  DA A   5      -3.531  -5.654 -14.384  1.00 32.41           C
ATOM     90  C1'  DA A   5      -3.435  -4.334 -15.130  1.00 28.44           C
ATOM     91  N9   DA A   5      -3.904  -3.143 -14.449  1.00 28.37           N
ATOM     92  C8   DA A   5      -5.187  -2.933 -14.022  1.00 27.53           C
ATOM     93  N7   DA A   5      -5.401  -1.724 -13.565  1.00 29.33           N
ATOM     94  C5   DA A   5      -4.187  -1.082 -13.747  1.00 23.78           C
ATOM     95  C6   DA A   5      -3.761   0.226 -13.474  1.00 25.22           C
ATOM     96  N6   DA A   5      -4.519   1.150 -12.896  1.00 25.69           N
ATOM     97  N1   DA A   5      -2.485   0.535 -13.749  1.00 24.39           N
ATOM     98  C2   DA A   5      -1.712  -0.389 -14.320  1.00 24.89           C
ATOM     99  N3   DA A   5      -2.001  -1.641 -14.653  1.00 28.33           N
ATOM    100  C4   DA A   5      -3.268  -1.935 -14.326  1.00 27.45           C
ATOM    101  P    DA A   6      -1.382  -8.057 -15.083  1.00 33.49           P
ATOM    102  OP1  DA A   6      -0.596  -8.971 -15.989  1.00 35.26           O
ATOM    103  OP2  DA A   6      -2.097  -8.481 -13.890  1.00 34.48           O
ATOM    104  O5'  DA A   6      -0.480  -6.949 -14.401  1.00 31.72           O
ATOM    105  C5'  DA A   6       0.398  -6.138 -15.188  1.00 28.12           C
ATOM    106  C4'  DA A   6       1.219  -5.272 -14.269  1.00 22.57           C
ATOM    107  O4'  DA A   6       0.380  -4.203 -13.784  1.00 23.34           O
ATOM    108  C3'  DA A   6       1.783  -5.982 -13.049  1.00 23.61           C
ATOM    109  O3'  DA A   6       3.202  -5.785 -13.150  1.00 22.60           O
ATOM    110  C2'  DA A   6       1.110  -5.289 -11.881  1.00 22.21           C
ATOM    111  C1'  DA A   6       0.653  -3.958 -12.418  1.00 20.89           C
ATOM    112  N9   DA A   6      -0.561  -3.398 -11.831  1.00 21.71           N
ATOM    113  C8   DA A   6      -1.777  -4.017 -11.666  1.00 23.62           C
ATOM    114  N7   DA A   6      -2.693  -3.249 -11.139  1.00 23.57           N
ATOM    115  C5   DA A   6      -2.071  -2.016 -11.029  1.00 20.29           C
ATOM    116  C6   DA A   6      -2.506  -0.774 -10.519  1.00 20.33           C
ATOM    117  N6   DA A   6      -3.763  -0.525 -10.122  1.00 20.36           N
ATOM    118  N1   DA A   6      -1.604   0.233 -10.486  1.00 20.84           N
ATOM    119  C2   DA A   6      -0.341  -0.023 -10.868  1.00 21.15           C
ATOM    120  N3   DA A   6       0.174  -1.126 -11.378  1.00 22.91           N
ATOM    121  C4   DA A   6      -0.746  -2.101 -11.433  1.00 20.00           C
ATOM    122  P    DT A   7       4.283  -6.215 -12.051  1.00 23.53           P
ATOM    123  OP1  DT A   7       5.598  -6.398 -12.780  1.00 27.73           O
ATOM    124  OP2  DT A   7       3.774  -7.297 -11.205  1.00 24.18           O
ATOM    125  O5'  DT A   7       4.350  -4.948 -11.106  1.00 22.94           O
ATOM    126  C5'  DT A   7       4.668  -3.709 -11.633  1.00 21.30           C
ATOM    127  C4'  DT A   7       4.525  -2.656 -10.580  1.00 20.84           C
ATOM    128  O4'  DT A   7       3.138  -2.512 -10.296  1.00 19.94           O
ATOM    129  C3'  DT A   7       5.205  -2.966  -9.250  1.00 20.02           C
ATOM    130  O3'  DT A   7       6.280  -2.035  -9.099  1.00 23.74           O
ATOM    131  C2'  DT A   7       4.144  -2.717  -8.200  1.00 19.47           C
ATOM    132  C1'  DT A   7       3.048  -2.015  -8.962  1.00 20.12           C
ATOM    133  N1   DT A   7       1.641  -2.197  -8.524  1.00 20.27           N
ATOM    134  C2   DT A   7       0.957  -1.108  -8.030  1.00 18.61           C
ATOM    135  O2   DT A   7       1.430   0.017  -7.926  1.00 19.56           O
ATOM    136  N3   DT A   7      -0.344  -1.365  -7.721  1.00 18.89           N
ATOM    137  C4   DT A   7      -1.018  -2.563  -7.836  1.00 21.94           C
ATOM    138  O4   DT A   7      -2.200  -2.640  -7.497  1.00 23.57           O
ATOM    139  C5   DT A   7      -0.226  -3.674  -8.271  1.00 18.09           C
ATOM    140  C7   DT A   7      -0.860  -5.022  -8.351  1.00 19.35           C
ATOM    141  C6   DT A   7       1.065  -3.446  -8.562  1.00 17.66           C
ATOM    142  P    DT A   8       7.284  -1.980  -7.857  1.00 26.43           P
ATOM    143  OP1  DT A   8       8.611  -1.444  -8.278  1.00 28.45           O
ATOM    144  OP2  DT A   8       7.248  -3.298  -7.198  1.00 27.17           O
ATOM    145  O5'  DT A   8       6.613  -0.927  -6.882  1.00 25.09           O
ATOM    146  C5'  DT A   8       6.357   0.403  -7.340  1.00 24.67           C
ATOM    147  C4'  DT A   8       5.543   1.125  -6.301  1.00 23.10           C
ATOM    148  O4'  DT A   8       4.228   0.541  -6.229  1.00 23.60           O
ATOM    149  C3'  DT A   8       6.127   1.057  -4.884  1.00 25.21           C
ATOM    150  O3'  DT A   8       6.507   2.380  -4.493  1.00 28.93           O
ATOM    151  C2'  DT A   8       5.018   0.434  -4.050  1.00 23.32           C
ATOM    152  C1'  DT A   8       3.795   0.667  -4.883  1.00 22.06           C
ATOM    153  N1   DT A   8       2.713  -0.291  -4.689  1.00 19.79           N
ATOM    154  C2   DT A   8       1.466   0.223  -4.414  1.00 18.40           C
ATOM    155  O2   DT A   8       1.263   1.399  -4.157  1.00 20.56           O
ATOM    156  N3   DT A   8       0.484  -0.716  -4.337  1.00 19.20           N
ATOM    157  C4   DT A   8       0.588  -2.075  -4.597  1.00 18.45           C
ATOM    158  O4   DT A   8      -0.397  -2.789  -4.538  1.00 21.38           O
ATOM    159  C5   DT A   8       1.920  -2.549  -4.859  1.00 17.02           C
ATOM    160  C7   DT A   8       2.126  -4.006  -5.116  1.00 20.50           C
ATOM    161  C6   DT A   8       2.895  -1.634  -4.959  1.00 19.29           C
TER     245       DG A  12
ATOM    325  P    DA B  17     -10.220   1.260  -1.207  1.00 27.94           P
ATOM    326  OP1  DA B  17     -11.370   2.143  -0.856  1.00 34.83           O
ATOM    327  OP2  DA B  17     -10.221   0.599  -2.553  1.00 31.17           O
ATOM    328  O5'  DA B  17      -8.842   2.020  -1.098  1.00 26.12           O
ATOM    329  C5'  DA B  17      -8.558   2.683   0.094  1.00 25.41           C
ATOM    330  C4'  DA B  17      -7.407   3.619  -0.107  1.00 26.38           C
ATOM    331  O4'  DA B  17      -6.208   2.886  -0.440  1.00 24.41           O
ATOM    332  C3'  DA B  17      -7.600   4.631  -1.214  1.00 27.57           C
ATOM    333  O3'  DA B  17      -6.972   5.834  -0.764  1.00 29.89           O
ATOM    334  C2'  DA B  17      -6.902   3.980  -2.406  1.00 26.29           C
ATOM    335  C1'  DA B  17      -5.771   3.225  -1.781  1.00 23.13           C
ATOM    336  N9   DA B  17      -5.444   1.986  -2.460  1.00 22.66           N
ATOM    337  C8   DA B  17      -6.295   0.942  -2.750  1.00 23.38           C
ATOM    338  N7   DA B  17      -5.700  -0.094  -3.288  1.00 20.62           N
ATOM    339  C5   DA B  17      -4.344   0.242  -3.234  1.00 20.59           C
ATOM    340  C6   DA B  17      -3.178  -0.447  -3.603  1.00 17.89           C
ATOM    341  N6   DA B  17      -3.184  -1.685  -4.072  1.00 20.22           N
ATOM    342  N1   DA B  17      -1.995   0.205  -3.497  1.00 19.61           N
ATOM    343  C2   DA B  17      -1.992   1.465  -3.030  1.00 20.38           C
ATOM    344  N3   DA B  17      -3.021   2.207  -2.621  1.00 20.80           N
ATOM    345  C4   DA B  17      -4.182   1.540  -2.774  1.00 19.17           C
ATOM    346  P    DA B  18      -6.994   7.132  -1.670  1.00 32.91           P
ATOM    347  OP1  DA B  18      -6.817   8.281  -0.798  1.00 37.55           O
ATOM    348  OP2  DA B  18      -8.060   7.037  -2.636  1.00 31.04           O
ATOM    349  O5'  DA B  18      -5.659   7.052  -2.535  1.00 30.20           O
ATOM    350  C5'  DA B  18      -4.377   7.074  -1.958  1.00 30.19           C
ATOM    351  C4'  DA B  18      -3.354   6.838  -3.036  1.00 28.09           C
ATOM    352  O4'  DA B  18      -3.424   5.481  -3.484  1.00 26.27           O
ATOM    353  C3'  DA B  18      -3.545   7.708  -4.286  1.00 29.73           C
ATOM    354  O3'  DA B  18      -2.469   8.627  -4.273  1.00 34.73           O
ATOM    355  C2'  DA B  18      -3.566   6.715  -5.433  1.00 27.32           C
ATOM    356  C1'  DA B  18      -3.010   5.448  -4.841  1.00 24.83           C
ATOM    357  N9   DA B  18      -3.488   4.196  -5.410  1.00 23.72           N
ATOM    358  C8   DA B  18      -4.794   3.799  -5.530  1.00 20.51           C
ATOM    359  N7   DA B  18      -4.937   2.581  -5.985  1.00 22.85           N
ATOM    360  C5   DA B  18      -3.636   2.147  -6.189  1.00 20.87           C
ATOM    361  C6   DA B  18      -3.111   0.950  -6.675  1.00 19.34           C
ATOM    362  N6   DA B  18      -3.852  -0.099  -7.028  1.00 21.88           N
ATOM    363  N1   DA B  18      -1.767   0.849  -6.776  1.00 19.77           N
ATOM    364  C2   DA B  18      -1.023   1.872  -6.374  1.00 21.42           C
ATOM    365  N3   DA B  18      -1.392   3.050  -5.910  1.00 22.19           N
ATOM    366  C4   DA B  18      -2.734   3.129  -5.836  1.00 21.41           C
ATOM    367  P    DT B  19      -2.064   9.546  -5.497  1.00 40.82           P
ATOM    368  OP1  DT B  19      -1.281  10.615  -4.939  1.00 44.52           O
ATOM    369  OP2  DT B  19      -3.292   9.787  -6.271  1.00 44.69           O
ATOM    370  O5'  DT B  19      -1.119   8.619  -6.355  1.00 30.72           O
ATOM    371  C5'  DT B  19       0.059   8.093  -5.804  1.00 29.16           C
ATOM    372  C4'  DT B  19       0.704   7.195  -6.832  1.00 26.15           C
ATOM    373  O4'  DT B  19      -0.129   6.045  -7.087  1.00 26.00           O
ATOM    374  C3'  DT B  19       0.941   7.859  -8.188  1.00 25.98           C
ATOM    375  O3'  DT B  19       2.343   7.877  -8.376  1.00 30.07           O
ATOM    376  C2'  DT B  19       0.207   6.968  -9.181  1.00 26.77           C
ATOM    377  C1'  DT B  19       0.036   5.665  -8.443  1.00 25.87           C
ATOM    378  N1   DT B  19      -1.122   4.839  -8.816  1.00 24.60           N
ATOM    379  C2   DT B  19      -0.906   3.556  -9.283  1.00 22.21           C
ATOM    380  O2   DT B  19       0.197   3.084  -9.451  1.00 22.06           O
ATOM    381  N3   DT B  19      -2.038   2.833  -9.519  1.00 22.04           N
ATOM    382  C4   DT B  19      -3.339   3.262  -9.380  1.00 21.81           C
ATOM    383  O4   DT B  19      -4.247   2.495  -9.615  1.00 24.16           O
ATOM    384  C5   DT B  19      -3.499   4.613  -8.891  1.00 22.25           C
ATOM    385  C7   DT B  19      -4.879   5.143  -8.663  1.00 23.26           C
ATOM    386  C6   DT B  19      -2.396   5.327  -8.640  1.00 22.85           C
ATOM    387  P    DT B  20       3.005   8.456  -9.725  1.00 32.03           P
ATOM    388  OP1  DT B  20       4.339   8.958  -9.284  1.00 35.31           O
ATOM    389  OP2  DT B  20       2.027   9.351 -10.442  1.00 33.99           O
ATOM    390  O5'  DT B  20       3.144   7.102 -10.543  1.00 31.33           O
ATOM    391  C5'  DT B  20       3.894   5.979 -10.032  1.00 28.60           C
ATOM    392  C4'  DT B  20       3.851   4.840 -11.020  1.00 28.63           C
ATOM    393  O4'  DT B  20       2.494   4.361 -11.145  1.00 26.47           O
ATOM    394  C3'  DT B  20       4.300   5.211 -12.437  1.00 31.59           C
ATOM    395  O3'  DT B  20       5.260   4.256 -12.875  1.00 39.07           O
ATOM    396  C2'  DT B  20       3.027   5.147 -13.257  1.00 26.06           C
ATOM    397  C1'  DT B  20       2.211   4.120 -12.529  1.00 24.42           C
ATOM    398  N1   DT B  20       0.757   4.123 -12.660  1.00 23.79           N
ATOM    399  C2   DT B  20       0.138   2.932 -12.972  1.00 25.04           C
ATOM    400  O2   DT B  20       0.741   1.921 -13.262  1.00 24.66           O
ATOM    401  N3   DT B  20      -1.229   2.977 -12.959  1.00 25.84           N
ATOM    402  C4   DT B  20      -2.022   4.071 -12.671  1.00 25.98           C
ATOM    403  O4   DT B  20      -3.234   3.948 -12.646  1.00 28.14           O
ATOM    404  C5   DT B  20      -1.311   5.298 -12.387  1.00 22.81           C
ATOM    405  C7   DT B  20      -2.094   6.540 -12.092  1.00 27.47           C
ATOM    406  C6   DT B  20       0.028   5.263 -12.401  1.00 26.29           C
TER     490       DG B  24
  """
  params_text = """\
  reference_model {
    reference_group {
      reference = chain 'A'
      selection = chain 'A'
      file_name = "ref.pdb"
    }
    reference_group {
      reference = chain 'B'
      selection = chain 'B'
      file_name = "ref.pdb"
    }
  }
  """
  ref_file = open("ref.pdb", 'w')
  ref_file.write(pdb_str_original)
  ref_file.close()
  log = StringIO()
  # log = sys.stdout
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=flex.split_lines(pdb_str_original),
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()
  for include_chains in [True, False]:
    def_pars = reference_model_params
    pars = iotbx.phil.parse(params_text)
    all_pars = None
    if include_chains:
      all_pars = def_pars.fetch(pars).extract()
      all_pars.reference_model.enabled = True
    else:
      all_pars = def_pars.extract()
      all_pars.reference_model.enabled = True
      all_pars.reference_model.file = "ref.pdb"
    rm = reference_model(
           model=model,
           reference_file_list=['ref.pdb'],
           params=all_pars.reference_model,
           log=log)
    rm.show_reference_summary(log=log)
    # more dihedrals via sugar
    np=90
    assert rm.get_n_proxies() == np, \
        "Expecting %s proxies, got %d" % (np, rm.get_n_proxies())
    log_strings = log.getvalue().split("\n")
    for needed_string in [
        " DA A   5  <=====>   DA A   5",
        " DA A   6  <=====>   DA A   6",
        " DT A   7  <=====>   DT A   7",
        " DT A   8  <=====>   DT A   8",
        " DA B  17  <=====>   DA B  17",
        " DA B  18  <=====>   DA B  18",
        " DT B  19  <=====>   DT B  19",
        " DT B  20  <=====>   DT B  20",
        ]:
      assert needed_string in log_strings, "'%s' not in log!" % needed_string

def exercise_3chains_self(mon_lib_srv, ener_lib):
  """
  Test reference model, 3 chains, reference is the same and selections are
  supposed to be mixed, like ref=(A or Bref or C) sel=(Aref, B, Cref)  """
  pdb_str_original = """\
CRYST1  129.069   83.165   84.393  90.00  90.00  90.00 P 1
ATOM      1  N   GLN A   1     118.638  78.165  29.859  1.00 32.70      A    N
ATOM      2  CA  GLN A   1     118.742  77.022  30.759  1.00 34.10      A    C
ATOM      3  CB  GLN A   1     117.844  77.222  31.984  1.00 34.85      A    C
ATOM      4  CG  GLN A   1     118.008  76.159  33.064  1.00 36.31      A    C
ATOM      5  CD  GLN A   1     117.167  76.441  34.295  1.00 37.02      A    C
ATOM      6  OE1 GLN A   1     116.477  77.458  34.372  1.00 36.35      A    O
ATOM      7  NE2 GLN A   1     117.221  75.538  35.268  1.00 38.45      A    N
ATOM      8  C   GLN A   1     118.377  75.725  30.039  1.00 35.48      A    C
ATOM      9  O   GLN A   1     119.251  75.008  29.552  1.00 35.59      A    O
ATOM     10  N   VAL A   2     117.083  75.432  29.969  1.00 36.49      A    N
ATOM     11  CA  VAL A   2     116.607  74.213  29.327  1.00 37.63      A    C
ATOM     12  CB  VAL A   2     115.168  73.893  29.751  1.00 38.86      A    C
ATOM     13  CG1 VAL A   2     114.654  72.672  29.006  1.00 39.72      A    C
ATOM     14  CG2 VAL A   2     115.096  73.682  31.255  1.00 39.93      A    C
ATOM     15  C   VAL A   2     116.698  74.327  27.811  1.00 36.86      A    C
ATOM     16  O   VAL A   2     116.042  75.175  27.207  1.00 36.25      A    O
ATOM     17  N   GLN A   3     117.506  73.466  27.200  1.00 37.05      A    N
ATOM     18  CA  GLN A   3     117.678  73.479  25.752  1.00 36.50      A    C
ATOM     19  CB  GLN A   3     118.915  74.294  25.365  1.00 35.23      A    C
ATOM     20  CG  GLN A   3     118.741  75.798  25.481  1.00 33.97      A    C
ATOM     21  CD  GLN A   3     119.955  76.561  24.994  1.00 32.60      A    C
ATOM     22  OE1 GLN A   3     120.999  75.974  24.707  1.00 32.65      A    O
ATOM     23  NE2 GLN A   3     119.823  77.879  24.893  1.00 31.40      A    N
ATOM     24  C   GLN A   3     117.794  72.072  25.176  1.00 37.49      A    C
ATOM     25  O   GLN A   3     118.315  71.161  25.827  1.00 38.37      A    O
ATOM     26  N   LEU A   4     117.302  71.907  23.951  1.00 37.44      A    N
ATOM     27  CA  LEU A   4     117.424  70.652  23.217  1.00 38.27      A    C
ATOM     28  CB  LEU A   4     116.107  69.872  23.227  1.00 39.05      A    C
ATOM     29  CG  LEU A   4     115.539  69.386  24.562  1.00 39.82      A    C
ATOM     30  CD1 LEU A   4     114.720  70.467  25.254  1.00 39.45      A    C
ATOM     31  CD2 LEU A   4     114.704  68.133  24.355  1.00 40.68      A    C
ATOM     32  C   LEU A   4     117.854  70.938  21.778  1.00 37.77      A    C
ATOM     33  O   LEU A   4     117.369  71.884  21.157  1.00 37.07      A    O
ATOM     34  N   LYS A   5     118.763  70.124  21.249  1.00 38.27      A    N
ATOM     35  CA  LYS A   5     119.265  70.333  19.895  1.00 37.96      A    C
ATOM     36  CB  LYS A   5     120.574  71.123  19.931  1.00 37.03      A    C
ATOM     37  CG  LYS A   5     121.114  71.505  18.561  1.00 36.53      A    C
ATOM     38  CD  LYS A   5     122.352  72.380  18.680  1.00 35.38      A    C
ATOM     39  CE  LYS A   5     122.875  72.783  17.311  1.00 34.77      A    C
ATOM     40  NZ  LYS A   5     124.069  73.666  17.413  1.00 33.45      A    N
ATOM     41  C   LYS A   5     119.467  69.009  19.161  1.00 39.11      A    C
ATOM     42  O   LYS A   5     120.063  68.075  19.695  1.00 39.98      A    O
ATOM     43  N   GLU A   6     118.969  68.936  17.931  1.00 39.24      A    N
ATOM     44  CA  GLU A   6     119.059  67.716  17.133  1.00 40.35      A    C
ATOM     45  CB  GLU A   6     117.809  67.541  16.267  1.00 40.77      A    C
ATOM     46  CG  GLU A   6     116.518  67.353  17.046  1.00 40.94      A    C
ATOM     47  CD  GLU A   6     115.905  68.664  17.505  1.00 39.97      A    C
ATOM     48  OE1 GLU A   6     116.574  69.714  17.400  1.00 39.09      A    O
ATOM     49  OE2 GLU A   6     114.743  68.644  17.962  1.00 40.09      A    O
ATOM     50  C   GLU A   6     120.296  67.709  16.241  1.00 40.24      A    C
ATOM     51  O   GLU A   6     120.601  68.697  15.574  1.00 39.29      A    O
ATOM     52  N   SER A   7     121.001  66.584  16.234  1.00 41.20      A    N
ATOM     53  CA  SER A   7     122.160  66.398  15.370  1.00 41.10      A    C
ATOM     54  CB  SER A   7     123.431  66.187  16.197  1.00 41.10      A    C
ATOM     55  OG  SER A   7     123.701  67.308  17.023  1.00 40.11      A    O
ATOM     56  C   SER A   7     121.930  65.211  14.444  1.00 42.23      A    C
ATOM     57  O   SER A   7     121.956  64.061  14.881  1.00 43.35      A    O
ATOM     58  N   GLY A   8     121.697  65.499  13.167  1.00 42.03      A    N
ATOM     59  CA  GLY A   8     121.422  64.464  12.188  1.00 43.17      A    C
ATOM     60  C   GLY A   8     122.302  64.540  10.956  1.00 42.71      A    C
ATOM     61  O   GLY A   8     123.142  65.433  10.843  1.00 41.37      A    O
ATOM     62  N   PRO A   9     122.117  63.592  10.024  1.00 43.74      A    N
ATOM     63  CD  PRO A   9     121.218  62.434  10.175  1.00 45.40      A    C
ATOM     64  CA  PRO A   9     122.906  63.513   8.790  1.00 43.28      A    C
ATOM     65  CB  PRO A   9     122.804  62.037   8.414  1.00 44.75      A    C
ATOM     66  CG  PRO A   9     121.465  61.631   8.919  1.00 46.12      A    C
ATOM     67  C   PRO A   9     122.361  64.390   7.666  1.00 42.79      A    C
ATOM     68  O   PRO A   9     123.122  64.826   6.800  1.00 41.70      A    O
ATOM     69  N   GLY A  10     121.055  64.637   7.679  1.00 43.56      A    N
ATOM     70  CA  GLY A  10     120.425  65.454   6.659  1.00 43.32      A    C
ATOM     71  C   GLY A  10     119.971  64.668   5.445  1.00 44.47      A    C
ATOM     72  O   GLY A  10     118.831  64.801   5.000  1.00 45.36      A    O
TER
ATOM   1645  N   ASP B   1      94.462  51.713  21.314  1.00 38.68      B    N
ATOM   1646  CA  ASP B   1      94.907  52.727  20.365  1.00 39.77      B    C
ATOM   1647  CB  ASP B   1      94.995  52.139  18.956  1.00 40.01      B    C
ATOM   1648  CG  ASP B   1      95.754  53.035  18.000  1.00 41.28      B    C
ATOM   1649  OD1 ASP B   1      96.565  53.857  18.476  1.00 41.67      B    O
ATOM   1650  OD2 ASP B   1      95.544  52.913  16.774  1.00 41.91      B    O
ATOM   1651  C   ASP B   1      93.966  53.928  20.386  1.00 40.53      B    C
ATOM   1652  O   ASP B   1      92.746  53.766  20.440  1.00 40.33      B    O
ATOM   1653  N   ILE B   2      94.537  55.130  20.341  1.00 41.47      B    N
ATOM   1654  CA  ILE B   2      93.756  56.358  20.452  1.00 42.47      B    C
ATOM   1655  CB  ILE B   2      94.214  57.190  21.663  1.00 42.52      B    C
ATOM   1656  CG2 ILE B   2      93.396  58.470  21.778  1.00 43.79      B    C
ATOM   1657  CG1 ILE B   2      94.108  56.362  22.947  1.00 41.45      B    C
ATOM   1658  CD1 ILE B   2      94.570  57.093  24.191  1.00 41.67      B    C
ATOM   1659  C   ILE B   2      93.841  57.180  19.168  1.00 43.93      B    C
ATOM   1660  O   ILE B   2      94.928  57.405  18.636  1.00 44.25      B    O
ATOM   1661  N   VAL B   3      92.687  57.629  18.680  1.00 45.00      B    N
ATOM   1662  CA  VAL B   3      92.612  58.382  17.432  1.00 46.83      B    C
ATOM   1663  CB  VAL B   3      91.712  57.671  16.414  1.00 47.07      B    C
ATOM   1664  CG1 VAL B   3      91.747  58.390  15.074  1.00 49.24      B    C
ATOM   1665  CG2 VAL B   3      92.136  56.224  16.262  1.00 45.38      B    C
ATOM   1666  C   VAL B   3      92.111  59.799  17.681  1.00 48.52      B    C
ATOM   1667  O   VAL B   3      91.117  60.002  18.375  1.00 48.59      B    O
ATOM   1668  N   MET B   4      92.799  60.776  17.099  1.00 49.99      B    N
ATOM   1669  CA  MET B   4      92.460  62.181  17.296  1.00 51.83      B    C
ATOM   1670  CB  MET B   4      93.658  62.944  17.866  1.00 51.35      B    C
ATOM   1671  CG  MET B   4      94.339  62.254  19.040  1.00 48.92      B    C
ATOM   1672  SD  MET B   4      93.360  62.268  20.552  1.00 48.39      B    S
ATOM   1673  CE  MET B   4      93.448  63.998  20.987  1.00 49.97      B    C
ATOM   1674  C   MET B   4      92.005  62.828  15.992  1.00 54.57      B    C
ATOM   1675  O   MET B   4      92.701  62.758  14.979  1.00 55.32      B    O
ATOM   1676  N   SER B   5      90.836  63.461  16.023  1.00 56.25      B    N
ATOM   1677  CA  SER B   5      90.302  64.140  14.847  1.00 59.25      B    C
ATOM   1678  CB  SER B   5      89.052  63.420  14.334  1.00 59.12      B    C
ATOM   1679  OG  SER B   5      89.335  62.070  14.010  1.00 56.93      B    O
ATOM   1680  C   SER B   5      89.975  65.596  15.161  1.00 61.44      B    C
ATOM   1681  O   SER B   5      89.374  65.889  16.191  1.00 60.75      B    O
ATOM   1682  N   GLN B   6      90.363  66.508  14.274  1.00 64.02      B    N
ATOM   1683  CA  GLN B   6      90.108  67.930  14.492  1.00 65.54      B    C
ATOM   1684  CB  GLN B   6      91.422  68.711  14.520  1.00 64.95      B    C
ATOM   1685  CG  GLN B   6      92.351  68.313  15.648  1.00 62.00      B    C
ATOM   1686  CD  GLN B   6      93.557  69.221  15.762  1.00 61.07      B    C
ATOM   1687  OE1 GLN B   6      94.648  68.778  16.125  1.00 58.90      B    O
ATOM   1688  NE2 GLN B   6      93.368  70.500  15.458  1.00 62.20      B    N
ATOM   1689  C   GLN B   6      89.182  68.526  13.435  1.00 68.28      B    C
ATOM   1690  O   GLN B   6      89.240  68.157  12.261  1.00 69.90      B    O
ATOM   1691  N   SER B   7      88.332  69.456  13.861  1.00 68.57      B    N
ATOM   1692  CA  SER B   7      87.413  70.129  12.949  1.00 70.24      B    C
ATOM   1693  CB  SER B   7      86.049  69.433  12.944  1.00 69.94      B    C
ATOM   1694  OG  SER B   7      86.154  68.099  12.477  1.00 69.41      B    O
ATOM   1695  C   SER B   7      87.253  71.595  13.333  1.00 69.83      B    C
ATOM   1696  O   SER B   7      87.048  71.909  14.503  1.00 68.66      B    O
ATOM   1697  N   PRO B   8      87.340  72.500  12.345  1.00 70.37      B    N
ATOM   1698  CD  PRO B   8      87.103  73.940  12.555  1.00 69.33      B    C
ATOM   1699  CA  PRO B   8      87.579  72.191  10.932  1.00 71.70      B    C
ATOM   1700  CB  PRO B   8      86.996  73.407  10.217  1.00 71.69      B    C
ATOM   1701  CG  PRO B   8      87.247  74.525  11.170  1.00 69.92      B    C
ATOM   1702  C   PRO B   8      89.059  72.030  10.600  1.00 71.13      B    C
ATOM   1703  O   PRO B   8      89.910  72.310  11.444  1.00 69.82      B    O
ATOM   1704  N   SER B   9      89.354  71.585   9.382  1.00 71.87      B    N
ATOM   1705  CA  SER B   9      90.732  71.404   8.940  1.00 70.91      B    C
ATOM   1706  CB  SER B   9      90.775  70.616   7.629  1.00 71.45      B    C
ATOM   1707  OG  SER B   9      89.977  71.234   6.633  1.00 71.72      B    O
ATOM   1708  C   SER B   9      91.432  72.749   8.770  1.00 69.14      B    C
ATOM   1709  O   SER B   9      92.628  72.878   9.037  1.00 67.60      B    O
ATOM   1710  N   SER B  10      90.677  73.747   8.324  1.00 69.19      B    N
ATOM   1711  CA  SER B  10      91.201  75.095   8.150  1.00 67.65      B    C
ATOM   1712  CB  SER B  10      92.017  75.200   6.860  1.00 66.83      B    C
ATOM   1713  OG  SER B  10      91.215  74.927   5.723  1.00 67.75      B    O
ATOM   1714  C   SER B  10      90.055  76.097   8.134  1.00 68.08      B    C
ATOM   1715  O   SER B  10      88.958  75.786   7.670  1.00 69.49      B    O
TER
ATOM   3353  N   GLN C   1      27.855   6.390  79.393  1.00 55.82      C    N
ATOM   3354  CA  GLN C   1      27.377   6.759  78.009  1.00 57.48      C    C
ATOM   3355  CB  GLN C   1      26.126   5.903  77.650  1.00 57.65      C    C
ATOM   3356  CG  GLN C   1      24.762   6.447  78.162  1.00 55.80      C    C
ATOM   3357  CD  GLN C   1      23.623   5.432  77.999  1.00 55.35      C    C
ATOM   3358  OE1 GLN C   1      22.972   5.032  78.969  1.00 53.50      C    O
ATOM   3359  NE2 GLN C   1      23.365   5.000  76.745  1.00 56.64      C    N
ATOM   3360  C   GLN C   1      27.097   8.250  77.886  1.00 56.46      C    C
ATOM   3361  O   GLN C   1      26.949   8.930  78.891  1.00 54.64      C    O
ATOM   3362  N   VAL C   2      27.019   8.808  76.660  1.00 57.14      C    N
ATOM   3363  CA  VAL C   2      26.719  10.217  76.428  1.00 55.41      C    C
ATOM   3364  CB  VAL C   2      27.931  10.987  75.916  1.00 55.45      C    C
ATOM   3365  CG1 VAL C   2      27.552  12.336  75.269  1.00 53.41      C    C
ATOM   3366  CG2 VAL C   2      28.864  11.254  77.110  1.00 54.84      C    C
ATOM   3367  C   VAL C   2      25.580  10.337  75.447  1.00 54.98      C    C
ATOM   3368  O   VAL C   2      25.617   9.750  74.367  1.00 56.38      C    O
ATOM   3369  N   GLN C   3      24.516  11.078  75.796  1.00 52.87      C    N
ATOM   3370  CA  GLN C   3      23.345  11.230  74.957  1.00 51.91      C    C
ATOM   3371  CB  GLN C   3      22.235  10.208  75.319  1.00 52.11      C    C
ATOM   3372  CG  GLN C   3      22.651   8.728  75.133  1.00 54.78      C    C
ATOM   3373  CD  GLN C   3      21.498   7.771  75.463  1.00 54.32      C    C
ATOM   3374  OE1 GLN C   3      20.976   7.743  76.584  1.00 52.76      C    O
ATOM   3375  NE2 GLN C   3      21.093   6.936  74.478  1.00 54.60      C    N
ATOM   3376  C   GLN C   3      22.755  12.621  75.095  1.00 49.15      C    C
ATOM   3377  O   GLN C   3      22.938  13.298  76.106  1.00 48.04      C    O
ATOM   3378  N   LEU C   4      22.017  13.068  74.065  1.00 47.90      C    N
ATOM   3379  CA  LEU C   4      21.268  14.320  74.040  1.00 45.28      C    C
ATOM   3380  CB  LEU C   4      21.893  15.309  73.052  1.00 44.69      C    C
ATOM   3381  CG  LEU C   4      23.323  15.776  73.331  1.00 45.36      C    C
ATOM   3382  CD1 LEU C   4      23.813  16.685  72.216  1.00 44.60      C    C
ATOM   3383  CD2 LEU C   4      23.407  16.483  74.670  1.00 44.38      C    C
ATOM   3384  C   LEU C   4      19.813  14.054  73.665  1.00 43.91      C    C
ATOM   3385  O   LEU C   4      19.518  13.665  72.535  1.00 44.12      C    O
ATOM   3386  N   GLN C   5      18.907  14.263  74.615  1.00 42.26      C    N
ATOM   3387  CA  GLN C   5      17.488  14.011  74.382  1.00 40.51      C    C
ATOM   3388  CB  GLN C   5      16.852  13.347  75.606  1.00 39.88      C    C
ATOM   3389  CG  GLN C   5      17.529  12.057  76.033  1.00 42.29      C    C
ATOM   3390  CD  GLN C   5      17.525  11.010  74.938  1.00 43.78      C    C
ATOM   3391  OE1 GLN C   5      18.566  10.448  74.598  1.00 46.37      C    O
ATOM   3392  NE2 GLN C   5      16.350  10.738  74.382  1.00 42.03      C    N
ATOM   3393  C   GLN C   5      16.760  15.304  74.047  1.00 38.05      C    C
ATOM   3394  O   GLN C   5      16.865  16.282  74.776  1.00 37.00      C    O
ATOM   3395  N   GLN C   6      16.020  15.312  72.945  1.00 37.12      C    N
ATOM   3396  CA  GLN C   6      15.338  16.529  72.520  1.00 34.96      C    C
ATOM   3397  CB  GLN C   6      15.647  16.828  71.051  1.00 35.40      C    C
ATOM   3398  CG  GLN C   6      17.121  17.076  70.774  1.00 37.30      C    C
ATOM   3399  CD  GLN C   6      17.370  17.631  69.387  1.00 37.26      C    C
ATOM   3400  OE1 GLN C   6      18.353  17.285  68.732  1.00 38.80      C    O
ATOM   3401  NE2 GLN C   6      16.478  18.503  68.932  1.00 35.44      C    N
ATOM   3402  C   GLN C   6      13.830  16.452  72.732  1.00 32.51      C    C
ATOM   3403  O   GLN C   6      13.230  15.380  72.646  1.00 32.39      C    O
ATOM   3404  N   SER C   7      13.227  17.604  73.013  1.00 30.41      C    N
ATOM   3405  CA  SER C   7      11.789  17.692  73.226  1.00 27.70      C    C
ATOM   3406  CB  SER C   7      11.431  19.008  73.916  1.00 25.77      C    C
ATOM   3407  OG  SER C   7      11.822  20.115  73.122  1.00 25.94      C    O
ATOM   3408  C   SER C   7      11.031  17.575  71.909  1.00 26.72      C    C
ATOM   3409  O   SER C   7      11.633  17.430  70.846  1.00 28.21      C    O
ATOM   3410  N   GLY C   8       9.707  17.636  71.984  1.00 24.02      C    N
ATOM   3411  CA  GLY C   8       8.882  17.597  70.792  1.00 22.64      C    C
ATOM   3412  C   GLY C   8       8.215  16.260  70.529  1.00 22.13      C    C
ATOM   3413  O   GLY C   8       8.319  15.342  71.342  1.00 22.76      C    O
ATOM   3414  N   PRO C   9       7.522  16.141  69.385  1.00 20.87      C    N
ATOM   3415  CD  PRO C   9       6.878  14.882  68.970  1.00 20.19      C    C
ATOM   3416  CA  PRO C   9       7.356  17.200  68.381  1.00 20.00      C    C
ATOM   3417  CB  PRO C   9       6.793  16.446  67.173  1.00 19.29      C    C
ATOM   3418  CG  PRO C   9       6.076  15.287  67.766  1.00 18.19      C    C
ATOM   3419  C   PRO C   9       6.401  18.308  68.821  1.00 17.07      C    C
ATOM   3420  O   PRO C   9       5.485  18.062  69.605  1.00 14.88      C    O
ATOM   3421  N   GLU C  10       6.626  19.516  68.314  1.00 16.93      C    N
ATOM   3422  CA  GLU C  10       5.849  20.679  68.723  1.00 14.34      C    C
ATOM   3423  CB  GLU C  10       6.778  21.768  69.269  1.00 15.71      C    C
ATOM   3424  CG  GLU C  10       7.498  21.386  70.554  1.00 17.24      C    C
ATOM   3425  CD  GLU C  10       6.559  21.288  71.742  1.00 14.84      C    C
ATOM   3426  OE1 GLU C  10       5.821  22.262  71.999  1.00 12.41      C    O
ATOM   3427  OE2 GLU C  10       6.556  20.236  72.416  1.00 15.30      C    O
ATOM   3428  C   GLU C  10       5.000  21.240  67.583  1.00 12.13      C    C
ATOM   3429  O   GLU C  10       5.457  21.355  66.443  1.00 13.28      C    O
TER
END
  """
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=flex.split_lines(pdb_str_original),
                                    source_info=None))
  model.process(make_restraints=True)
  pdb_h = model.get_hierarchy()
  ref_h = pdb_h.deep_copy()
  # pdb_h.atoms().reset_i_seq()
  # ref_h.atoms().reset_i_seq()

  log = StringIO()
  # log = sys.stdout
  def_pars = reference_model_params
  all_pars = def_pars.fetch().extract()
  all_pars.reference_model.use_starting_model_as_reference=True
  all_pars.reference_model.enabled = True
  rm = reference_model(
         model=model,
         reference_hierarchy_list=\
            [model.get_hierarchy()],
         params=all_pars.reference_model,
         log=log)
  rm.show_reference_summary(log=log)
  assert rm.get_n_proxies() == 141, \
      "Expecting 141 proxies, got %d" % rm.get_n_proxies()
  log_strings = log.getvalue().split("\n")
  # print "========"
  # print "\n".join(log_strings)
  # print "========"
  for needed_string in [
      "GLY A   8  <=====>  GLY A   8",
      "PRO A   9  <=====>  PRO A   9",
      "GLY A  10  <=====>  GLY A  10",
      "ASP B   1  <=====>  ASP B   1",
      "ILE B   2  <=====>  ILE B   2",
      "SER B  10  <=====>  SER B  10",
      "GLN C   1  <=====>  GLN C   1",
      "VAL C   2  <=====>  VAL C   2",
      ]:
    assert needed_string in log_strings, "'%s' not in log!" % needed_string


def run(args):
  t0 = time.time()
  import mmtbx.monomer_library
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  exercise_reference_model(args, mon_lib_srv, ener_lib)
  exercise_multiple_to_one(args, mon_lib_srv, ener_lib)
  exercise_multiple_ncs_groups_found(mon_lib_srv, ener_lib)
  exercise_cutted_residue(mon_lib_srv, ener_lib)
  exercise_dna(mon_lib_srv, ener_lib)
  exercise_3chains_self(mon_lib_srv, ener_lib)
  print("OK. Time: %8.3f"%(time.time()-t0))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
