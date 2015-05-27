from __future__ import division
from mmtbx import monomer_library
from mmtbx.torsion_restraints.reference_model import reference_model, reference_model_params
from mmtbx.torsion_restraints import utils
from mmtbx.validation.rotalyze import rotalyze
from cctbx.array_family import flex
import iotbx.phil
import iotbx.utils
import iotbx.pdb
from libtbx.test_utils import show_diff
import libtbx.load_env
import cStringIO
import sys, os, time

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

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
reference_model
{
  include \
    scope mmtbx.torsion_restraints.reference_model.reference_model_params
}
""", process_includes=True)

def exercise_reference_model(args, mon_lib_srv, ener_lib):
  log = cStringIO.StringIO()
  master_phil = get_master_phil()
  input_objects = iotbx.utils.process_command_line_inputs(
    args=args,
    master_phil=master_phil,
    input_types=("mtz", "pdb", "cif"))
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  master_phil_str_overrides = """
  reference_model {
    fix_outliers=False
  }
  """
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  work_params = master_phil.fetch(sources=phil_objects).extract()
  pdb_hierarchy = iotbx.pdb.input(
    source_info=None,
    lines=flex.split_lines(model_raw_records)).construct_hierarchy()
  reference_hierarchy_list = []
  tmp_hierarchy = iotbx.pdb.input(
    source_info=None,
    lines=flex.split_lines(reference_raw_records)).construct_hierarchy()

  reference_hierarchy_list.append(tmp_hierarchy)
  rm = reference_model(
         pdb_hierarchy=pdb_hierarchy,
         reference_hierarchy_list=reference_hierarchy_list,
         params=work_params.reference_model,
         log=log)
  assert len(rm.reference_dihedral_proxies) == 7

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
    pdb_hierarchy=pdb_hierarchy)
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
    pdb_hierarchy=pdb_hierarchy)
  assert i_seq_element_hash == \
    {0: ' N', 1: ' C', 2: ' C', 3: ' O', 4: ' C', 5: ' C', 6: ' O', 7: ' N',
     8: ' N', 9: ' C', 10: ' C', 11: ' O', 12: ' C', 13: ' C', 14: ' C',
     15: ' C'}

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
  assert len(dihedral_hash) == 7
  reference_dihedral_proxies = rm.reference_dihedral_proxies.deep_copy()
  assert reference_dihedral_proxies is not None
  assert len(reference_dihedral_proxies) == len(dihedral_hash)
  for rdp in reference_dihedral_proxies:
    assert rdp.limit == work_params.reference_model.limit

  r1 = rotalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  out1 = cStringIO.StringIO()
  r1.show_old_output(out=out1)
  r2 = rotalyze(pdb_hierarchy=ref_pdb_hierarchy, outliers_only=False)
  out2 = cStringIO.StringIO()
  r2.show_old_output(out=out2)

  assert not show_diff(out1.getvalue(), """\
 C 236  ASN:1.00:0.2:227.3:80.2:::OUTLIER:OUTLIER
 C 237  LEU:1.00:0.0:209.6:357.2:::OUTLIER:OUTLIER
""")

  assert not show_diff(out2.getvalue(), """\
 C 236  ASN:1.00:39.1:203.2:43.6:::Favored:t0
 C 237  LEU:1.00:60.8:179.1:57.3:::Favored:tp
""")

  xray_structure = pdb_hierarchy.extract_xray_structure()
  rm.set_rotamer_to_reference(
    xray_structure=xray_structure,
    mon_lib_srv=mon_lib_srv,
    quiet=True)
  pdb_hierarchy.adopt_xray_structure(xray_structure)
  r2 = rotalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  out3 = cStringIO.StringIO()
  r2.show_old_output(out=out3)
  assert not show_diff(out3.getvalue(), """\
 C 236  ASN:1.00:39.1:203.2:43.6:::Favored:t0
 C 237  LEU:1.00:60.8:179.1:57.3:::Favored:tp
""")

  match_map = rm.match_map['ref1']
  assert match_map == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11,
   12: 12, 13: 13, 14: 14, 15: 15}

  master_phil_str_overrides = """
  reference_model.reference_group {
    reference= chain B and resseq 246:247
    selection= chain C and resid 236:237
  }
  """
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  work_params_alt = master_phil.fetch(sources=phil_objects).extract()
  rm = reference_model(
         pdb_hierarchy=pdb_hierarchy,
         reference_hierarchy_list=reference_hierarchy_list_alt_seq,
         params=work_params_alt.reference_model,
         log=log)
  match_map = rm.match_map
  assert match_map['ref1'] == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11,
   12: 12, 13: 13, 14: 14, 15: 15}

  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_hierarchy = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
  reference_file_list = []
  reference_file_list.append(pdb_file)
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  master_phil_str_overrides = """
  reference_model {
    fix_outliers=False
  }
  """
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  work_params = master_phil.fetch(sources=phil_objects).extract()
  rm = reference_model(
         pdb_hierarchy=pdb_hierarchy,
         reference_file_list=reference_file_list,
         params=work_params.reference_model,
         log=log)
  reference_dihedral_proxies = rm.reference_dihedral_proxies
  standard_weight = 0
  for dp in reference_dihedral_proxies:
    if dp.weight == 1.0:
      standard_weight += 1
  assert standard_weight == 1181
  if (not libtbx.env.has_module(name="ksdssp")):
    print "Skipping KSDSSP tests: ksdssp module not available."
  else:
    master_phil_str_overrides = """
    reference_model {
      secondary_structure_only = True
    }
    """
    phil_objects = [
      iotbx.phil.parse(input_string=master_phil_str_overrides)]
    work_params_ss = master_phil.fetch(sources=phil_objects).extract()
    rm.params = work_params_ss.reference_model
    rm.get_reference_dihedral_proxies()
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
  log = cStringIO.StringIO()
  # log = sys.stdout
  # orig_file = open("start.pdb", "w")
  # orig_file.write(pdb_str_original)
  # orig_file.close()
  def_pars = reference_model_params
  params_text = """\
 reference_model {
    file = "C550C_M14471_modrefiner_pdbset1.pdb"
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
  inp_hierarchy = iotbx.pdb.input(
      source_info=None,
      lines=flex.split_lines(pdb_str_original)).construct_hierarchy()
  rm = reference_model(
         pdb_hierarchy=inp_hierarchy,
         reference_file_list=['ref.pdb'],
         mon_lib_srv=mon_lib_srv,
         ener_lib=ener_lib,
         params=all_pars,
         log=log)
  # rm.show_reference_summary(log=log)
  new_h = inp_hierarchy.deep_copy()
  xray_structure = new_h.extract_xray_structure()
  rm.set_rotamer_to_reference(
    xray_structure=xray_structure)
  new_h.adopt_xray_structure(xray_structure)
  r1 = rotalyze(pdb_hierarchy=new_h, outliers_only=False)
  assert r1.n_outliers == 0
  # new_h.write_pdb_file(file_name="final.pdb")


def run(args):
  t0 = time.time()
  import mmtbx.monomer_library
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  exercise_reference_model(args, mon_lib_srv, ener_lib)
  exercise_multiple_to_one(args, mon_lib_srv, ener_lib)
  print "OK. Time: %8.3f"%(time.time()-t0)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
