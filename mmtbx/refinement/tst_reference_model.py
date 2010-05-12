import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from mmtbx.refinement import reference_model
from mmtbx.validation.rotalyze import rotalyze
import iotbx.phil
import iotbx.utils
import sys

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
""".splitlines()

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
""".splitlines()

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
reference_model
{
    file = None
      .type = str
    selection = None
      .type = str
    sigma = 1.0
      .type = float
    limit = 15.0
      .type = float
    slack = 0
      .type = float
    hydrogens = False
      .type = bool
    main_chain = True
      .type = bool
    side_chain = True
      .type = bool
}
""")

def exercise_reference_model(args, mon_lib_srv, ener_lib):
  master_phil = get_master_phil()
  input_objects = iotbx.utils.process_command_line_inputs(
    args=args,
    master_phil=master_phil,
    input_types=("mtz", "pdb", "cif"))
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  work_params = work_phil.extract()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=model_raw_records,
    for_dihedral_reference=True)
  geometry = processed_pdb_file.geometry_restraints_manager()
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  xray_structure=processed_pdb_file.xray_structure()
  processed_pdb_file_ref = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=reference_raw_records,
    for_dihedral_reference=True)
  geometry_ref = processed_pdb_file_ref.geometry_restraints_manager()
  sites_cart_ref = processed_pdb_file_ref.all_chain_proxies.sites_cart

  i_seq_name_hash = reference_model.build_name_hash(
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  assert i_seq_name_hash == \
    {0: ' N   ASN C 236 ', 1: ' CA  ASN C 236 ', 2: ' C   ASN C 236 ',
     3: ' O   ASN C 236 ', 4: ' CB  ASN C 236 ', 5: ' CG  ASN C 236 ',
     6: ' OD1 ASN C 236 ', 7: ' ND2 ASN C 236 ', 8: ' N   LEU C 237 ',
     9: ' CA  LEU C 237 ', 10: ' C   LEU C 237 ', 11: ' O   LEU C 237 ',
     12: ' CB  LEU C 237 ', 13: ' CG  LEU C 237 ', 14: ' CD1 LEU C 237 ',
     15: ' CD2 LEU C 237 '}

  i_seq_element_hash = reference_model.build_element_hash(
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  assert i_seq_element_hash == \
    {0: ' N', 1: ' C', 2: ' C', 3: ' O', 4: ' C', 5: ' C', 6: ' O', 7: ' N',
     8: ' N', 9: ' C', 10: ' C', 11: ' O', 12: ' C', 13: ' C', 14: ' C',
     15: ' C'}

  dihedral_hash = reference_model.build_dihedral_hash(
    geometry=geometry_ref,
    sites_cart=sites_cart_ref,
    pdb_hierarchy=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy,
    include_hydrogens=False,
    include_main_chain=True,
    include_side_chain=True)
  assert dihedral_hash == \
    {' CA  ASN C 236  C   ASN C 236  N   LEU C 237  CA  LEU C 237 ': -171.86697238051445,
     ' C   ASN C 236  N   LEU C 237  CA  LEU C 237  C   LEU C 237 ': -137.00198604910656,
     ' CA  LEU C 237  CB  LEU C 237  CG  LEU C 237  CD2 LEU C 237 ': -178.92930350265104,
     ' CA  ASN C 236  CB  ASN C 236  CG  ASN C 236  OD1 ASN C 236 ': 43.580408297861439,
     ' N   ASN C 236  CA  ASN C 236  C   ASN C 236  N   LEU C 237 ': 29.747422877599458,
     ' CA  ASN C 236  N   ASN C 236  C   ASN C 236  CB  ASN C 236 ': -123.34365764907594,
     ' N   LEU C 237  CA  LEU C 237  CB  LEU C 237  CG  LEU C 237 ': 179.11361522972217,
     ' CA  LEU C 237  N   LEU C 237  C   LEU C 237  CB  LEU C 237 ': -123.06126402446726,
     ' N   ASN C 236  CA  ASN C 236  CB  ASN C 236  CG  ASN C 236 ': -156.76768556981202}

  reference_dihedral_proxies = reference_model.get_home_dihedral_proxies(
    work_params=work_params.reference_model,
    geometry=geometry,
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    geometry_ref=geometry_ref,
    sites_cart_ref=sites_cart_ref,
    pdb_hierarchy_ref=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy)
  assert reference_dihedral_proxies is not None
  assert len(reference_dihedral_proxies) == len(dihedral_hash)
  for rdp in reference_dihedral_proxies:
    assert rdp.limit == work_params.reference_model.limit

  r = rotalyze()
  rot_list_model, coot_model = r.analyze_pdb(
                                   hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  rot_list_reference, coot_reference = r.analyze_pdb(
                               hierarchy=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy)

  assert rot_list_model == """\
C 236 ASN:1.2:227.3:80.2:::t30
C 237 LEU:0.0:209.6:357.2:::OUTLIER"""

  assert rot_list_reference == """\
C 236 ASN:41.4:203.2:43.6:::t30
C 237 LEU:52.8:179.1:57.3:::tp"""

  reference_model.set_rotamer_to_reference(
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    pdb_hierarchy_ref=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy,
    xray_structure=xray_structure,
    quiet=True)
  rot_list_model, coot_model = r.analyze_pdb(
                                   hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  pdb_string = xray_structure.as_pdb_file()
  sites_cart = xray_structure.sites_cart()
  for atom in processed_pdb_file.all_chain_proxies.pdb_hierarchy.atoms():
    atom.set_xyz(sites_cart[atom.i_seq])
  rot_list_model, coot_model = r.analyze_pdb(
                                   hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  assert rot_list_model == """\
C 236 ASN:1.2:227.3:80.2:::t30
C 237 LEU:52.8:179.1:57.3:::tp"""

def run(args):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  exercise_reference_model(args, mon_lib_srv, ener_lib)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
