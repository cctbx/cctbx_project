import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from mmtbx.refinement.reference_model import reference_model
from mmtbx.validation.rotalyze import rotalyze
import iotbx.phil
import iotbx.utils
from iotbx import file_reader
import libtbx.load_env
import sys, os

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
""".splitlines()

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
""".splitlines()

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
reference_model
{
    file = None
      .type = path
      .short_caption = PDB file
      .style = bold file_type:pdb
    selection = None
      .type = str
      .short_caption = Atom selection
      .input_size = 400
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
    fix_outliers = True
      .type = bool
    auto_align = False
      .type = bool
    secondary_structure_only = False
      .type = bool
    reference_group
      .multiple=True
      .optional=True
      .short_caption=Reference group
      .style = noauto auto_align
    {
      reference= None
        .type=str
        .optional=True
        .short_caption=Reference selection
        .input_size=400
        .style = selection bold
      selection= None
        .type=str
        .short_caption=Restrained selection
        .input_size=400
        .style = selection bold
    }
    alignment
      .help = Set of parameters for sequence alignment. Defaults are good for most \
          of cases
      .short_caption = Sequence alignment
      .style = box auto_align
    {
      alignment_style =  local *global
        .type = choice
      gap_opening_penalty = 1
        .type = float
      gap_extension_penalty = 1
        .type = float
      similarity_matrix =  blosum50  dayhoff *identity
        .type = choice
    }
    alignment_group
     .multiple=True
     .optional=True
     .short_caption=Reference group
     .style = noauto auto_align menu_item parent_submenu:reference_model
    {
       reference=None
         .type=str
         .short_caption=Reference selection
         .style = selection
       selection=None
         .type=str
         .short_caption=Restrained selection
         .style = selection
    }
}
""")

def exercise_reference_model(args, mon_lib_srv, ener_lib):
  rm = reference_model()
  master_phil = get_master_phil()
  #print master_phil.show()
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
  pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy
  processed_pdb_file_ref = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=reference_raw_records,
    for_dihedral_reference=True)
  processed_pdb_file_ref_alt_seq = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=reference_raw_records_alt_seq,
    for_dihedral_reference=True)
  processed_pdb_file_ref_match = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=reference_raw_records_match,
    for_dihedral_reference=True)
  pdb_hierarchy_ref_alt_seq = processed_pdb_file_ref_alt_seq.all_chain_proxies.pdb_hierarchy
  pdb_hierarchy_ref_match = processed_pdb_file_ref_match.all_chain_proxies.pdb_hierarchy
  geometry_ref = processed_pdb_file_ref.geometry_restraints_manager()
  sites_cart_ref = processed_pdb_file_ref.all_chain_proxies.sites_cart
  pdb_hierarchy_ref=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy
  i_seq_name_hash = rm.build_name_hash(
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  assert i_seq_name_hash == \
    {0: ' N   ASN C 236 ', 1: ' CA  ASN C 236 ', 2: ' C   ASN C 236 ',
     3: ' O   ASN C 236 ', 4: ' CB  ASN C 236 ', 5: ' CG  ASN C 236 ',
     6: ' OD1 ASN C 236 ', 7: ' ND2 ASN C 236 ', 8: ' N   LEU C 237 ',
     9: ' CA  LEU C 237 ', 10: ' C   LEU C 237 ', 11: ' O   LEU C 237 ',
     12: ' CB  LEU C 237 ', 13: ' CG  LEU C 237 ', 14: ' CD1 LEU C 237 ',
     15: ' CD2 LEU C 237 '}

  i_seq_element_hash = rm.build_element_hash(
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  assert i_seq_element_hash == \
    {0: ' N', 1: ' C', 2: ' C', 3: ' O', 4: ' C', 5: ' C', 6: ' O', 7: ' N',
     8: ' N', 9: ' C', 10: ' C', 11: ' O', 12: ' C', 13: ' C', 14: ' C',
     15: ' C'}

  dihedral_hash = rm.build_dihedral_hash(
    geometry=geometry_ref,
    sites_cart=sites_cart_ref,
    pdb_hierarchy=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy,
    include_hydrogens=False,
    include_main_chain=True,
    include_side_chain=True)
  assert len(dihedral_hash) == 9

  reference_dihedral_proxies = rm.get_home_dihedral_proxies(
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
  master_phil_str_overrides = """
  reference_model.reference_group {
    reference= chain B and resseq 246:247
    selection= chain C and resid 236:237
  }
  """
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  work_params_alt = master_phil.fetch(sources=phil_objects).extract()
  match_map = rm.process_reference_groups(
                pdb_hierarchy,
                pdb_hierarchy_ref_alt_seq,
                work_params_alt.reference_model)
  assert match_map == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11,
   12: 12, 13: 13, 14: 14, 15: 15}
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

  rm.set_rotamer_to_reference(
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

  cbetadev_hash = rm.build_cbetadev_hash(
                    pdb_hierarchy=processed_pdb_file_ref.all_chain_proxies.pdb_hierarchy)
  assert cbetadev_hash == \
    {' ASN C 236': '  0.015', ' LEU C 237': '  0.038'}
  match_map = rm.process_reference_groups(pdb_hierarchy=pdb_hierarchy,
                                                       pdb_hierarchy_ref=pdb_hierarchy_ref,
                                                       params=work_params.reference_model)
  assert match_map == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11,
   12: 12, 13: 13, 14: 14, 15: 15}
  master_phil_str_overrides = """
  reference_model {
    auto_align = True
  }
  """
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  work_params_match = master_phil.fetch(sources=phil_objects).extract()
  match_map = rm.process_reference_groups(pdb_hierarchy=pdb_hierarchy,
                                                       pdb_hierarchy_ref=pdb_hierarchy_ref_match,
                                                       params=work_params_match.reference_model)
  assert match_map == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}

  master_phil_str_overrides = """
  reference_model {
    auto_align = True
    alignment_group {
      reference = resid 270
      selection = resid 236
    }
  }
  """
  phil_objects = [
    iotbx.phil.parse(input_string=master_phil_str_overrides)]
  work_params_match = master_phil.fetch(sources=phil_objects).extract()
  match_map = rm.process_reference_groups(pdb_hierarchy=pdb_hierarchy,
                                                       pdb_hierarchy_ref=pdb_hierarchy_ref_match,
                                                       params=work_params_match.reference_model)
  assert match_map == \
  {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}

  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb").file_object
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    pdb_inp = pdb_in,
    for_dihedral_reference=True)
  geometry = processed_pdb_file.geometry_restraints_manager()
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  xray_structure=processed_pdb_file.xray_structure()
  pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy
  reference_dihedral_proxies = rm.get_home_dihedral_proxies(
    work_params=work_params.reference_model,
    geometry=geometry,
    pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    geometry_ref=geometry,
    sites_cart_ref=sites_cart,
    pdb_hierarchy_ref=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  standard_weight = 0
  rm.show_reference_summary()
  for dp in reference_dihedral_proxies:
    if dp.weight == 1.0:
      standard_weight += 1
  assert standard_weight == 1403
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
    reference_dihedral_proxies = rm.get_home_dihedral_proxies(
      work_params=work_params_ss.reference_model,
      geometry=geometry,
      pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      geometry_ref=geometry,
      sites_cart_ref=sites_cart,
      pdb_hierarchy_ref=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
    ss_weight = 0
    for dp in reference_dihedral_proxies:
      if dp.weight == 1.0:
        ss_weight += 1
    assert ss_weight == 837

def run(args):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  exercise_reference_model(args, mon_lib_srv, ener_lib)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
