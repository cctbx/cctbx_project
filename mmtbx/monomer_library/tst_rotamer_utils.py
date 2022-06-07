from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.server
import iotbx.pdb.amino_acid_codes
import cctbx.geometry_restraints
import scitbx.math
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import show_diff
from libtbx.utils import sequence_index_dict, format_cpu_times
import libtbx.load_env
import string
import sys, os
from six.moves import zip
from six.moves import range
op = os.path

protein_pdb_files = libtbx.env.find_in_repositories(
  relative_path="phenix_regression/protein_pdb_files",
  optional=False)

reference_pdb_file_names = """\
ala_chain_all_h_1ozo_v3.ent
arg_chain_all_h_1o8t_v3.ent
asn_chain_all_h_1o8t_v3.ent
asp_chain_all_h_1jjx_v3.ent
cys_chain_all_h_1rfa_v3.ent
gln_chain_all_h_1o8t_v3.ent
glu_chain_all_h_1bm4_v3.ent
gly_chain_all_h_1ozo_v3.ent
his_chain_all_h_1g7e_v3.ent
ile_chain_all_h_1ozo_v3.ent
leu_chain_all_h_1ozo_v3.ent
lys_chain_all_h_1o8t_v3.ent
met_chain_all_h_1ozo_v3.ent
mse_chain_all_h_1ozo_v3.ent
phe_chain_all_h_1hdj_v3.ent
pro_chain_all_h_1a03_v3.ent
ser_chain_all_h_1o8t_v3.ent
thr_chain_all_h_1o8t_v3.ent
trp_chain_all_h_1cx1_v3.ent
tyr_chain_all_h_1cx1_v3.ent
val_chain_all_h_1ozo_v3.ent
""".splitlines()
def __init_reference_pdb_file_name_lookup():
  result = {}
  for file_name in reference_pdb_file_names:
    result[file_name[:3].upper()] = file_name
  return result
reference_pdb_file_name_lookup = __init_reference_pdb_file_name_lookup()

def report_tors(
      mon_lib_atom_names,
      sites_cart,
      tor_atom_ids_by_tor_id,
      target_angles,
      verbose):
  n_mismatches = 0
  i_seq_by_atom_name = sequence_index_dict(seq=mon_lib_atom_names)
  for tor_id,atom_ids in tor_atom_ids_by_tor_id.items():
    js = [i_seq_by_atom_name.get(ai) for ai in atom_ids]
    if (js.count(None) != 0):
      angle_model = None
    else:
      d_sites = [sites_cart[j] for j in js]
      angle_model = scitbx.math.dihedral_angle(sites=d_sites, deg=True)
    target_angle = target_angles.get(tor_id)
    if (angle_model is not None and target_angle is not None):
      if (cctbx.geometry_restraints.angle_delta_deg(
            angle_1=angle_model,
            angle_2=target_angle) > 1.e-5):
        annotation = "MISMATCH"
        n_mismatches += 1
      else:
        annotation = "OK_target"
    else:
      annotation = "no_target"
    if (verbose): print(tor_id, atom_ids, angle_model, annotation)
  assert n_mismatches == 0, n_mismatches

def remove_atom_ids_not_handled(pdb_hierarchy, atom_ids_not_handled):
  if (atom_ids_not_handled is None): return
  ag = pdb_hierarchy.only_atom_group()
  for atom in ag.atoms():
    if (atom.name.strip() in atom_ids_not_handled):
      ag.remove_atom(atom=atom)

def exercise_server_rotamer_iterator(mon_lib_srv, pdb_hierarchy, verbose):
  resname = pdb_hierarchy.only_residue().resname
  atom_ids_not_handled = {
    "ASP": ["HD2"],
    "GLU": ["HE2"]}.get(resname)
  pdb_atoms = pdb_hierarchy.only_residue().atoms()
  comp_comp_id = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
  if (atom_ids_not_handled is not None):
    rotamer_iterator = comp_comp_id.rotamer_iterator(
      atom_names=pdb_atoms.extract_name(),
      sites_cart=pdb_atoms.extract_xyz())
    assert not show_diff(str(rotamer_iterator.problem_message)[3:-3],
      ": rotamer_info does not handle these atoms: ")
    remove_atom_ids_not_handled(
      pdb_hierarchy=pdb_hierarchy,
      atom_ids_not_handled=atom_ids_not_handled)
    pdb_atoms = pdb_hierarchy.only_residue().atoms()
  #
  pdb_atoms.set_occ(new_occ=flex.double(pdb_atoms.size(), 1))
  pdb_atoms.set_b(new_b=flex.double(pdb_atoms.size(), 0))
  rg = pdb_hierarchy.only_residue_group()
  rg.resseq = 1
  rg.icode = " "
  #
  assert pdb_hierarchy.only_atom_group().altloc == ""
  for strip_hydrogens in [False, True]:
    if (strip_hydrogens):
      ag = pdb_hierarchy.only_atom_group()
      for atom in ag.atoms():
        if (atom.element == " H"):
          ag.remove_atom(atom=atom)
      pdb_atoms = pdb_hierarchy.only_residue().atoms()
    pdb_atoms.reset_i_seq()
    rotamer_iterator = comp_comp_id.rotamer_iterator(
      atom_names=pdb_atoms.extract_name(),
      sites_cart=pdb_atoms.extract_xyz(),fine_sampling=False)
    if (rotamer_iterator.problem_message):
      raise RuntimeError(rotamer_iterator.problem_message)
    if (rotamer_iterator.rotamer_info is not None):
      if (strip_hydrogens):
        assert len(rotamer_iterator.mon_lib_atom_names) \
            == comp_comp_id.chem_comp.number_atoms_nh
      comp_atom_name_set = set([atom.atom_id
        for atom in comp_comp_id.atom_list])
      for name in rotamer_iterator.mon_lib_atom_names:
        if (name not in comp_atom_name_set):
          raise RuntimeError("Missing comp atom: %s %s" % (resname, name))
      unused_rotamer_info_tor_ids = []
      for tor_id in rotamer_iterator.rotamer_info.tor_ids:
        if (tor_id not in rotamer_iterator.i_q_packed_by_tor_id):
          unused_rotamer_info_tor_ids.append(tor_id)
      if (len(unused_rotamer_info_tor_ids) != 0):
        if (verbose):
          print("Info: unused rotamer_info.tor_ids:", \
            " ".join(unused_rotamer_info_tor_ids))
        assert strip_hydrogens
      if (strip_hydrogens):
        rotamers_sub_dir = "rotamers_no_h"
      else:
        rotamers_sub_dir = "rotamers_with_h"
      rotamers_sep_sub_dir = rotamers_sub_dir + "_sep"
      if (not os.path.isdir(rotamers_sub_dir)):
        os.mkdir(rotamers_sub_dir)
      if (not os.path.isdir(rotamers_sep_sub_dir)):
        os.mkdir(rotamers_sep_sub_dir)
      remark_strings = []
      atom_strings = []
      atom_serial_first_value = 1
      for i_rotamer,(rotamer,rotamer_sites_cart) in enumerate(rotamer_iterator):
        report_tors(
          mon_lib_atom_names=rotamer_iterator.mon_lib_atom_names,
          sites_cart=rotamer_sites_cart,
          tor_atom_ids_by_tor_id
            =rotamer_iterator.rotamer_tor_atom_ids_by_tor_id,
          target_angles=dict(zip(
            rotamer_iterator.rotamer_info.tor_ids,
            rotamer.angles)),
          verbose=verbose)
        #
        pdb_hierarchy_work = pdb_hierarchy.deep_copy()
        pdb_atoms_work = pdb_hierarchy_work.atoms()
        pdb_atoms_work.set_xyz(new_xyz=rotamer_sites_cart)
        pdb_atoms_work.reset_serial(first_value=1)
        file_name = "%s/%s_%s.pdb" % (
          rotamers_sep_sub_dir, resname, rotamer.id)
        if (verbose): print("Writing file:", file_name)
        with  open(file_name, "w") as f:
          print("REMARK %s %s" % (resname, rotamer.id), file=f)
          f.write(pdb_hierarchy_work.as_pdb_string(append_end=True))
        #
        rotamer_sites_cart += matrix.col((4,4,4)) * i_rotamer
        pdb_atoms_work.set_xyz(new_xyz=rotamer_sites_cart)
        pdb_atoms_work.reset_serial(first_value=atom_serial_first_value)
        atom_serial_first_value += pdb_atoms_work.size()
        nl = ''.join(['%i' % i for i in range(10)])
        chain_id = (string.ascii_uppercase + string.ascii_lowercase + nl)[i_rotamer]
        pdb_hierarchy_work.only_chain().id = chain_id
        remark_strings.append(
          "REMARK %s %s = chain %s" % (resname, rotamer.id, chain_id))
        atom_strings.append(pdb_hierarchy_work.as_pdb_string(append_end=False))
      file_name = "%s/%s.pdb" % (rotamers_sub_dir, resname)
      if (verbose): print("Writing file:", file_name)
      with open(file_name, "w") as f:
        for s in remark_strings:
          print(s, file=f)
        for s in atom_strings:
          f.write(s)
        print("END", file=f)

def compare_dihedrals(
      mon_lib_srv,
      amino_acid_resnames,
      pdb_dir,
      file_name_extension,
      verbose):
  for resname in amino_acid_resnames:
    if (resname in ["PRO", 'PYL', 'SEC']):
      # compatible semi emp files not available (not important enough to
      # warrant extra effort)
      continue
    rotamer_info = mon_lib_srv.get_comp_comp_id_direct(
      comp_id=resname).rotamer_info()
    if (rotamer_info is None): continue
    for rotamer in rotamer_info.rotamer:
      file_name = "%s_%s%s" % (resname, rotamer.id, file_name_extension)
      if (verbose): print(file_name)
      path = op.join(pdb_dir, file_name)
      pdb_inp = iotbx.pdb.input(file_name=path)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      ag = pdb_hierarchy.only_atom_group()
      for atom in ag.atoms():
        if (   atom.name in [" H2 ",  " HC1"]
            or (resname == "ASP" and atom.name == "HD21")
            or (resname == "GLU" and atom.name == "HE21")):
          ag.remove_atom(atom=atom)
      comp_comp_id = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
      pdb_atoms = pdb_hierarchy.only_residue().atoms()
      rotamer_iterator = comp_comp_id.rotamer_iterator(
        atom_names=pdb_atoms.extract_name(),
        sites_cart=pdb_atoms.extract_xyz())
      assert rotamer_iterator.problem_message is None
      angle_start_by_tor_id = rotamer_iterator.angle_start_by_tor_id
      for tor_id,angle_tab in zip(rotamer_info.tor_ids, rotamer.angles):
        angle_pdb = angle_start_by_tor_id[tor_id]
        if (verbose): print(tor_id, angle_tab, angle_pdb)
        if (cctbx.geometry_restraints.angle_delta_deg(
              angle_1=angle_tab,
              angle_2=angle_pdb) > 0.5):
          if (verbose
              or resname not in ["ARG", "ASN", "GLN"]):
            # Keeping all hydrogen dihedrals in ARG, ASN, GLN at 180
            # after discussions.
            print("Mismatch", resname, rotamer.id, tor_id, \
              "pdb: %.0f" % angle_pdb, \
              "tab: %.0f" % angle_tab)
      if (verbose): print()

def exercise_termini(mon_lib_srv, pdb_file_name):
  pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  resname = pdb_hierarchy.only_residue().resname
  atom_ids_not_handled_by_resname = {
    "ASP": ["HD1", "HD2", "2HD"],
    "GLU": ["HE1", "HE2", "2HE"]}
  remove_atom_ids_not_handled(
    pdb_hierarchy=pdb_hierarchy,
    atom_ids_not_handled=atom_ids_not_handled_by_resname.get(resname))
  pdb_atoms = pdb_hierarchy.only_residue().atoms()
  comp_comp_id = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
  rotamer_iterator = comp_comp_id.rotamer_iterator(
    atom_names=pdb_atoms.extract_name(),
    sites_cart=pdb_atoms.extract_xyz())
  msg = rotamer_iterator.problem_message
  if (msg is not None):
    raise RuntimeError("rotamer_iterator.problem_message: %s" % msg)

def exercise_pro_missing_hd1(mon_lib_srv):
  pdb_inp = iotbx.pdb.input(source_info=None, lines="""\
ATOM    110  N   PRO A 263       0.453 -20.680 -39.256  1.00 53.34           N
ATOM    111  CA  PRO A 263       0.444 -22.054 -39.751  1.00 50.42           C
ATOM    112  C   PRO A 263       0.860 -22.998 -38.645  1.00 52.10           C
ATOM    113  O   PRO A 263       1.693 -22.614 -37.817  1.00 48.32           O
ATOM    114  CB  PRO A 263       1.491 -22.052 -40.887  1.00 53.30           C
ATOM    115  CG  PRO A 263       2.012 -20.645 -40.990  1.00 57.05           C
ATOM    116  CD  PRO A 263       1.586 -19.897 -39.782  1.00 53.45           C
ATOM    117  HA  PRO A 263      -0.437 -22.302 -40.100  1.00 60.51           H
ATOM    118  HB2 PRO A 263       2.210 -22.664 -40.664  1.00 63.96           H
ATOM    119  HB3 PRO A 263       1.066 -22.318 -41.718  1.00 63.96           H
ATOM    120  HG2 PRO A 263       2.980 -20.669 -41.043  1.00 68.45           H
ATOM    121  HG3 PRO A 263       1.645 -20.229 -41.786  1.00 68.45           H
ATOM    122  HD2 PRO A 263       1.267 -19.021 -40.049  1.00 64.14           H
""")
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  residue = pdb_hierarchy.only_residue()
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    comp_id=residue.resname,
    atom_names=residue.atoms().extract_name(),
    sites_cart=residue.atoms().extract_xyz())
  assert not show_diff(
    rotamer_iterator.problem_message,
    'resname=PRO: missing atom "HD1" for tor_id "hh3"')

def run(args):
  verbose = False
  semi_emp_rotamer_pdb_dirs = []
  for arg in args:
    if (arg == "--verbose"):
      verbose = True
    else:
      assert op.isdir(arg)
      semi_emp_rotamer_pdb_dirs.append(arg)
  mon_lib_srv = mmtbx.monomer_library.server.server()
  amino_acid_resnames = sorted(
    iotbx.pdb.amino_acid_codes.one_letter_given_three_letter.keys())
  for resname in amino_acid_resnames:
    if (verbose): print("resname:", resname)
    if resname in ["UNK", 'PYL', 'SEC']:
      # skipping UNK residue because there is no rotamers available for it
      continue
    pdb_inp = iotbx.pdb.input(
      file_name=op.join(
        protein_pdb_files, reference_pdb_file_name_lookup[resname]))
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    exercise_server_rotamer_iterator(
      mon_lib_srv=mon_lib_srv,
      pdb_hierarchy=pdb_hierarchy,
      verbose=verbose)
    if (verbose): print()
  if (len(semi_emp_rotamer_pdb_dirs) == 0):
    pdb_dir = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/semi_emp_rotamer_pdb")
    if (pdb_dir is None):
      print("Skipping compare_dihedrals(): semi_emp_rotamer_pdb not available.")
    else:
      semi_emp_rotamer_pdb_dirs.append(pdb_dir)
  for pdb_dir in semi_emp_rotamer_pdb_dirs:
    compare_dihedrals(
      mon_lib_srv=mon_lib_srv,
      amino_acid_resnames=amino_acid_resnames,
      pdb_dir=pdb_dir,
      file_name_extension=".uhf_631dp.pdb",
      verbose=verbose)
  for file_name in os.listdir(protein_pdb_files):
    if (not file_name.endswith(".ent")): continue
    if (verbose): print(file_name)
    exercise_termini(
      mon_lib_srv=mon_lib_srv,
      pdb_file_name=op.join(protein_pdb_files, file_name))
    if (verbose): print()
  exercise_pro_missing_hd1(mon_lib_srv=mon_lib_srv)
  print(format_cpu_times())
  print("OK")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
