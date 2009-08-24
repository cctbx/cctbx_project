import mmtbx.monomer_library.server
import iotbx.pdb.amino_acid_codes
import cctbx.geometry_restraints
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import show_diff
from libtbx.utils import sequence_index_dict
import libtbx.load_env
import string
import sys, os
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
      d = cctbx.geometry_restraints.dihedral(
        sites=d_sites, angle_ideal=0, weight=1)
      angle_model = d.angle_model
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
    if (verbose): print tor_id, atom_ids, angle_model, annotation
  assert n_mismatches == 0, n_mismatches

def exercise_server_rotamer_iterator(mon_lib_srv, resname, verbose):
  pdb_inp = iotbx.pdb.input(
    file_name=op.join(
      protein_pdb_files, reference_pdb_file_name_lookup[resname]))
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  atom_ids_not_handled = {
    "ASP": ["HD2"],
    "GLU": ["HE2"]}.get(resname)
  pdb_atoms = pdb_hierarchy.only_residue().atoms()
  comp_comp_id = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
  if (atom_ids_not_handled is not None):
    rotamer_iterator = mon_lib_srv.rotamer_iterator(
      comp_comp_id=comp_comp_id,
      atom_names=pdb_atoms.extract_name(),
      sites_cart=pdb_atoms.extract_xyz())
    assert not show_diff(str(rotamer_iterator.problem_message)[3:-3],
      ": rotamer_info does not handle these atoms: ")
    ag = pdb_hierarchy.only_atom_group()
    for atom in ag.atoms():
      if (atom.name.strip() in atom_ids_not_handled):
        ag.remove_atom(atom=atom)
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
    rotamer_iterator = mon_lib_srv.rotamer_iterator(
      comp_comp_id=comp_comp_id,
      atom_names=pdb_atoms.extract_name(),
      sites_cart=pdb_atoms.extract_xyz())
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
          print "Info: unused rotamer_info.tor_ids:", \
            " ".join(unused_rotamer_info_tor_ids)
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
        if (verbose): print "Writing file:", file_name
        f = open(file_name, "w")
        print >> f, "REMARK %s %s" % (resname, rotamer.id)
        f.write(pdb_hierarchy_work.as_pdb_string(append_end=True))
        del f
        #
        rotamer_sites_cart += matrix.col((4,4,4)) * i_rotamer
        pdb_atoms_work.set_xyz(new_xyz=rotamer_sites_cart)
        pdb_atoms_work.reset_serial(first_value=atom_serial_first_value)
        atom_serial_first_value += pdb_atoms_work.size()
        chain_id = (string.uppercase + string.lowercase)[i_rotamer]
        pdb_hierarchy_work.only_chain().id = chain_id
        remark_strings.append(
          "REMARK %s %s = chain %s" % (resname, rotamer.id, chain_id))
        atom_strings.append(pdb_hierarchy_work.as_pdb_string(append_end=False))
      file_name = "%s/%s.pdb" % (rotamers_sub_dir, resname)
      if (verbose): print "Writing file:", file_name
      f = open(file_name, "w")
      for s in remark_strings:
        print >> f, s
      for s in atom_strings:
        f.write(s)
      print >> f, "END"
      del f

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = "--verbose" in args
  mon_lib_srv = mmtbx.monomer_library.server.server()
  amino_acid_resnames = sorted(
    iotbx.pdb.amino_acid_codes.one_letter_given_three_letter.keys())
  for resname in amino_acid_resnames:
    if (verbose): print "resname:", resname
    exercise_server_rotamer_iterator(
      mon_lib_srv=mon_lib_srv,
      resname=resname,
      verbose=verbose)
    if (verbose): print
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
