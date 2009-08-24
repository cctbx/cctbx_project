import mmtbx.monomer_library.server
import mmtbx.monomer_library.rotamer_utils
import iotbx.pdb.atom_name_interpretation
import iotbx.pdb.amino_acid_codes
import cctbx.geometry_restraints
import scitbx.rigid_body
import scitbx.graph.tardy_tree
from scitbx.array_family import flex
from scitbx import matrix
import libtbx.phil
from libtbx.test_utils import Exception_expected, show_diff
from libtbx.str_utils import show_string
from libtbx.utils import sequence_index_dict
import libtbx.load_env
import math
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

def report_tors(comp, residue_sites, matched_mon_lib_atom_names, targets):
  lookup = {}
  for j,atom_id in enumerate(matched_mon_lib_atom_names):
    lookup[atom_id] = j
  for tor in comp.tor_list:
    atom_ids = tor.atom_ids()
    js = [lookup.get(ai) for ai in atom_ids]
    if (js.count(None) != 0):
      angle_model = None
    else:
      d_sites = [residue_sites[j] for j in js]
      d = cctbx.geometry_restraints.dihedral(
        sites=d_sites, angle_ideal=0, weight=1)
      angle_model = d.angle_model
    target = targets.get(tor.id)
    if (angle_model is not None and target is not None):
      if (cctbx.geometry_restraints.angle_delta_deg(
            angle_1=angle_model,
            angle_2=target) > 1.e-5):
        annotation = "MISMATCH"
      else:
        annotation = "OK_target"
    else:
      annotation = "no_target"
    print tor.id, atom_ids, angle_model, annotation

def generate_rotamers(comp, rotamer_info, bonds_to_omit, strip_hydrogens):
  resname = comp.chem_comp.id
  comp_atom_names = set([atom.atom_id for atom in comp.atom_list])
  pdb_inp = iotbx.pdb.input(
    file_name=op.join(
      protein_pdb_files, reference_pdb_file_name_lookup[resname]))
  pdb_atoms = pdb_inp.atoms()
  pdb_atoms.reset_i_seq()
  matched_atom_names = iotbx.pdb.atom_name_interpretation.interpreters[
    resname].match_atom_names(atom_names=pdb_atoms.extract_name())
  names = matched_atom_names.unexpected
  if (len(names) != 0):
    raise RuntimeError("%s: unexpected atoms: %s" % (
      resname, " ".join(sorted(names))))
  names = matched_atom_names.missing_atom_names(ignore_hydrogen=True)
  if (len(names) != 0):
    raise RuntimeError("%s: missing atoms: %s" % (
      resname, " ".join(sorted(names))))
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  ag = pdb_hierarchy.only_atom_group()
  if (strip_hydrogens):
    for atom in ag.atoms():
      if (atom.element == " H"):
        ag.remove_atom(atom=atom)
  elif (    rotamer_info is not None
        and rotamer_info.atom_ids_not_handled is not None):
    pdb_residue = pdb_hierarchy.only_residue()
    atom_ids_not_handled = set(rotamer_info.atom_ids_not_handled)
    if (len(atom_ids_not_handled) != 0):
      for atom in ag.atoms():
        if (atom.name.strip() in atom_ids_not_handled):
          ag.remove_atom(atom=atom)
  pdb_residue = pdb_hierarchy.only_residue()
  pdb_atoms = pdb_residue.atoms()
  matched_mon_lib_atom_names = flex.select(
    sequence=matched_atom_names.mon_lib_names(),
    permutation=pdb_atoms.extract_i_seq())
  for p,m in zip(pdb_atoms, matched_mon_lib_atom_names):
    print 'atom name mapping: pdb="%s" -> %s' % (p.name, m)
  if (strip_hydrogens):
    assert len(matched_mon_lib_atom_names) == comp.chem_comp.number_atoms_nh
  comp_atom_name_set = set([atom.atom_id for atom in comp.atom_list])
  for name in matched_mon_lib_atom_names:
    if (name not in comp_atom_name_set):
      raise RuntimeError(
        "Missing comp atom: %s %s" % (pdb_residue.resname, name))
  #
  pdb_atoms.reset_i_seq()
  pdb_atoms.set_occ(new_occ=flex.double(pdb_atoms.size(), 1))
  pdb_atoms.set_b(new_b=flex.double(pdb_atoms.size(), 0))
  rg = pdb_hierarchy.only_residue_group()
  rg.resseq = 1
  rg.icode = " "
  assert pdb_hierarchy.only_atom_group().altloc == ""
  #
  if (rotamer_info is None):
    c = None
  else:
    c = rotamer_info.constrain_dihedrals_with_sigma_less_than_or_equal_to
  tardy_model = mmtbx.monomer_library.rotamer_utils.tardy_model(
    comp_comp_id=comp,
    input_atom_names=pdb_atoms.extract_name(),
    mon_lib_atom_names=matched_mon_lib_atom_names,
    sites_cart=pdb_atoms.extract_xyz(),
    bonds_to_omit=bonds_to_omit,
    constrain_dihedrals_with_sigma_less_than_or_equal_to=c)
  if (rotamer_info is None):
    return None
  #
  comp_tor_by_id = {}
  for tor in comp.tor_list:
    assert tor.id not in comp_tor_by_id
    comp_tor_by_id[tor.id] = tor
  rotmer_info_tor_ids = set(rotamer_info.tor_ids)
  rotamer_tor_by_id = {}
  for tor_atom_ids in rotamer_info.tor_atom_ids:
    assert len(tor_atom_ids) == 5
    tor_id = tor_atom_ids[0]
    assert tor_id in rotmer_info_tor_ids
    assert tor_id not in rotamer_tor_by_id
    rotamer_tor_by_id[tor_id] = tuple(tor_atom_ids[1:])
  rotamer_tor_atom_ids_by_tor_id = {}
  for tor_id in rotamer_info.tor_ids:
    atom_ids = rotamer_tor_by_id.get(tor_id)
    if (atom_ids is not None):
      rotamer_tor_atom_ids_by_tor_id[tor_id] = atom_ids
    else:
      comp_tor = comp_tor_by_id.get(tor_id)
      if (comp_tor is not None):
        rotamer_tor_atom_ids_by_tor_id[tor_id] = comp_tor.atom_ids()
      else:
        raise RuntimeError(
          "rotamer_info.tor_id %s is unknown." % show_string(tor_id))
  #
  tor_id_by_rotatable_bond_atom_names = {}
  for tor_id,atom_ids in rotamer_tor_atom_ids_by_tor_id.items():
    atom_names = tuple(sorted(atom_ids[1:3]))
    assert atom_names not in tor_id_by_rotatable_bond_atom_names
    tor_id_by_rotatable_bond_atom_names[atom_names] = tor_id
  #
  tor_id_i_q_packed_matches = {}
  number_of_trees = 0
  for i_body,he in enumerate(
                     tardy_model.tardy_tree.cluster_manager.hinge_edges):
    if (he[0] == -1):
      number_of_trees += 1
      continue
    hinge_atom_names = [tardy_model.labels[i].strip() for i in he]
    atom_names = tuple(sorted(hinge_atom_names))
    tor_id = tor_id_by_rotatable_bond_atom_names.get(atom_names)
    if (tor_id is None):
      raise RuntimeError(
        "rotatable bond atoms %s - %s (as defined by tardy_tree):"
        " no match in rotamer_info.tor_ids" % tuple(hinge_atom_names))
    tor_id_i_q_packed_matches[tor_id] = i_body - 1
  assert number_of_trees == 1
  #
  unused_rotamer_info_tor_ids = []
  for tor_id in rotamer_info.tor_ids:
    if (tor_id not in tor_id_i_q_packed_matches):
      unused_rotamer_info_tor_ids.append(tor_id)
  if (len(unused_rotamer_info_tor_ids) != 0):
    print "Info: unused rotamer_info.tor_ids:", \
      " ".join(unused_rotamer_info_tor_ids)
    assert strip_hydrogens
  #
  tors_start = {}
  atom_indices = sequence_index_dict(seq=matched_mon_lib_atom_names)
  for tor_id in tor_id_i_q_packed_matches.keys():
    tor_atom_ids = rotamer_tor_atom_ids_by_tor_id[tor_id]
    ai = [atom_indices.get(atom_id) for atom_id in tor_atom_ids]
    assert ai.count(None) == 0
    d_sites = [tardy_model.sites[i] for i in ai]
    d = cctbx.geometry_restraints.dihedral(
      sites=d_sites, angle_ideal=0, weight=1)
    assert tor_id not in tors_start
    tors_start[tor_id] = d.angle_model
  #
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
  for i_rotamer,rotamer in enumerate(rotamer_info.rotamer):
    q_packed_work = flex.double(tardy_model.q_packed_size, 0)
    for tor_id,angle in zip(rotamer_info.tor_ids, rotamer.angles):
      i_q_packed = tor_id_i_q_packed_matches.get(tor_id)
      if (i_q_packed is not None and angle is not None):
        q_packed_work[i_q_packed] = math.radians(angle - tors_start[tor_id])
    tardy_model.unpack_q(q_packed=q_packed_work)
    rotamer_sites = tardy_model.sites_moved()
    rotamer_sites += matrix.col((4,4,4)) * i_rotamer
    pdb_atoms.set_xyz(new_xyz=rotamer_sites)
    report_tors(
      comp=comp,
      residue_sites=rotamer_sites,
      matched_mon_lib_atom_names=matched_mon_lib_atom_names,
      targets=dict(zip(rotamer_info.tor_ids, rotamer.angles)))
    #
    file_name = "%s/%s_%s.pdb" % (
      rotamers_sep_sub_dir, pdb_residue.resname, rotamer.id)
    print "Writing file:", file_name
    f = open(file_name, "w")
    print >> f, "REMARK %s %s" % (pdb_residue.resname, rotamer.id)
    pdb_atoms.reset_serial(first_value=1)
    f.write(pdb_hierarchy.as_pdb_string(append_end=True))
    del f
    #
    pdb_atoms.reset_serial(first_value=atom_serial_first_value)
    atom_serial_first_value += pdb_atoms.size()
    chain_id = (string.uppercase + string.lowercase)[i_rotamer]
    pdb_hierarchy.only_chain().id = chain_id
    remark_strings.append(
      "REMARK %s %s = chain %s" % (pdb_residue.resname, rotamer.id, chain_id))
    atom_strings.append(pdb_hierarchy.as_pdb_string(append_end=False))
  file_name = "%s/%s.pdb" % (rotamers_sub_dir, pdb_residue.resname)
  print "Writing file:", file_name
  f = open(file_name, "w")
  for s in remark_strings:
    print >> f, s
  for s in atom_strings:
    f.write(s)
  print >> f, "END"
  del f

def process_rotamer_info(rotamer_info_master_phil, comp):
  assert len(comp.rotamer_info) < 2
  if (len(comp.rotamer_info) == 0):
    print "No rotamer_info."
    rotamer_info = None
  else:
    rotamer_info_phil = rotamer_info_master_phil.fetch(
      source=libtbx.phil.parse(
        input_string=comp.rotamer_info[0].phil_str))
    rotamer_info_phil.show()
    rotamer_info = rotamer_info_phil.extract()
    n_missing_frequencies = 0
    for rotamer in rotamer_info.rotamer:
      assert rotamer.id is not None
      assert len(rotamer.id.strip()) == len(rotamer.id)
      assert len(rotamer.id.split()) == 1
      if (rotamer.frequency is None):
        if (rotamer.frequency_annotation != "for more uniform sampling"):
          n_missing_frequencies += 1
      else:
        assert rotamer.frequency > 0
        assert rotamer.frequency < 1
      assert rotamer.angles is not None
      assert len(rotamer.angles) == len(rotamer_info.tor_ids)
      for angle in rotamer.angles:
        assert angle is None or -180 < angle <= 180
    if (n_missing_frequencies != 0):
      print "Warning: number of missing frequencies:", n_missing_frequencies
  return rotamer_info

def exercise_server_rotamer_iterator(mon_lib_srv, comp_comp_id):
  resname = comp_comp_id.chem_comp.id
  pdb_inp = iotbx.pdb.input(
    file_name=op.join(
      protein_pdb_files, reference_pdb_file_name_lookup[resname]))
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  atom_ids_not_handled = {
    "ASP": ["HD2"],
    "GLU": ["HE2"]}.get(comp_comp_id.chem_comp.id)
  pdb_atoms = pdb_hierarchy.only_residue().atoms()
  if (atom_ids_not_handled is not None):
    try:
      mon_lib_srv.rotamer_iterator(
        comp_comp_id=comp_comp_id,
        atom_names=pdb_atoms.extract_name(),
        sites_cart=pdb_atoms.extract_xyz(),
        skip_if_atom_name_problems=False)
    except RuntimeError, e:
      assert not show_diff(str(e)[3:-3],
        ": rotamer_info does not handle these atoms: ")
    else: raise Exception_expected
    rotamer_iterator = mon_lib_srv.rotamer_iterator(
      comp_comp_id=comp_comp_id,
      atom_names=pdb_atoms.extract_name(),
      sites_cart=pdb_atoms.extract_xyz(),
      skip_if_atom_name_problems=True)
    assert rotamer_iterator is None
    ag = pdb_hierarchy.only_atom_group()
    for atom in ag.atoms():
      if (atom.name.strip() in atom_ids_not_handled):
        ag.remove_atom(atom=atom)
    pdb_atoms = pdb_hierarchy.only_residue().atoms()
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    comp_comp_id=comp_comp_id,
    atom_names=pdb_atoms.extract_name(),
    sites_cart=pdb_atoms.extract_xyz(),
    skip_if_atom_name_problems=False)
  if (rotamer_iterator is not None):
    for obj in rotamer_iterator():
      print obj

def process(mon_lib_srv, resname):
  print "resname:", resname
  comp = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
  tor_dict = {}
  for tor in comp.tor_list:
    atom_names = tuple(sorted([tor.atom_id_2, tor.atom_id_3]))
    tor_dict.setdefault(atom_names, []).append(tor)
    print tor.id, ", ".join([
      tor.atom_id_1, tor.atom_id_2, tor.atom_id_3, tor.atom_id_4]), \
      tor.value_angle_esd
  for atom_ids,tors in tor_dict.items():
    if (len(tors) != 1):
      print "Info: redundant tors:", ", ".join([tor.id for tor in tors])
  #
  rotamer_info = process_rotamer_info(
    rotamer_info_master_phil=mon_lib_srv.rotamer_info_master_phil(),
    comp=comp)
  bonds_to_omit = mmtbx.monomer_library.rotamer_utils.extract_bonds_to_omit(
    rotamer_info=rotamer_info)
  for strip_hydrogens in [True, False]:
    generate_rotamers(
      comp=comp,
      rotamer_info=rotamer_info,
      bonds_to_omit=bonds_to_omit,
      strip_hydrogens=strip_hydrogens)
  if (len(comp.rotamer_info) != 0):
    exercise_server_rotamer_iterator(
      mon_lib_srv=mon_lib_srv, comp_comp_id=comp)
  print

def run(args):
  assert len(args) == 0
  mon_lib_srv = mmtbx.monomer_library.server.server()
  amino_acid_resnames = sorted(
    iotbx.pdb.amino_acid_codes.one_letter_given_three_letter.keys())
  for resname in amino_acid_resnames:
    process(
      mon_lib_srv=mon_lib_srv,
      resname=resname)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
