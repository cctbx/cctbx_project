import mmtbx.alignment
from iotbx.pdb import amino_acid_codes
from libtbx import group_args
import cctbx.geometry_restraints
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.refinement import fit_rotamers
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
from libtbx import Auto
import libtbx.load_env
from libtbx.utils import format_exception, Sorry
from mmtbx import secondary_structure
import sys, re

reference_group_params = iotbx.phil.parse("""
 reference_group
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
""")

class reference_model(object):

  def selection(self, string, cache):
    return cache.selection(
      string=string)

  def iselection(self, string, cache=None):
    return self.selection(string=string, cache=cache).iselection()

  def phil_atom_selection_multiple(
        self,
        cache,
        string_list,
        allow_none=False,
        allow_auto=False,
        raise_if_empty_selection=True):
    result = []
    for string in string_list:
      if (string is None):
        if (allow_none): return None
        raise Sorry('Atom selection cannot be None:\n  =None')
      elif (string is Auto):
        if (allow_auto): return Auto
        raise Sorry('Atom selection cannot be Auto:\n  %s=Auto')
      try:
          result.append(self.selection(string=string, cache=cache).iselection())
      except KeyboardInterrupt: raise
      except Exception, e: # keep e alive to avoid traceback
        fe = format_exception()
        raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
          'reference_group', string, fe))
      if (raise_if_empty_selection and result.count(True) == 0):
        raise Sorry('Empty atom selection:\n  %s=%s' % (
          'reference_group', string))
    return result

  def phil_atom_selections_as_i_seqs_multiple(self,
                                              cache,
                                              string_list):
    result = []
    iselection = self.phil_atom_selection_multiple(
          cache=cache,
          string_list=string_list,
          raise_if_empty_selection=False)
    for i in iselection:
      if (i.size() == 0):
        raise Sorry("No atom selected")
      for atom in i:
        result.append(atom)
    return result

  def is_residue_in_selection(self, i_seqs, selection):
    for i_seq in i_seqs:
      if i_seq not in selection:
        return False
    return True

  def get_i_seqs(self, atoms):
    i_seqs = []
    for atom in atoms:
      i_seqs.append(atom.i_seq)
    return i_seqs

  def extract_sequence_and_sites(self, pdb_hierarchy, selection):
    seq = []
    result = []
    counter = 0
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          if(len(rg.unique_resnames())==1):
            resname = rg.unique_resnames()[0]
            olc=amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
            atoms = rg.atoms()
            i_seqs = self.get_i_seqs(atoms)
            if(olc!="X") and self.is_residue_in_selection(i_seqs, selection):
              seq.append(olc)
              result.append(group_args(i_seq = counter, rg = rg))
              counter += 1
    return "".join(seq), result

  def _alignment(self, pdb_hierarchy,
                    pdb_hierarchy_ref,
                    params,
                    selections,
                    log=sys.stdout):
    res_match_hash = {}
    model_mseq_res_hash = {}
    model_seq, model_structures = self.extract_sequence_and_sites(
      pdb_hierarchy=pdb_hierarchy,
      selection=selections[0])
    ref_mseq_res_hash = {}
    ref_seq, ref_structures = self.extract_sequence_and_sites(
      pdb_hierarchy = pdb_hierarchy_ref,
      selection=selections[1])
    for struct in model_structures:
      model_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
    for struct in ref_structures:
      ref_mseq_res_hash[struct.i_seq] = struct.rg.atoms()[0].pdb_label_columns()[4:]
    align_obj = mmtbx.alignment.align(
      seq_a                 = model_seq,
      seq_b                 = ref_seq,
      gap_opening_penalty   = params.alignment.gap_opening_penalty,
      gap_extension_penalty = params.alignment.gap_extension_penalty,
      similarity_function   = params.alignment.similarity_matrix,
      style                 = params.alignment.alignment_style)
    alignment = align_obj.extract_alignment()
    matches = alignment.matches()
    exact_match_selections = alignment.exact_match_selections()
    exact_a = tuple(exact_match_selections[0])
    exact_b = tuple(exact_match_selections[1])
    for i, i_seq in enumerate(alignment.i_seqs_a):
      if i_seq != None:
        if alignment.i_seqs_b[i] != None and matches[i] in ['*','|']:
          res_match_hash[model_mseq_res_hash[i_seq]] = \
            ref_mseq_res_hash[alignment.i_seqs_b[i]]
    print >> log, "  --> aligning model sequence to reference sequence"
    alignment.pretty_print(block_size  = 50,
                           n_block     = 1,
                           top_name    = "model",
                           bottom_name = "ref",
                           out         = log)
    return res_match_hash

  def process_reference_groups(self,
                               pdb_hierarchy,
                               pdb_hierarchy_ref,
                               params,
                               log=sys.stdout):
    model_iseq_hash = self.build_iseq_hash(pdb_hierarchy=pdb_hierarchy)
    model_name_hash = self.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    ref_iseq_hash = self.build_iseq_hash(pdb_hierarchy=pdb_hierarchy_ref)
    sel_cache = pdb_hierarchy.atom_selection_cache()
    sel_cache_ref = pdb_hierarchy_ref.atom_selection_cache()
    match_map = {}
    #check for auto alignment compatability
    if params.auto_align == True:
      try:
        assert len(params.reference_group) == 0
      except:
        raise Sorry("""
  Cannot use reference_group selections with automatic alignment.
  Please use alignment_group selections.  See documentation for details."
  """)
      if len(params.alignment_group) == 0:
        ref_list = ['ALL']
        selection_list = ['ALL']
        sel_atoms = (self.phil_atom_selections_as_i_seqs_multiple(
                     cache=sel_cache,
                     string_list=selection_list))
        sel_atoms_ref = (self.phil_atom_selections_as_i_seqs_multiple(
                         cache=sel_cache_ref,
                         string_list=ref_list))
        selections = (sel_atoms, sel_atoms_ref)
        residue_match_map = self._alignment(pdb_hierarchy=pdb_hierarchy,
                                              pdb_hierarchy_ref=pdb_hierarchy_ref,
                                              params=params,
                                              selections=selections,
                                              log=log)
        for i_seq in sel_atoms:
          key = model_name_hash[i_seq]
          atom = key[0:4]
          res_key = key[4:]
          try:
            match_key = atom+residue_match_map[res_key]
            match_map[i_seq] = ref_iseq_hash[match_key]
          except:
            continue

      else:
        for ag in params.alignment_group:
          sel_atoms = (self.phil_atom_selections_as_i_seqs_multiple(
                       cache=sel_cache,
                       string_list=[ag.selection]))
          sel_atoms_ref = (self.phil_atom_selections_as_i_seqs_multiple(
                           cache=sel_cache_ref,
                           string_list=[ag.reference]))
          selections = (sel_atoms, sel_atoms_ref)
          residue_match_map = self._alignment(pdb_hierarchy=pdb_hierarchy,
                                                pdb_hierarchy_ref=pdb_hierarchy_ref,
                                                params=params,
                                                selections=selections,
                                                log=log)
          for i_seq in sel_atoms:
            key = model_name_hash[i_seq]
            atom = key[0:4]
            res_key = key[4:]
            try:
              match_key = atom+residue_match_map[res_key]
              match_map[i_seq] = ref_iseq_hash[match_key]
            except:
              continue
    else:
      if len(params.reference_group) == 0:
        ref_list = ['ALL']
        selection_list = ['ALL']
        sel_atoms = self.phil_atom_selections_as_i_seqs_multiple(
                        cache=sel_cache,
                        string_list=selection_list)
        for i_seq in sel_atoms:
          key = model_name_hash[i_seq]
          try:
            match_map[i_seq] = ref_iseq_hash[key]
          except:
            continue

      #specified reference groups
      for rg in params.reference_group:
        model_chain = None
        ref_chain = None
        model_res_min = None
        model_res_max = None
        ref_res_min = None
        ref_res_max = None
        sel_atoms = (self.phil_atom_selections_as_i_seqs_multiple(
                      cache=sel_cache,
                      string_list=[rg.selection]))
        sel_atoms_ref = (self.phil_atom_selections_as_i_seqs_multiple(
                      cache=sel_cache_ref,
                      string_list=[rg.reference]))
        sel_model = re.split(r"AND|OR|NOT",rg.selection.upper())
        sel_ref = re.split(r"AND|OR|NOT",rg.reference.upper())
        for sel in sel_model:
          if sel.strip().startswith("CHAIN"):
            if model_chain is None:
              model_chain = sel.strip().split(' ')[-1]
            else:
              raise Sorry("Cannot specify more than one chain per selection")
          if sel.strip().startswith("RESSEQ") or sel.strip().startswith("RESID"):
            res = sel.strip().split(' ')[-1].split(':')
            if len(res) > 1:
              if model_res_min is None and model_res_max is None:
                model_res_min = res[0]
                model_res_max = res[1]
              else:
                raise Sorry("Cannot specify more than one residue or residue range per selection")
            elif len(res) == 1:
              if model_res_min is None and model_res_max is None:
                model_res_min = res[0]
                model_res_max = res[0]
            else:
              raise Sorry("Do not understand residue selection")
        for sel in sel_ref:
          if sel.strip().startswith("CHAIN"):
            if ref_chain is None:
              ref_chain = sel.strip().split(' ')[-1]
            else:
              raise Sorry("Cannot specify more than one chain per selection")
          if sel.strip().startswith("RESSEQ") or sel.strip().startswith("RESID"):
            res = sel.strip().split(' ')[-1].split(':')
            if len(res) > 1:
              if ref_res_min is None and ref_res_max is None:
                ref_res_min = res[0]
                ref_res_max = res[1]
              else:
                raise Sorry("Cannot specify more than one residue or residue range per selection")
            elif len(res) == 1:
              if ref_res_min is None and model_res_max is None:
                ref_res_min = res[0]
                ref_res_max = res[0]
            else:
              raise Sorry("Do not understand residue selection")
        #check consistency
        assert (ref_chain is None and model_chain is None) or \
               (ref_chain is not None and model_chain is not None)
        assert (ref_res_min is None and ref_res_max is None \
                and model_res_min is None and model_res_max is None) or \
                (ref_res_min is not None and ref_res_max is not None \
                and model_res_min is not None and model_res_max is not None)
        #calculate residue offset
        offset = 0
        if (ref_res_min is not None and ref_res_max is not None \
            and model_res_min is not None and model_res_max is not None):
          offset = int(model_res_min) - int(ref_res_min)
          assert offset == (int(model_res_max) - int(ref_res_max))
        for i_seq in sel_atoms:
          key = model_name_hash[i_seq]
          if ref_chain is not None:
            if len(ref_chain)==1:
              ref_chain = ' '+ref_chain
            key = re.sub(r"(.{5}\D{3})(.{2})(.{4})",r"\1"+ref_chain+r"\3",key)
          if offset != 0:
            resnum = key[10:14]
            new_num = "%4d" % (int(resnum) - offset)
            key = re.sub(r"(.{5}\D{3})(.{2})(.{4})",r"\1"+ref_chain+new_num,key)
          try:
            assert ref_iseq_hash[key] in sel_atoms_ref
            match_map[i_seq] = ref_iseq_hash[key]
          except:
            continue
    return match_map

  def build_name_hash(self, pdb_hierarchy):
    i_seq_name_hash = dict()
    for atom in pdb_hierarchy.atoms():
      i_seq_name_hash[atom.i_seq]=atom.pdb_label_columns()
    return i_seq_name_hash

  def build_iseq_hash(self, pdb_hierarchy):
    name_i_seq_hash = dict()
    for atom in pdb_hierarchy.atoms():
      name_i_seq_hash[atom.pdb_label_columns()]=atom.i_seq
    return name_i_seq_hash

  def build_element_hash(self, pdb_hierarchy):
    i_seq_element_hash = dict()
    for atom in pdb_hierarchy.atoms():
      i_seq_element_hash[atom.i_seq]=atom.element
    return i_seq_element_hash

  def build_cbetadev_hash(self, pdb_hierarchy):
    cb = cbetadev()
    cbetadev_hash = dict()
    cbeta_out = cb.analyze_pdb(hierarchy=pdb_hierarchy)
    for line in cbeta_out[0].splitlines():
      temp = line.split(':')
      dev = temp[5]
      if dev == "dev":
        continue
      key = temp[1].upper()+temp[2].upper()+temp[3]+temp[4].rstrip()
      cbetadev_hash[key] = dev
    return cbetadev_hash

  def build_dihedral_hash(self,
                          geometry=None,
                          sites_cart=None,
                          pdb_hierarchy=None,
                          include_hydrogens=False,
                          include_main_chain=True,
                          include_side_chain=True):
    if not include_hydrogens:
      i_seq_element_hash = self.build_element_hash(pdb_hierarchy=pdb_hierarchy)
    i_seq_name_hash = self.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    dihedral_hash = dict()

    for dp in geometry.dihedral_proxies:
      try:
        #check for H atoms if required
        if not include_hydrogens:
          for i_seq in dp.i_seqs:
            if i_seq_element_hash[i_seq] == " H":
              raise StopIteration()
        #ignore backbone dihedrals
        if not include_main_chain:
          sc_atoms = False
          for i_seq in dp.i_seqs:
            if i_seq_name_hash[i_seq][0:4] not in [' CA ', ' N  ', ' C  ', ' O  ']:
              sc_atoms = True
              break
          if not sc_atoms:
            raise StopIteration()
        if not include_side_chain:
          sc_atoms = False
          for i_seq in dp.i_seqs:
            if i_seq_name_hash[i_seq][0:4] not in [' CA ', ' N  ', ' C  ', ' O  ']:
              sc_atoms = True
              break
          if sc_atoms:
            raise StopIteration()
        key = ""
        for i_seq in dp.i_seqs:
          key = key+i_seq_name_hash[i_seq]
        di = cctbx.geometry_restraints.dihedral(sites_cart=sites_cart, proxy=dp)
        dihedral_hash[key] = di.angle_model
      except StopIteration:
        pass

    #add dihedral for CB
    cbetadev_hash = self.build_cbetadev_hash(pdb_hierarchy=pdb_hierarchy)
    for cp in geometry.chirality_proxies:
      c_beta = True
      key = ""
      CAxyz = None
      Cxyz = None
      Nxyz = None
      CBxyz = None
      for i_seq in cp.i_seqs:
        key = key+i_seq_name_hash[i_seq]
        if i_seq_name_hash[i_seq][0:4] not in [' CA ', ' N  ', ' C  ', ' CB ']:
          c_beta = False
        if i_seq_name_hash[i_seq][0:4] == ' CA ':
          CAxyz = sites_cart[i_seq]
        elif i_seq_name_hash[i_seq][0:4] == ' C  ':
          Cxyz = sites_cart[i_seq]
        elif i_seq_name_hash[i_seq][0:4] == ' N  ':
          Nxyz = sites_cart[i_seq]
        elif i_seq_name_hash[i_seq][0:4] == ' CB ':
          CBxyz = sites_cart[i_seq]
          try:
            if float(cbetadev_hash[i_seq_name_hash[i_seq][4:14]]) >= 0.25:
              c_beta = False
              print "skipping C-beta restraint for %s" % i_seq_name_hash[i_seq][4:14]
          except:
              c_beta = False
      if c_beta:
        assert CAxyz is not None
        assert Cxyz is not None
        assert Nxyz is not None
        assert CBxyz is not None
        sites = [Cxyz, Nxyz, CAxyz, CBxyz]
        d = cctbx.geometry_restraints.dihedral(
          sites=sites,
          angle_ideal=0,
          weight=1)
        dihedral_hash[key] = d.angle_model
    return dihedral_hash

  def get_home_dihedral_proxies(self,
                                work_params,
                                geometry,
                                pdb_hierarchy,
                                geometry_ref,
                                sites_cart_ref,
                                pdb_hierarchy_ref,
                                log=sys.stdout):
    ss_selection = None
    residue_match_hash = {}
    reference_dihedral_proxies = cctbx.geometry_restraints.shared_dihedral_proxy()
    sigma = work_params.sigma
    limit = work_params.limit
    i_seq_name_hash = self.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    i_seq_name_hash_ref = self.build_name_hash(pdb_hierarchy=pdb_hierarchy_ref)
    reference_dihedral_hash = self.build_dihedral_hash(
                           geometry=geometry_ref,
                           sites_cart=sites_cart_ref,
                           pdb_hierarchy=pdb_hierarchy_ref,
                           include_hydrogens=work_params.hydrogens,
                           include_main_chain=work_params.main_chain,
                           include_side_chain=work_params.side_chain)
    match_map = self.process_reference_groups(
                               pdb_hierarchy=pdb_hierarchy,
                               pdb_hierarchy_ref=pdb_hierarchy_ref,
                               params=work_params,
                               log=log)
    if work_params.secondary_structure_only:
      if (not libtbx.env.has_module(name="ksdssp")):
        raise RuntimeError(
          "ksdssp module is not configured, cannot generate secondary structure reference")
      ref_ss_m = secondary_structure.manager(
                   pdb_hierarchy=pdb_hierarchy_ref,
                   xray_structure=pdb_hierarchy_ref.extract_xray_structure(),
                   sec_str_from_pdb_file=None)
      ref_ss_m.find_automatically()
      pdb_str = pdb_hierarchy_ref.as_pdb_string()
      (records, stderr) = secondary_structure.run_ksdssp_direct(pdb_str)
      sec_str_from_pdb_file = iotbx.pdb.secondary_structure.process_records(
                                records=records,
                                allow_none=True)
      if sec_str_from_pdb_file != None:
        overall_helix_selection = sec_str_from_pdb_file.overall_helix_selection()
        overall_sheet_selection = sec_str_from_pdb_file.overall_sheet_selection()
        overall_selection = overall_helix_selection +' or ' + overall_sheet_selection
        sel_cache_ref = pdb_hierarchy_ref.atom_selection_cache()
        ss_selection = (self.phil_atom_selections_as_i_seqs_multiple(
                        cache=sel_cache_ref,
                        string_list=[overall_selection]))
    for dp in geometry.dihedral_proxies:
      key = ""
      key_work = ""
      for i_seq in dp.i_seqs:
        key_work = key_work + i_seq_name_hash[i_seq]
        try:
          key = key+i_seq_name_hash_ref[match_map[i_seq]]
        except:
          continue
      try:
        reference_angle = reference_dihedral_hash[key]
        if key[5:14] == key[20:29] and \
           key[5:14] == key[35:44] and \
           key[5:14] == key[50:59] and \
           key_work[5:14] == key_work[20:29] and \
           key_work[5:14] == key_work[35:44] and \
           key_work[5:14] == key_work[50:59]:
          residue_match_hash[key_work[5:14]] = key[5:14]
      except:
        continue
      if work_params.secondary_structure_only and ss_selection != None:
          if match_map[dp.i_seqs[0]] in ss_selection and \
             match_map[dp.i_seqs[1]] in ss_selection and \
             match_map[dp.i_seqs[2]] in ss_selection and \
             match_map[dp.i_seqs[3]] in ss_selection:
            dp_add = cctbx.geometry_restraints.dihedral_proxy(
              i_seqs=dp.i_seqs,
              angle_ideal=reference_angle,
              weight=1/(1.0**2),
              limit=30.0)
            reference_dihedral_proxies.append(dp_add)
          else:
            dp_add = cctbx.geometry_restraints.dihedral_proxy(
              i_seqs=dp.i_seqs,
              angle_ideal=reference_angle,
              weight=1/(5.0**2),
              limit=15.0)
            reference_dihedral_proxies.append(dp_add)
      else:
        dp_add = cctbx.geometry_restraints.dihedral_proxy(
          i_seqs=dp.i_seqs,
          angle_ideal=reference_angle,
          weight=1/sigma**2,
          limit=limit)
        reference_dihedral_proxies.append(dp_add)

    for cp in geometry.chirality_proxies:
      key = ""
      CAsite = None
      Csite = None
      Nsite = None
      CBsite = None
      for i_seq in cp.i_seqs:
        try:
          key = key+i_seq_name_hash_ref[match_map[i_seq]]
        except:
          continue
        if i_seq_name_hash_ref[match_map[i_seq]][0:4] == ' CA ':
          CAsite = i_seq
        elif i_seq_name_hash_ref[match_map[i_seq]][0:4] == ' CB ':
          CBsite = i_seq
        elif i_seq_name_hash_ref[match_map[i_seq]][0:4] == ' C  ':
          Csite = i_seq
        elif i_seq_name_hash_ref[match_map[i_seq]][0:4] == ' N  ':
          Nsite = i_seq

      try:
        reference_angle = reference_dihedral_hash[key]
      except:
        continue
      if CAsite is None or Csite is None or CBsite is None or Nsite is None:
        continue
      i_seqs = [Csite, Nsite, CAsite, CBsite]
      if work_params.secondary_structure_only and ss_selection != None:
          if match_map[i_seqs[0]] in ss_selection and \
             match_map[i_seqs[1]] in ss_selection and \
             match_map[i_seqs[2]] in ss_selection and \
             match_map[i_seqs[3]] in ss_selection:
            dp_add = cctbx.geometry_restraints.dihedral_proxy(
              i_seqs=i_seqs,
              angle_ideal=reference_angle,
              weight=1/(1.0**2),
              limit=30.0)
            reference_dihedral_proxies.append(dp_add)
          else:
            dp_add = cctbx.geometry_restraints.dihedral_proxy(
              i_seqs=i_seqs,
              angle_ideal=reference_angle,
              weight=1/(5.0**2),
              limit=15.0)
            reference_dihedral_proxies.append(dp_add)
      else:
        dp_add = cctbx.geometry_restraints.dihedral_proxy(
          i_seqs=i_seqs,
          angle_ideal=reference_angle,
          weight=1/sigma**2,
          limit=limit)
        reference_dihedral_proxies.append(dp_add)
    self.residue_match_hash = residue_match_hash
    return reference_dihedral_proxies

  def show_reference_summary(self, log=sys.stdout):
    print >> log, "--------------------------------------------------------"
    print >> log, "Reference Model Matching Summary:"
    print >> log, "Model:              Reference:"
    keys = self.residue_match_hash.keys()
    def get_key_chain_num(res):
      return res[4:]
    keys.sort(key=get_key_chain_num)
    for key in keys:
      print >> log, "%s  <=====>  %s" % (key, self.residue_match_hash[key])
    print >> log, "--------------------------------------------------------"

  def add_reference_dihedral_proxies(self, geometry, reference_dihedral_proxies):
    geometry.reference_dihedral_proxies=reference_dihedral_proxies

  def set_rotamer_to_reference(self,
                               pdb_hierarchy,
                               pdb_hierarchy_ref,
                               xray_structure,
                               log=None,
                               quiet=False):
    if(log is None): log = sys.stdout
    print >> log, "  --> pre-correcting rotamer outliers"
    r = rotalyze()
    sa = SidechainAngles(False)
    mon_lib_srv = mmtbx.monomer_library.server.server()
    rot_list_model, coot_model = r.analyze_pdb(hierarchy=pdb_hierarchy)
    rot_list_reference, coot_reference = r.analyze_pdb(hierarchy=pdb_hierarchy_ref)
    model_hash = {}
    model_chis = {}
    reference_hash = {}
    reference_chis = {}
    for line in rot_list_model.splitlines():
      res, rotamericity, chi1, chi2, chi3, chi4, name = line.split(':')
      model_hash[res]=name

    for line in rot_list_reference.splitlines():
      res, rotamericity, chi1, chi2, chi3, chi4, name = line.split(':')
      reference_hash[res]=name

    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = r.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
                  model_chis[key] = chis
              except:
                print >> log, \
                  '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)

    for model in pdb_hierarchy_ref.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = r.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
                  reference_chis[key] = chis
              except:
                print >> log, \
                  '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)

    sites_cart_start = xray_structure.sites_cart()
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            key = '%s%5s %s' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
            try:
              if model_hash[key] == 'OUTLIER' and reference_hash[key] != 'OUTLIER':
                axis_and_atoms_to_rotate=fit_rotamers.axes_and_atoms_aa_specific(
                      residue=atom_group,
                      mon_lib_srv=mon_lib_srv,
                      remove_clusters_with_all_h=False,
                      log=None)
                m_chis = model_chis[key]
                r_chis = reference_chis[key]
                assert len(m_chis) == len(r_chis)
                assert len(m_chis) == len(axis_and_atoms_to_rotate)
                counter = 0
                residue_iselection = atom_group.atoms().extract_i_seq()
                sites_cart_residue = xray_structure.sites_cart().select(residue_iselection)
                for aa in axis_and_atoms_to_rotate:
                  axis = aa[0]
                  atoms = aa[1]
                  atom_group.atoms().set_xyz(new_xyz=sites_cart_residue)
                  new_xyz = flex.vec3_double()
                  angle_deg = r_chis[counter] - m_chis[counter]
                  if angle_deg < 0:
                    angle_deg += 360.0
                  for atom in atoms:
                    new_xyz = fit_rotamers.rotate_point_around_axis(
                                                        axis_point_1=sites_cart_residue[axis[0]],
                                                        axis_point_2=sites_cart_residue[axis[1]],
                                                        point=sites_cart_residue[atom],
                                                        angle_deg=angle_deg)
                    sites_cart_residue[atom] = new_xyz
                  sites_cart_start = sites_cart_start.set_selected(
                        residue_iselection, sites_cart_residue)
                  counter += 1
                xray_structure.set_sites_cart(sites_cart_start)
            except:
              pass
