from __future__ import division
import mmtbx.alignment
from iotbx.pdb import amino_acid_codes
from libtbx import group_args
import cctbx.geometry_restraints
from mmtbx.validation import rotalyze
from mmtbx.utils import rotatable_bonds
from mmtbx.rotamer.sidechain_angles import SidechainAngles
import mmtbx.monomer_library
from cctbx.array_family import flex
import iotbx.phil
import libtbx.load_env
from libtbx.utils import Sorry
from mmtbx import secondary_structure
from scitbx.matrix import rotate_point_around_axis
from libtbx.str_utils import make_sub_header
from mmtbx.torsion_restraints import utils
import sys, re

TOP_OUT_FLAG = True

reference_model_params = iotbx.phil.parse("""
 file = None
   .type = path
   .short_caption = Reference model
   .style = bold file_type:pdb hidden
   .multiple = True
 use_starting_model_as_reference = False
   .type = bool
   .short_caption = use starting model as reference
 sigma = 1.0
   .type = float
 limit = 15.0
   .type = float
 hydrogens = False
   .type = bool
 main_chain = True
   .type = bool
 side_chain = True
   .type = bool
 fix_outliers = True
   .type = bool
 strict_rotamer_matching = False
   .type = bool
 auto_shutoff_for_ncs = False
   .type = bool
 SSM_alignment = True
   .type = bool
 similarity = .80
   .type = float
   .short_caption = Sequence similarity cutoff
 secondary_structure_only = False
   .type = bool
 reference_group
  .multiple=True
  .optional=True
  .short_caption=Reference group
  .style = noauto auto_align menu_item parent_submenu:reference_model
{
  reference=None
    .type=atom_selection
    .short_caption=Reference selection
  selection=None
    .type=atom_selection
    .short_caption=Restrained selection
  file_name=None
    .type=path
    .optional=True
    .short_caption = Reference model for this restraint group
    .style = bold hidden
    .help = this is to used internally to disambiguate cases where multiple \
            reference models contain the same chain ID. This normally does \
            not need to be set by the user
}
""")

class reference_model(object):

  def __init__(self,
               pdb_hierarchy,
               reference_hierarchy_list=None,
               reference_file_list=None,
               mon_lib_srv=None,
               ener_lib=None,
               has_hd=False,
               params=None,
               selection=None,
               log=None):
    import time
    t0 = time.time()
    assert [reference_hierarchy_list,
            reference_file_list].count(None) == 1
    if(log is None):
      self.log = sys.stdout
    else:
      self.log = log
    self.params = params
    self.selection = selection
    self.mon_lib_srv = mon_lib_srv
    self.ener_lib = ener_lib
    pdb_hierarchy.reset_i_seq_if_necessary()
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    if self.selection is None:
      self.selection = flex.bool(len(sites_cart), True)
    self.pdb_hierarchy = pdb_hierarchy
    t1 = time.time()
    if reference_hierarchy_list is None:
      reference_hierarchy_list = \
        utils.process_reference_files(
          reference_file_list=reference_file_list,
          log=self.log)
    t2 = time.time()
    if reference_file_list is None:
      reference_file_list = []
      ref_counter = 1
      for hierarchy in reference_hierarchy_list:
        key = "ref%d" % ref_counter
        reference_file_list.append(key)
        ref_counter += 1
        hierarchy.reset_i_seq_if_necessary()
    #
    # this takes 20% of constructor time.
    self.dihedral_proxies_ref = utils.get_reference_dihedral_proxies(
        reference_hierarchy_list=reference_hierarchy_list,
        reference_file_list=reference_file_list,
        mon_lib_srv=self.mon_lib_srv,
        ener_lib=self.ener_lib,
        log=self.log)
    self.i_seq_name_hash = utils.build_name_hash(
                             pdb_hierarchy=self.pdb_hierarchy)
    t3 = time.time()

    #reference model components
    self.sites_cart_ref = {}
    self.pdb_hierarchy_ref = {}
    self.i_seq_name_hash_ref = {}
    self.reference_dihedral_hash = {}
    self.reference_file_list = reference_file_list

    #triage reference model files
    for file, hierarchy in zip(reference_file_list,
                               reference_hierarchy_list):
      self.sites_cart_ref[file] = hierarchy.atoms().extract_xyz()
      self.pdb_hierarchy_ref[file] = hierarchy
      self.i_seq_name_hash_ref[file] = \
        utils.build_name_hash(
          pdb_hierarchy=hierarchy)
      self.reference_dihedral_hash[file] = \
        self.build_dihedral_hash(
          dihedral_proxies=self.dihedral_proxies_ref[file],
          sites_cart=self.sites_cart_ref[file],
          pdb_hierarchy=hierarchy,
          include_hydrogens=self.params.hydrogens,
          include_main_chain=self.params.main_chain,
          include_side_chain=self.params.side_chain)
    t4 = time.time()
    self.match_map = None
    self.proxy_map = None
    self.build_reference_dihedral_proxy_hash()
    t5 = time.time()
    #
    # This takes 80% of constructor time!!!
    self.get_reference_dihedral_proxies()
    t6 = time.time()
    # print "Timing in ref dih proxies constructor: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f. Total: %.3f" % (
    #     t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t6-t0)

  def top_out_function(self, x, weight, top):
    return top*(1-exp(-weight*x**2/top))

  def top_out_gradient(self, x, weight, top):
    return (2*weight*x)*exp(-(weight*x**2)/top)

  def top_out_curvature(self, x, weight, top):
    return (2*weight*(top - 2*weight*x**2))/top**2*exp(-(weight*x**2)/top)

  def extract_sequence_and_sites(self, pdb_hierarchy, selection):
    assert 0, "Seems that nobody is using this"
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
            i_seqs = atoms.extract_i_seq()
            if(olc!="X") and utils.is_residue_in_selection(i_seqs, selection):
              seq.append(olc)
              result.append(group_args(i_seq = counter, rg = rg))
              counter += 1
    return "".join(seq), result

  def process_reference_groups(self,
                               pdb_hierarchy,
                               pdb_hierarchy_ref,
                               params,
                               log=None):
    if(log is None): log = sys.stdout
    model_iseq_hash = utils.build_i_seq_hash(pdb_hierarchy=pdb_hierarchy)
    model_name_hash = {v:k for k, v in model_iseq_hash.iteritems()}
    sel_cache = pdb_hierarchy.atom_selection_cache()
    ref_iseq_hash = {}
    sel_cache_ref = {}
    for key in pdb_hierarchy_ref:
      ref_hier = pdb_hierarchy_ref[key]
      ref_iseq_hash[key] = \
        utils.build_i_seq_hash(pdb_hierarchy=ref_hier)
      sel_cache_ref[key] = ref_hier.atom_selection_cache()
    match_map = {}
    for file in self.reference_file_list:
      match_map[file] = {}
    if len(params.reference_group) == 0:
      matched_model_chains = []
      reference_matches = get_matching_chains(
        pdb_hierarchy=pdb_hierarchy,
        pdb_hierarchy_ref=pdb_hierarchy_ref,
        reference_file_list=self.reference_file_list,
        params=params,
        log=log)
      new_reference_groups = ""
      for file in reference_matches.keys():
        matches = reference_matches[file]
        for match in matches:
          model_chain = match[0]
          ref_chain = match[1]
          if model_chain in matched_model_chains:
            raise Sorry("multiple reference models for one chain not "+\
                        "currently supported")
          matched_model_chains.append(model_chain)
          new_reference_groups += "  reference_group {\n"
          new_reference_groups += "   reference = "+ref_chain+"\n"
          new_reference_groups += "   selection = "+model_chain+"\n"
          new_reference_groups += "   file_name = "+file+"\n"
          new_reference_groups += "  }\n"
        new_reference_groups_as_phil = iotbx.phil.parse(
          new_reference_groups)
        params.reference_group = \
          reference_model_params.\
            fetch(new_reference_groups_as_phil)\
            .extract().\
            reference_group
    else:
      for rg in params.reference_group:
        if rg.file_name is None:
          if len(self.reference_file_list) > 1:
            raise Sorry("Ambiguous reference group selection - please "+\
                        "specify file name for reference group")
          #only one reference file specified
          rg.file_name = self.reference_file_list[0]
    for rg in params.reference_group:
      model_chain = None
      ref_chain = None
      model_res_min = None
      model_res_max = None
      ref_res_min = None
      ref_res_max = None
      #check for selection sanity
      sel_model = re.split(r"AND|OR|NOT",rg.selection.upper())
      sel_ref = re.split(r"AND|OR|NOT",rg.reference.upper())
      for sel in sel_model:
        if sel.strip().strip('(').strip(')').startswith("CHAIN"):
          if model_chain is None:
            model_chain = sel.strip().strip('(').strip(')').split(' ')[-1].\
              strip("'").strip('"')
          else:
            raise Sorry("Cannot specify more than one chain per selection")
        if sel.strip().strip('(').strip(')').startswith("RESSEQ") or \
           sel.strip().strip('(').strip(')').startswith("RESID"):
          res = sel.strip().strip('(').strip(')').split(' ')[-1].split(':')
          if len(res) > 1:
            if model_res_min is None and model_res_max is None:
              model_res_min = res[0]
              model_res_max = res[1]
            else:
              raise Sorry("Cannot specify more than one residue "+\
                "or residue range per selection")
          elif len(res) == 1:
            if model_res_min is None and model_res_max is None:
              model_res_min = res[0]
              model_res_max = res[0]
          else:
            raise Sorry("Do not understand residue selection")
      for sel in sel_ref:
        if sel.strip().strip('(').strip(')').startswith("CHAIN"):
          if ref_chain is None:
            ref_chain = sel.strip().strip('(').strip(')').split(' ')[-1].\
              strip("'").strip('"')
            if ref_chain.startswith("'") or ref_chain.startswith('"'):
              ref_chain = ref_chain.lstrip()
            if ref_chain.endswith("'") or ref_chain.endswith('"'):
              ref_chain = ref_chain.rstrip()
          else:
            raise Sorry("Cannot specify more than one chain per selection")
        if sel.strip().strip('(').strip(')').startswith("RESSEQ") or \
           sel.strip().strip('(').strip(')').startswith("RESID"):
          res = sel.strip().strip('(').strip(')').split(' ')[-1].split(':')
          if len(res) > 1:
            if ref_res_min is None and ref_res_max is None:
              ref_res_min = res[0]
              ref_res_max = res[1]
            else:
              raise Sorry("Cannot specify more than one residue "+\
                "or residue range per selection")
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
      #prep for SSM alignment
      file = rg.file_name
      sel = sel_cache.selection(string=rg.selection)
      sel_atoms = sel.iselection()
      sel_ref = sel_cache_ref[file].selection(string=rg.reference)
      mod_h = utils.hierarchy_from_selection(
                pdb_hierarchy=pdb_hierarchy,
                selection = sel,
                log=log).models()[0].chains()[0]
      ref_h = utils.hierarchy_from_selection(
                pdb_hierarchy=pdb_hierarchy_ref[file],
                selection = sel_ref,
                log=log).models()[0].chains()[0]
      ssm = None
      if params.SSM_alignment:
        try: #do SSM alignment
          ssm, ssm_align = utils._ssm_align(
                      reference_chain = ref_h,
                      moving_chain = mod_h)
        except RuntimeError, e:
          if (str(e) != "can't make graph for first structure" and \
              str(e) != "can't make graph for second structure" and \
              str(e) != "secondary structure does not match"):
            raise e
          else:
            print >> log, "SSM alignment failed..."
      if ssm != None:
        for pair in ssm_align.pairs:
          model_res = pair[0]
          ref_res = pair[1]
          if model_res is None or ref_res is None:
            continue
          temp_model_atoms = {}
          temp_ref_atoms = {}
          for atom in model_res.atoms():
            atom_temp = atom.pdb_label_columns()+atom.segid
            temp_model_atoms[atom.name] = model_iseq_hash[atom_temp]
          for atom in ref_res.atoms():
            atom_temp = atom.pdb_label_columns()+atom.segid
            temp_ref_atoms[atom.name] = ref_iseq_hash[file][atom_temp]
          for key in temp_model_atoms.keys():
            ref_atom = temp_ref_atoms.get(key)
            if ref_atom != None:
              match_map[file][temp_model_atoms[key]] = temp_ref_atoms[key]

      else: #ssm not selected, or failed
        print >> log, "trying simple matching..."
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
            elif len(ref_chain)==0:
              ref_chain = '  '
            key = \
              re.sub(r"(.{5}\D{3})(.{2})(.{4})",r"\1"+ref_chain+r"\3",key)
          if offset != 0:
            resnum = key[10:14]
            new_num = "%4d" % (int(resnum) - offset)
            key = \
              re.sub(r"(.{5}\D{3})(.{2})(.{4})",r"\1"+ref_chain+new_num,key)
          cur_ref = ref_iseq_hash[file].get(key)
          if cur_ref is None:
            #try no alternate
            key = \
              re.sub(r"(.{4}).(\D{3})(.{2})(.{4})",
                     r"\1"+" "+r"\2"+ref_chain+r"\4",key)
          cur_ref = ref_iseq_hash[file].get(key)
          if cur_ref is None:
            continue
          else:
            match_map[file][i_seq] = ref_iseq_hash[file][key]
    return match_map

  def build_reference_dihedral_proxy_hash(self):
    self.reference_dihedral_proxy_hash = {}
    for ref in self.dihedral_proxies_ref.keys():
      self.reference_dihedral_proxy_hash[ref] = {}
      proxies = self.dihedral_proxies_ref[ref]
      for dp in proxies:
        key = ""
        for i_seq in dp.i_seqs:
          key += self.i_seq_name_hash_ref[ref][i_seq]
        self.reference_dihedral_proxy_hash[ref][key] = dp

  def build_dihedral_hash(self,
                          dihedral_proxies=None,
                          sites_cart=None,
                          pdb_hierarchy=None,
                          include_hydrogens=False,
                          include_main_chain=True,
                          include_side_chain=True):
    if not include_hydrogens:
      i_seq_element_hash = \
        utils.build_element_hash(pdb_hierarchy=pdb_hierarchy)
    i_seq_name_hash = \
      utils.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    dihedral_hash = dict()

    for dp in dihedral_proxies:
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
            if i_seq_name_hash[i_seq][0:4] \
              not in [' CA ', ' N  ', ' C  ', ' O  ']:
              sc_atoms = True
              break
          if not sc_atoms:
            raise StopIteration()
        if not include_side_chain:
          sc_atoms = False
          for i_seq in dp.i_seqs:
            if i_seq_name_hash[i_seq][0:4] \
              not in [' CA ', ' N  ', ' C  ', ' O  ']:
              sc_atoms = True
              break
          if sc_atoms:
            raise StopIteration()
        key = ""
        for i_seq in dp.i_seqs:
          key = key+i_seq_name_hash[i_seq]
        di = \
          cctbx.geometry_restraints.dihedral(
            sites_cart=sites_cart,
            proxy=dp)
        dihedral_hash[key] = di.angle_model
      except StopIteration:
        pass
    return dihedral_hash

  def get_reference_dihedral_proxies(self):
    residue_match_hash = {}
    complete_dihedral_proxies = utils.get_complete_dihedral_proxies(
                                  pdb_hierarchy=self.pdb_hierarchy,
                                  mon_lib_srv=self.mon_lib_srv,
                                  ener_lib=self.ener_lib,
                                  log=self.log)
    self.reference_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
    sigma = self.params.sigma
    limit = self.params.limit
    match_map = self.process_reference_groups(
                               pdb_hierarchy=self.pdb_hierarchy,
                               pdb_hierarchy_ref=self.pdb_hierarchy_ref,
                               params=self.params,
                               log=self.log)
    self.match_map = match_map
    ref_ss_m = None
    ss_selection = None
    if self.params.secondary_structure_only:
      if (not libtbx.env.has_module(name="ksdssp")):
        raise RuntimeError(
          "ksdssp module is not configured, "+\
          "cannot generate secondary structure reference")
      ref_ss_m = {}
      ss_selection = {}
      for file in self.reference_file_list:
        ref_ss_m[file] = secondary_structure.manager(
          pdb_hierarchy=self.pdb_hierarchy_ref[file],
          sec_str_from_pdb_file=None)
        pdb_str = self.pdb_hierarchy_ref[file].as_pdb_string()
        (records, stderr) = secondary_structure.run_ksdssp_direct(pdb_str)
        sec_str_from_pdb_file = iotbx.pdb.secondary_structure.process_records(
                                  records=records,
                                  allow_none=True)
        if sec_str_from_pdb_file != None:
          overall_helix_selection = \
            sec_str_from_pdb_file.overall_helix_selection()
          overall_sheet_selection = \
            sec_str_from_pdb_file.overall_sheet_selection()
          overall_selection = \
            overall_helix_selection +' or ' + overall_sheet_selection
          sel_cache_ref = self.pdb_hierarchy_ref[file].atom_selection_cache()
          bsel = sel_cache_ref.selection(string=overall_selection)
          if bsel.all_eq(False):
            raise Sorry("No atom selected")
          ss_selection[file] = bsel
    for dp in complete_dihedral_proxies:
      key_work = ""
      complete = True
      for i_seq in dp.i_seqs:
        if not self.selection[i_seq]:
          complete = False
      if not complete:
        continue
      for i_seq in dp.i_seqs:
        key_work = key_work + self.i_seq_name_hash[i_seq]
      #find matching key
      key = None
      file_match = None
      for file in self.reference_file_list:
        if key is not None:
          continue
        else:
          key = ""
          ref_match = True
          for i_seq in dp.i_seqs:
            if ref_match:
              map_part = match_map[file].get(i_seq)
              if map_part is not None:
                key_part = \
                  self.i_seq_name_hash_ref[file].get(map_part)
                if key_part is None:
                  ref_match = False
                  key = None
                  file_match = None
                else:
                  key = key+key_part
              else:
                ref_match = False
                key = None
                file_match = None
            if key is not None:
              file_match = file
      try:
        reference_angle = self.reference_dihedral_hash[file_match][key]
        if key[5:18] == key[24:37] and \
           key[5:18] == key[43:56] and \
           key_work[5:18] == key_work[24:37] and \
           key_work[5:18] == key_work[43:56]:
          residue_match_hash[key_work[5:14]] = (file_match, key[5:14])
      except Exception:
        continue
      w_limit = limit
      w_weight = 1/sigma**2
      if (self.params.secondary_structure_only and ss_selection is not None
          and ss_selection[file_match] is not None):
        limit2 = 15.0
        w_weight = 0.04
        if (ss_selection[file_match][match_map[file_match][dp.i_seqs[0]]] and
           ss_selection[file_match][match_map[file_match][dp.i_seqs[1]]] and
           ss_selection[file_match][match_map[file_match][dp.i_seqs[2]]] and
           ss_selection[file_match][match_map[file_match][dp.i_seqs[3]]]):
          limit2 = 30.0
          w_weight = 1
        dp_add = cctbx.geometry_restraints.dihedral_proxy(
            i_seqs=dp.i_seqs,
            angle_ideal=reference_angle,
            weight=w_weight,
            limit=w_limit,
            top_out=TOP_OUT_FLAG)
        self.reference_dihedral_proxies.append(dp_add)
      else:
        dp_add = cctbx.geometry_restraints.dihedral_proxy(
            i_seqs=dp.i_seqs,
            angle_ideal=reference_angle,
            weight=1/sigma**2,
            limit=limit,
            top_out=TOP_OUT_FLAG)
        self.reference_dihedral_proxies.append(dp_add)
    self.residue_match_hash = residue_match_hash

  def show_reference_summary(self, log=None):
    if(log is None): log = sys.stdout
    print >> log, "--------------------------------------------------------"
    print >> log, "Reference Model Matching Summary:"
    keys = self.residue_match_hash.keys()
    def get_key_chain_num(res):
      return res[4:]
    keys.sort(key=get_key_chain_num)
    for file in self.reference_file_list:
      print >> log, "\nreference file: %s\n" % file
      print >> log, "Model:              Reference:"
      for key in keys:
        if self.residue_match_hash[key][0] == file:
          print >> log, "%s  <=====>  %s" % \
            (key, self.residue_match_hash[key][1])
    print >> log, "\nTotal # of matched residue pairs: %d" % len(keys)
    print >> log, "Total # of reference model restraints: %d" % \
      len(self.reference_dihedral_proxies)
    print >> log, "--------------------------------------------------------"

  def add_reference_dihedral_proxies(self, geometry):
    geometry.reference_dihedral_proxies= \
      self.reference_dihedral_proxies

  def set_rotamer_to_reference(self,
                               xray_structure,
                               mon_lib_srv=None,
                               log=None,
                               quiet=False):
    if self.mon_lib_srv is None:
      self.mon_lib_srv = mon_lib_srv
    assert isinstance(self.mon_lib_srv, mmtbx.monomer_library.server.server)
    if(log is None): log = self.log
    make_sub_header(
      "Correcting rotamer outliers to match reference model",
      out=log)
    sa = SidechainAngles(False)
    r = rotalyze.rotalyze(pdb_hierarchy=self.pdb_hierarchy)
    rot_list_reference = {}
    coot_reference = {}
    for key in self.pdb_hierarchy_ref.keys():
      hierarchy = self.pdb_hierarchy_ref[key]
      rot_list_reference[key] = \
        rotalyze.rotalyze(pdb_hierarchy=hierarchy)
    model_hash = {}
    model_chis = {}
    reference_hash = {}
    reference_chis = {}
    model_outliers = 0
    for rot in r.results:
      model_hash[rot.id_str()] = rot.rotamer_name
      if rot.rotamer_name == "OUTLIER":
        model_outliers += 1

    for key in rot_list_reference.keys():
      reference_hash[key] = {}
      for rot in rot_list_reference[key].results:
        reference_hash[key][rot.id_str()] = rot.rotamer_name

    print >> log, "** evaluating rotamers for working model **"
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
            all_dict = rotalyze.construct_complete_sidechain(residue_group)
            for atom_group in residue_group.atom_groups():
              try:
                atom_dict = all_dict.get(atom_group.altloc)
                chis = sa.measureChiAngles(atom_group, atom_dict)
                if chis is not None:
                  key = utils.id_str(
                          chain_id=chain.id,
                          resseq=residue_group.resseq,
                          resname=atom_group.resname,
                          icode=residue_group.icode,
                          altloc=atom_group.altloc)
                  model_chis[key] = chis
              except Exception:
                print >> log, \
                  '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                      chain.id, residue_group.resid(),
                      atom_group.altloc+atom_group.resname)
    if model_outliers == 0:
      print >> log, "No rotamer outliers detected in working model"
      return
    else:
      print >> log, "Number of rotamer outliers: %d" % model_outliers

    print >> log, "\n** evaluating rotamers for reference model **"
    for file in self.pdb_hierarchy_ref.keys():
      hierarchy = self.pdb_hierarchy_ref[file]
      reference_chis[file] = {}
      for model in hierarchy.models():
        for chain in model.chains():
          for residue_group in chain.residue_groups():
              all_dict = rotalyze.construct_complete_sidechain(residue_group)
              for atom_group in residue_group.atom_groups():
                try:
                  atom_dict = all_dict.get(atom_group.altloc)
                  chis = sa.measureChiAngles(atom_group, atom_dict)
                  if chis is not None:
                    key = utils.id_str(
                            chain_id=chain.id,
                            resseq=residue_group.resseq,
                            resname=atom_group.resname,
                            icode=residue_group.icode,
                            altloc=atom_group.altloc)
                    reference_chis[file][key] = chis
                except Exception:
                  print >> log, \
                    '  %s%5s %s is missing some sidechain atoms, **skipping**' % (
                        chain.id, residue_group.resid(),
                        atom_group.altloc+atom_group.resname)

    print >> log, "\n** fixing outliers **"
    sites_cart_start = xray_structure.sites_cart()
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          if len(residue_group.conformers()) > 1:
            print >> log, "  %s%5s %s has multiple conformations, **skipping**" % (
              chain.id, residue_group.resid(),
              " "+residue_group.atom_groups()[0].resname)
            continue
          for conformer in residue_group.conformers():
            for residue in conformer.residues():
              if residue.resname == "PRO":
                continue
              key = utils.id_str(
                      chain_id=chain.id,
                      resseq=residue_group.resseq,
                      resname=residue_group.atom_groups()[0].resname,
                      icode=residue_group.icode,
                      altloc=conformer.altloc)
              if len(chain.id) == 1:
                chain_id = " "+chain.id
              else:
                chain_id = chain.id
              file_key = '%s%s%s' %(residue.resname,
                                    chain_id,
                                    residue_group.resid())
              file_key = file_key.strip()
              file_match = self.residue_match_hash.get(file_key)
              if file_match is not None:
                file = file_match[0]
              else:
                continue
              model_rot = model_hash.get(key)
              reference_rot = reference_hash[file].get(self.one_key_to_another(file_match[1]))
              m_chis = model_chis.get(key)
              r_chis = reference_chis[file].get(self.one_key_to_another(file_match[1]))
              assert len(m_chis) == len(r_chis)
              if model_rot is not None and reference_rot is not None and \
                 m_chis is not None and r_chis is not None:
                if (model_rot == 'OUTLIER' and \
                    reference_rot != 'OUTLIER'): # or \
                    #atom_group.resname in ["LEU", "VAL", "THR"]:
                  self.change_residue_rotamer_in_place(
                      sites_cart_start,residue, m_chis,r_chis,self.mon_lib_srv)
                  xray_structure.set_sites_cart(sites_cart_start)

                elif self.params.strict_rotamer_matching and \
                  (model_rot != 'OUTLIER' and reference_rot != 'OUTLIER'):
                  if model_rot != reference_rot:
                    self.change_residue_rotamer_in_place(
                        sites_cart_start,residue, m_chis,r_chis,self.mon_lib_srv)
                    xray_structure.set_sites_cart(sites_cart_start)

  def one_key_to_another(self,key):
    # Work-around function, don't have time to dig into 10 different cache
    # types... Probably better data structure could be suggested to handle
    # all needed information.

    # sp = key.split()
    # assert len(sp) == 3
    var1 = " %s  %s" % (key[4:], key[:3])
    return var1


  def change_residue_rotamer_in_place(self,sites_cart, residue,
      m_chis, r_chis, mon_lib_srv):
    axis_and_atoms_to_rotate= \
      rotatable_bonds.axes_and_atoms_aa_specific(
          residue=residue,
          mon_lib_srv=mon_lib_srv,
          remove_clusters_with_all_h=True,
          log=None)
    if axis_and_atoms_to_rotate is None:
      return
    assert len(m_chis) == len(axis_and_atoms_to_rotate)
    counter = 0
    residue_iselection = residue.atoms().extract_i_seq()
    sites_cart_residue = sites_cart.select(residue_iselection)
    for aa in axis_and_atoms_to_rotate:
      axis = aa[0]
      atoms = aa[1]
      residue.atoms().set_xyz(new_xyz=sites_cart_residue)
      new_xyz = flex.vec3_double()
      angle_deg = r_chis[counter] - m_chis[counter]
      if angle_deg < 0:
        angle_deg += 360.0
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
                    axis_point_1=sites_cart_residue[axis[0]],
                    axis_point_2=sites_cart_residue[axis[1]],
                    point=sites_cart_residue[atom],
                    angle=angle_deg, deg=True)
        sites_cart_residue[atom] = new_xyz
      sites_cart = sites_cart.set_selected(
            residue_iselection, sites_cart_residue)
      counter += 1

  def remove_restraints_with_ncs_matches(self,
                                         ncs_dihedral_proxies,
                                         ncs_match_hash):
    proxy_list = []
    remaining_proxies = cctbx.geometry_restraints.shared_dihedral_proxy()
    remaining_match_hash = {}
    for dp in ncs_dihedral_proxies:
      proxy_list.append(dp.i_seqs)
    print len(self.reference_dihedral_proxies)
    for dp in self.reference_dihedral_proxies:
      if dp.i_seqs not in proxy_list:
        remaining_proxies.append(dp)
    for key in self.residue_match_hash:
      found_match = False
      for key2 in ncs_match_hash:
        if key == key2:
          found_match = True
        else:
          for match in ncs_match_hash[key2]:
            if key == match:
              found_match = True
      if not found_match:
        remaining_match_hash[key] = self.residue_match_hash[key]
    print len(remaining_proxies)
    self.reference_dihedral_proxies = remaining_proxies
    self.residue_match_hash = remaining_match_hash
    print >> self.log, "\n**Removed reference restraints that overlap "+ \
                       "with torsion NCS restraints**\n"
    print >> self.log, "Updated Reference Model Restraints:"
    self.show_reference_summary()

def get_matching_chains(pdb_hierarchy,
                        pdb_hierarchy_ref,
                        reference_file_list,
                        params,
                        log):
  reference_matches = {}
  print >> log, "determining reference matches automatically..."
  chains = pdb_hierarchy.models()[0].chains()
  atom_labels = list(pdb_hierarchy.atoms_with_labels())
  segids = flex.std_string([ a.segid for a in atom_labels ])
  use_segid = not segids.all_eq('    ')
  am = utils.alignment_manager(
         pdb_hierarchy=pdb_hierarchy,
         use_segid=use_segid)
  for file in reference_file_list:
    pair_hash = {}
    reference_matches[file] = []
    hierarchy_ref = pdb_hierarchy_ref[file]
    chains_ref = hierarchy_ref.models()[0].chains()
    atom_labels_ref = list(hierarchy_ref.atoms_with_labels())
    segids_ref = flex.std_string([ a.segid for a in atom_labels_ref ])
    use_segid_ref = not segids_ref.all_eq('    ')
    am_ref = utils.alignment_manager(
               pdb_hierarchy=hierarchy_ref,
               use_segid=use_segid_ref)
    for i, chain_i in enumerate(chains):
      found_conformer = False
      for conformer in chain_i.conformers():
        if not conformer.is_protein() and not conformer.is_na():
          continue
        else:
          found_conformer = True
      if not found_conformer:
        continue
      segid_i = utils.get_unique_segid(chain_i)
      if segid_i == None:
        print >> log, \
          "chain %s has conflicting segid values - skipping" % chain_i.id
        continue
      if (use_segid) :
        chain_i_str = "chain '%s' and segid '%s'" % \
          (chain_i.id, segid_i)
      else :
        chain_i_str = "chain '%s'" % chain_i.id
      for j, chain_j in enumerate(chains_ref):
        found_conformer = False
        for conformer in chain_j.conformers():
          if not conformer.is_protein() and not conformer.is_na():
            continue
          else:
            found_conformer = True
        if not found_conformer:
          continue
        segid_j = utils.get_unique_segid(chain_j)
        if segid_j == None:
          continue
        if (use_segid_ref) :
          chain_j_str = "chain '%s' and segid '%s'" % (chain_j.id, segid_j)
        else :
          chain_j_str = "chain '%s'" % chain_j.id
        seq_pair = (am.sequences[chain_i_str],
                    am_ref.sequences[chain_j_str])
        seq_pair_padded = (am.padded_sequences[chain_i_str],
                           am_ref.padded_sequences[chain_j_str])
        struct_pair = (am.structures[chain_i_str],
                       am_ref.structures[chain_j_str])
        residue_match_map = \
          utils._alignment(
            sequences=seq_pair,
            padded_sequences=seq_pair_padded,
            structures=struct_pair,
            log=log)
        #require length of matches to be similar to shorter chain
        if ( (len(residue_match_map) / min(chain_i.residue_groups_size(),
                                           chain_j.residue_groups_size()))
              >= params.similarity ):
          key = (chain_i_str, chain_j_str)
          pair_key = (chain_i.id, segid_i)
          match_key = (chain_j.id, segid_j)
          if (not pair_key in pair_hash) :
            pair_hash[pair_key] = []
          pair_hash[pair_key].append(match_key)
    for key in pair_hash.keys():
      if (use_segid) :
        chain_str = "chain '%s' and segid '%s'" % (key[0], key[1])
      else :
        chain_str = "chain '%s'" % (key[0])
      exact_match = None
      first_match = None
      for match in pair_hash[key]:
        if (use_segid) :
          ref_str = "chain '%s' and segid '%s'" % \
            (match[0], match[1])
        else :
          ref_str = "chain '%s'" % (match[0])
        if first_match == None:
          first_match = ref_str
        if key[0] == match[0]:
          exact_match = ref_str
      if exact_match != None:
        match_pair = (chain_str, exact_match)
      else:
        match_pair = (chain_str, first_match)
      reference_matches[file].append(match_pair)
    reference_matches[file].sort()
  return reference_matches
