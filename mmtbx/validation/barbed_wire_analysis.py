from __future__ import absolute_import, division, print_function
import os, sys
import iotbx
from libtbx.utils import null_out

packing_quality_key = ["no packing", "underpacked", "marginal packing", "well packed"]  # indices 0,1,2,3

def loadModel(filename):  # from suitename/suites.py
  from iotbx.data_manager import DataManager
  dm = DataManager()  # Initialize the DataManager and call it dm
  dm.set_overwrite(True)  # tell the DataManager to overwrite files with the same name
  model = dm.get_model(filename)
  return model

class single_contact():
  def __init__(self, probe_line):
    x = probe_line.split(':')
    # Probe Unformatted Output:
    # name:pat:type:srcAtom:targAtom:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    # for condensed output we have:
    # name:pat:type:srcAtom:targAtom:*dotcount*:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    ###'name' is set by the user on the command line
    ###'pat' is one of 1->1, 1->2, or 2->1; where 1 is src and 2 is targ.
    ###'type' is one of wc, cc, so, bo, hb (wide/close contact, small/bad overlap, h-bond).
    ###'srcAtom' and 'targAtom' follow the pattern CNNNNITTT AAAAL, where C is chain, N is number, I is insertion code, T is residue type, A is atom name, and L is alternate conformation flag.
    ###'*dotcount*' is condensed-output-only, and gives the number of dots in the contact
    ###'min-gap' is the distance between atoms, minus their van der Waals radii; i.e., the distance of closest approach for their vdW surfaces. gap is the distance between vdW surfaces at the current dot. Negative values indicate overlap (clashes or H-bonds).
    ###'x','y','z' is a point on the vdW surface; 'spX','spY','spZ' is tip of spike, if any (same as x,y,z for contacts)
    ###'score' is "this dot's contribution to the [Probe] score" (scaled already? YES)
    ###'stype' and 'ttype' are heavy-atom element name (C, N, O, etc)

    self.name = x[0]
    self.pattern = x[1]
    self.interactiontype = x[2]

    self.srcAtom = x[3]
    self.srcChain = self.srcAtom[0:2].strip()
    self.srcNum = int(self.srcAtom[2:6].strip())
    self.srcNumStr = self.srcAtom[2:6]
    self.srcIns = self.srcAtom[6:7]  # .strip()
    self.srcResname = self.srcAtom[7:10].strip()
    ###if srcResname == 'HOH': continue  # skip waters
    self.srcAtomname = self.srcAtom[11:15]  # .strip()
    self.srcAlt = self.srcAtom[15:16].strip()

    self.trgAtom = x[4]
    # going to count dots per bond as a measure of strength instead
    self.trgChain = self.trgAtom[0:2].strip()
    self.trgNum = int(self.trgAtom[2:6].strip())
    self.trgNumStr = self.trgAtom[2:6]
    self.trgIns = self.trgAtom[6:7]  # .strip()
    self.trgResname = self.trgAtom[7:10].strip()
    ### if trgResname == 'HOH': continue #skip waters
    self.trgAtomname = self.trgAtom[11:15]  # .strip()
    self.trgAlt = self.trgAtom[15:16].strip()

    self.dotcount = x[5]
    self.mingap = x[6]

    ##mc_contact = self.contact_type(srcAtomname,trgAtomname)
    self.srcAtomid = ",".join(
      [self.srcChain.strip(), self.srcAlt.strip(), self.srcNumStr, self.srcIns, self.srcAtomname])
    self.srcResid = ','.join([self.srcChain.strip(), self.srcNumStr])
    self.trgAtomid = ",".join(
      [self.trgChain.strip(), self.trgAlt.strip(), self.trgNumStr, self.trgIns, self.trgAtomname])
    self.trgResid = ','.join([self.trgChain.strip(), self.trgNumStr])
    self.contactId = ",".join([self.srcAtomid, self.trgAtomid])


def do_probe(model):
  from libtbx import easy_run

  pdblines = model.get_hierarchy().as_pdb_string()
  probe_temp = open("probe_contacts_temp.pdb", "w")
  print(pdblines, file=probe_temp)
  probe_temp.close()
  cmd = "phenix.reduce -noflip probe_contacts_temp.pdb"
  reduce_run = easy_run.fully_buffered(cmd).stdout_lines
  probe_temp = open("probe_contacts_temp.pdb", "w")
  for line in reduce_run:
    print(line, file=probe_temp)
  probe_temp.close()
  cmd = "phenix.probe -u -con -self -mc ALL probe_contacts_temp.pdb"
  probe_run = easy_run.fully_buffered(cmd).stdout_lines

  contacts = []
  for line in probe_run:
    if not line.strip(): continue  # averts an IndexError problem with empty lines
    contacts.append(single_contact(line))

  sort_order = {'bo': 1, 'hb': 2, 'so': 3, 'cc': 4, 'wc': 5}
  contacts.sort(key=lambda x: sort_order[x.interactiontype], reverse=True)
  contact_dict = {}
  for contact in contacts:
    contact_dict[contact.contactId] = contact
  # a single atom pair can appear in multiple categories
  # this sort->dict means that each atom pair gets a single contact entry of the highest priority type
  # lower priority types get overwritten by higher (eg 'so' is stored first, then overwritten by 'bo' if present)
  os.remove("probe_contacts_temp.pdb")  # cleanup
  return contact_dict

class alphafold_chunk():
  def __init__(self, res_list, prediction_type):
    self.members = [r for r in res_list]
    self.prediction_type = prediction_type

  @property
  def start(self):
    return self.members[0]

  @property
  def end(self):
    return self.members[-1]

  def as_selection_string(self):
    # resid is "chain,resseq"
    chain = self.start.split(',')[0]
    start_res = self.start.split(',')[1].strip()
    end_res = self.end.split(',')[1].strip()
    return "(chain %s and resseq %s:%s)" % (chain, start_res, end_res)

  def add_to_start(self, resid):
    self.members.insert(0, resid)

  def add_to_end(self, resid):
    self.members.append(resid)

  def remove_from_start(self):
    return self.members.pop(0)

  def remove_from_end(self):
    return self.members.pop(-1)

class predicted_residue():
  def __init__(self, rg):
    self.chain = rg.parent().id
    self.resseq = rg.resseq
    self.resnum = int(self.resseq.strip())
    self.resid = ','.join([self.chain, self.resseq])  # assume no alternates, no inserts in predictions
    self.caxyz, self.plddt = self.find_ca_plddt(rg)
    self.heavy_atom_count = self.find_heavy_atom_count(rg)
    self.out_rama = None
    self.rama_high_psi = None
    self.rama_high_phi = None
    self.out_omega = None
    self.out_geom = None
    self.out_cablam = None
    self.ss_type = None  # "helix","sheet", or None for coil
    self.ss_index = None  # index of helix or sheet, used to tell if residues are part of same element
    self.packing_hb = []
    self.packing_vdw = []
    self.packing_bo = []
    self.severe_clashes = 0  # cutoff set in add_contacts, this may indicate knots/interpenetration
    self.packing_quality = None
    self.feedback = ''
    self.text_code = '      '

  def find_ca_plddt(self, rg):
    for atom in rg.atoms():
      if atom.name == " CA ":
        return atom.xyz, atom.b
    return None, None

  def find_heavy_atom_count(self, rg):
    # print(rg.atom_groups()[0].resname, len(rg.atoms()))
    return (len(rg.atoms()))

  def packing_contact_count(self):
    # sum of contacts with this residue
    # hb and bo types my be upweighted
    return len(self.packing_hb) + len(self.packing_vdw) + len(self.packing_bo)

  def has_any_validation_errors(self):
    if self.text_code[2:] == '----':
      # text_code is 'Lprocg' for all problems.  Any letter may be replaced with '-'
      # '----' in last 4 slots indicates no validation problems
      return False
    else:
      return True

  def as_text(self):
    return (' '.join([self.chain, self.resseq, self.text_code, self.feedback]))


class barbed_wire_analysis():
  def __init__(self, model):
    self.res_dict = {}  # keyed for easy lookup
    self.res_list = []  # ordered for finding sequence-related residues
    self.chunk_list = []

    hierarchy = model.get_hierarchy()
    self.load_residues(hierarchy)
    ##print('\n'.join(bwa.res_dict.keys()))

    self.add_secondary_structure(model, hierarchy)
    self.add_contacts(model)
    self.analyze_contacts()
    self.add_ramalyze(hierarchy)
    self.add_cablam(hierarchy)
    self.add_omegalyze(hierarchy)
    self.add_covalent_geometry(model)

    self.predictalyze()
    self.merge_analyses()
    self.merge_small_chunks()

  def load_residues(self, hierarchy):
    for chain in hierarchy.chains():
      for rg in chain.residue_groups():
        r = predicted_residue(rg)
        self.res_dict[r.resid] = r
        self.res_list.append(r)

  def add_secondary_structure(self, model, hierarchy):
    from mmtbx.secondary_structure import sec_str_master_phil_str
    from mmtbx.secondary_structure import manager as ss_manager
    asc = model.get_atom_selection_cache()
    sec_str_master_phil = iotbx.phil.parse(sec_str_master_phil_str)
    ss_params = sec_str_master_phil.fetch().extract()
    ss_params.secondary_structure.protein.search_method = "ksdssp"
    ssm = ss_manager(hierarchy,
                     atom_selection_cache=asc,
                     geometry_restraints_manager=None,
                     sec_str_from_pdb_file=None,
                     # params=None,
                     params=ss_params.secondary_structure,
                     was_initialized=False,
                     mon_lib_srv=None,
                     verbose=-1,
                     log=null_out(),
                     # log=sys.stdout,
                     )
    secstr = ssm.find_sec_str(hierarchy)
    for helix in secstr.helices:
      helix_id = helix.helix_id
      start_chain_id = helix.start_chain_id.strip()
      start_resseq = helix.start_resseq
      start_resid = ','.join([start_chain_id, start_resseq])
      r = self.res_dict[start_resid]
      i = self.res_list.index(r)
      res_range = self.res_list[i: i + helix.length]
      for r in res_range:
        r.ss_type = "helix"
        r.ss_index = helix_id

    for sheet in secstr.sheets:
      for strand in sheet.strands:
        sheet_id = strand.sheet_id
        # strand_id = strand.strand_id
        start_chain_id = strand.start_chain_id.strip()
        start_resseq = strand.start_resseq
        end_resseq = strand.end_resseq
        start_resid = ','.join([start_chain_id, start_resseq])
        r = self.res_dict[start_resid]
        i = self.res_list.index(r)
        strand_length = int(end_resseq.strip()) - int(start_resseq.strip()) + 1
        res_range = self.res_list[i: i + strand_length]
        for r in res_range:
          r.ss_type = "sheet"
          r.ss_index = sheet_id

  def are_same_ss_element(self, resid1, resid2):
    res1 = self.res_dict[resid1]
    res2 = self.res_dict[resid2]
    if not (res1.ss_type and res2.ss_type):
      # at least one coil - only sequence proximity matters for coil
      return False
    if not res1.ss_type == res2.ss_type:
      # different secondary structure types cannot be same element
      return False
    if not res1.ss_index == res2.ss_index:
      # same index means same element if same type (HELIX and SHEET can share indices)
      return False
    # same type and same index
    return True

  def add_ramalyze(self, hierarchy):
    from mmtbx.validation import ramalyze
    rama_results = ramalyze.ramalyze(pdb_hierarchy=hierarchy,
                                     outliers_only=False)
    for result in rama_results.results:
      resid = ','.join([result.chain_id, result.resseq])
      if result.is_outlier():
        self.res_dict[resid].out_rama = True
      if 80.0 < result.psi < 170.0:
        self.res_dict[resid].rama_high_psi = True
      if 20.0 < result.phi < 130.0:  # this region is indicative of general-case residue types
        self.res_dict[resid].rama_high_phi = True

  def add_omegalyze(self, hierarchy):
    from mmtbx.validation import omegalyze
    omega_results = omegalyze.omegalyze(pdb_hierarchy=hierarchy,
                                        nontrans_only=True)
    for result in omega_results.results:
      resid = ','.join([result.chain_id, result.resseq])
      if result.resname == "PRO":
        if result.omegalyze_type == "Twisted":  # (cis Pro is okay by default)
          self.res_dict[resid].out_omega = True
      else:
        if result.is_outlier():  # should always be true with nontrans_only=True run
          self.res_dict[resid].out_omega = True

  def add_contacts(self, model):
    contacts = do_probe(model)
    for c in contacts.values():
      if self.are_same_ss_element(c.srcResid, c.trgResid):
        # only interested in packing between different elements
        # will have to develop better method for large beta, eventually
        continue
      r = self.res_dict[c.srcResid]
      if abs(r.resnum - self.res_dict[c.trgResid].resnum) <= 4:
        # sequence distance between src and trg
        # local contacts do not count for packing
        # sequence dist cutoff subject to change
        continue
      if c.interactiontype == "bo":
        r.packing_bo.append([c.contactId])
        if float(c.mingap) >= 0.9:
          r.severe_clashes += 1
      elif c.interactiontype == "hb":
        r.packing_hb.append([c.contactId])
      else:  # so, cc, wc
        r.packing_vdw.append([c.contactId])

  def analyze_contacts(self):
    i, j = 0, 5  # window of 5 res
    while j < len(self.res_list):
      res_slice = self.res_list[i:j]
      total_packing = sum([r.packing_contact_count() for r in res_slice])
      total_heavy_atoms = sum([r.heavy_atom_count for r in res_slice])
      packing_ratio = total_packing / total_heavy_atoms
      r = res_slice[2]
      if r.ss_type == "sheet":
        if packing_ratio <= 0.1:
          r.packing_quality = 0
        elif packing_ratio <= 0.35:  # different for sheet
          r.packing_quality = 1
        elif packing_ratio <= 1:
          r.packing_quality = 2
        else:
          r.packing_quality = 3
      else:  # helix and coil
        if packing_ratio <= 0.1:
          r.packing_quality = 0
        elif packing_ratio <= 0.6:  # general case (helix/coil)
          r.packing_quality = 1
        elif packing_ratio <= 1:
          r.packing_quality = 2
        else:
          r.packing_quality = 3
        pass
      # print(res_slice[2].resseq, res_slice[2].ss_type, res_slice[2].ss_index, "%.3f" % (total_packing/total_heavy_atoms))#, res_slice[2].packing_vdw)
      i += 1
      j += 1

  def add_covalent_geometry(self, model):
    from mmtbx.validation.mp_validate_bonds import mp_angles
    from mmtbx.model import manager
    from libtbx.utils import null_out
    mc_atoms = [" CA ", " N  ", " C  ", " O  "]

    model.set_stop_for_unknowns(False)
    hierarchy = model.get_hierarchy()
    p = manager.get_default_pdb_interpretation_params()
    p.pdb_interpretation.allow_polymer_cross_special_position = True
    p.pdb_interpretation.flip_symmetric_amino_acids = False
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
    model.set_log(log = null_out())
    model.process(make_restraints=True, pdb_interpretation_params=p)
    geometry = model.get_restraints_manager().geometry
    atoms = hierarchy.atoms()
    angles = mp_angles(
      pdb_hierarchy=hierarchy,
      pdb_atoms=atoms,
      geometry_restraints_manager=geometry,
      outliers_only=True)
    for angle in angles.results:
      resid = ','.join([angle.atoms_info[1].chain_id, angle.atoms_info[1].resseq])
      atomnames = [angle.atoms_info[0].name, angle.atoms_info[1].name, angle.atoms_info[2].name]
      if atomnames[0] not in mc_atoms: continue
      if atomnames[1] not in mc_atoms: continue
      if atomnames[2] not in mc_atoms: continue
      self.res_dict[resid].out_geom = True

  def add_cablam(self, hierarchy):
    from mmtbx.validation import cablam
    cablam_results = cablam.cablamalyze(pdb_hierarchy=hierarchy,
                                        outliers_only=True,
                                        out=sys.stdout,
                                        quiet=True)
    for result in cablam_results.results:
      resid = ','.join([result.chain_id.strip(), result.resseq])
      if result.feedback.cablam_outlier or result.feedback.c_alpha_geom_outlier:  # ignore cablam_disfavored
        self.res_dict[resid].out_cablam = True

  def has_barbed_wire_errors(self, res_slice):
    text_code = ['-'] * 6  # Lproag
    rama = 0
    omega = 0
    cablam = 0
    geom = 0
    psi = 0
    for r in res_slice:
      if r.out_rama: rama += 1
      if r.out_cablam: cablam += 1
      if r.out_omega: omega += 1
      if r.out_geom: geom += 1
      if r.rama_high_psi: psi += 1
    score = 0
    if rama >= 1 and psi == 3:  # all in high-psi region, plus at least 1 outlier
      score += 1
      text_code[2] = 'r'
    if omega >= 2:
      score += 1
      text_code[3] = 'o'
    if cablam >= 2:
      score += 1
      text_code[4] = 'c'
    if geom >= 2:
      score += 1
      text_code[5] = 'g'
    if score >= 2:
      is_barbed_like = True
    else:
      is_barbed_like = False
    return is_barbed_like, text_code

  def predictalyze(self):
    i, j = 0, 3  # window of 3 res
    while j < len(self.res_list):
      res_slice = self.res_list[i:j]
      r = res_slice[1]
      if r.packing_quality is None:
        i += 1
        j += 1
        continue
      is_barbed_like, text_code = self.has_barbed_wire_errors(res_slice)

      if r.plddt >= 70:
        if r.packing_quality >= 1:
          r.feedback = "Predictive"
        else:
          text_code[1] = 'p'
          r.feedback = "Unpacked high pLDDT"
      else:  # plddt should get subdivided based on further analysis
        text_code[0] = 'L'
        if r.packing_quality >= 1:
          if is_barbed_like:
            r.feedback = "Unphysical"
          else:
            r.feedback = "Near-predictive"
        else:  # poor packing
          text_code[1] = 'p'
          if is_barbed_like:
            r.feedback = "Barbed wire"
          else:
            r.feedback = "Unpacked possible"
      r.text_code = ''.join(text_code)
      i += 1
      j += 1

  def merge_analyses(self):
    chunk_list = []
    prev_r = None
    current_chunk = []
    for r in self.res_list:
      if prev_r is None:
        current_chunk.append(r.resid)
        prev_r = r
        continue

      # TODO: also start new chunk if chain is new

      if r.feedback == prev_r.feedback:
        # if same type, extend current chunk
        current_chunk.append(r.resid)
      else:
        # if different, end this chunk and start a new one with this residue
        chunk_list.append(alphafold_chunk(current_chunk, prev_r.feedback))
        current_chunk = [r.resid]
      prev_r = r
    chunk_list.append(alphafold_chunk(current_chunk, r.feedback))
    self.chunk_list = chunk_list

  def merge_small_chunks(self):
    i = 0
    while i < len(self.chunk_list):
      c = self.chunk_list[i]
      if (i - 1 < 0) or (i + 1 >= len(self.chunk_list)):
        i += 1
        continue
      prev_c = self.chunk_list[i - 1]
      next_c = self.chunk_list[i + 1]
      if len(c.members) > 2 or c.prediction_type != "Unpacked possible":
        i += 1
        continue
      if prev_c.prediction_type == "Barbed wire":
        for r in c.members:
          if self.res_dict[r].has_any_validation_errors():
            prev_c.add_to_end(c.remove_from_start())
          else:
            break
      if len(c.members) == 0:
        i += 1
        continue
      else:
        if next_c.prediction_type == "Barbed wire":
          for r in reversed(c.members):
            if self.res_dict[r].has_any_validation_errors():
              next_c.add_to_start(c.remove_from_end())
            else:
              break
      i += 1
    self.remove_empty_chunks()
    self.merge_similar_chunks()
    self.remove_empty_chunks()

  def remove_empty_chunks(self):
    i = len(self.chunk_list) - 1
    while i >= 0:
      c = self.chunk_list[i]
      if len(c.members) == 0:
        self.chunk_list.pop(i)
      i -= 1

  def merge_similar_chunks(self):
    i = 0
    while i < len(self.chunk_list):
      c = self.chunk_list[i]
      if (i - 1 < 0):
        i += 1
        continue
      prev_c = self.chunk_list[i - 1]
      if c.prediction_type == prev_c.prediction_type:
        c.members = prev_c.members + c.members
        prev_c.members = []
      i += 1

  # ----------------------OUTPUT---------------------
  def as_text_residues(self, out=sys.stdout):
    for r in self.res_list:
      print(r.as_text(), file=out)

  def as_text_chunks(self, out=sys.stdout):
    for c in self.chunk_list:
      print(c.start, "to", c.end, c.prediction_type, len(c.members), file=out)

  def as_kinemage(self, out=sys.stdout):
    # colored ball at each CA, color based on current synthesis
    # label with bc--go- style text showing components of decision
    balls = []
    labels = []
    for r in self.res_list:
      if not r.feedback:
        continue
      ball_color = None
      if r.feedback == "Predictive":
        ball_color = "sky"
      elif r.feedback == "Unphysical":
        ball_color = "purple"
      elif r.feedback == "Barbed wire":
        ball_color = "hotpink"
      elif r.feedback == "Unpacked possible":
        ball_color = "gold"
      elif r.feedback == "Near-predictive":
        ball_color = 'green'
      elif r.feedback == "Unpacked high pLDDT":
        ball_color = 'gray'
      else:
        ball_color = 'brown'
      ballline = "{%s %s}%s %.3f %.3f %.3f" % (r.resid, r.text_code, ball_color, r.caxyz[0], r.caxyz[1], r.caxyz[2])
      labelline = "{  %s}%s %.3f %.3f %.3f" % (r.text_code, ball_color, r.caxyz[0], r.caxyz[1], r.caxyz[2])
      balls.append(ballline)
      labels.append(labelline)
    print("@group {bwa markup}", file=out)
    print("@balllist {bwa balls} radius=0.6", file=out)
    for line in balls:
      print(line, file=out)
    print("@labellist {bwa labels}", file=out)
    for line in labels:
      print(line, file=out)

  def as_selection_string(self):
    selection_list = []
    for c in self.chunk_list:
      if c.prediction_type in ["Predictive","Near-predictive"]:
        selection_list.append(c.as_selection_string())
    return " or ".join(selection_list)

  def as_selection_file(self, model, out=sys.stdout):
    hierarchy = model.get_hierarchy()
    atom_selection = self.as_selection_string()
    sele = hierarchy.apply_atom_selection(atom_selection)
    print(sele.as_pdb_string(), file=out)

  # ----------------------END OUTPUT---------------------
