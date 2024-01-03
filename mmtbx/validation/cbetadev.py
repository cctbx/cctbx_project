
"""
Validation of protein geometry by analysis of the positions of C-beta sidechain
atoms.  Significant deviations from ideality often indicate misfit rotamers
and/or necessity of mainchain movement, especially alternate conformations.

Reference:

Lovell SC, Davis IW, Arendall WB 3rd, de Bakker PI, Word JM, Prisant MG,
Richardson JS, Richardson DC.  Structure validation by Calpha geometry: phi,psi
and Cbeta deviation.  Proteins. 2003 Feb 15;50(3):437-50.
http://www.ncbi.nlm.nih.gov/pubmed/12557186
"""

from __future__ import absolute_import, division, print_function
from mmtbx.validation import residue, validation
from scitbx.matrix import col, dihedral_angle, rotate_point_around_axis
import sys
from scitbx.array_family import flex
from libtbx import group_args

from mmtbx.conformation_dependent_library import generate_protein_threes
from mmtbx.validation import cbd_utils
import json

relevant_atom_names = {
  " CA ": None, " N  ": None, " C  ": None, " CB ": None} # FUTURE: set

def get_phi_psi_dict(pdb_hierarchy):
  rc = {}
  for i, three in enumerate(generate_protein_threes(hierarchy=pdb_hierarchy,
                                                    geometry=None)):
    phi_psi_angles = three.get_phi_psi_angles()
    is_alt_conf = ' '
    relevant_atoms = {}
    for atom in three[1].atoms():
      if (atom.name in relevant_atom_names):
        if (len(atom.parent().altloc) != 0):
          is_alt_conf = atom.parent().altloc
          break
    id_str = '|%s:%s|' % (three[1].id_str(), is_alt_conf)
    rc[id_str] = phi_psi_angles
  return rc

class cbeta(residue):
  """
  Result class for protein C-beta deviation analysis (phenix.cbetadev).
  """
  __cbeta_attr__ = [
    "deviation",
    "dihedral_NABB",
    "ideal_xyz",
    "model_id"
  ]
  __slots__ = residue.__slots__ + __cbeta_attr__

  @staticmethod
  def header():
    return "%-20s  %5s" % ("Residue", "Dev.")

  def as_string(self):
    return "%-20s  %5.2f" % (self.id_str(), self.deviation)

  # Backwards compatibility
  def format_old(self):
    return "%s:%s:%2s:%4s%1s:%7.3f:%7.2f:%7.2f:%s:" % (self.altloc,
      self.resname.lower(), self.chain_id, self.resseq, self.icode,
      self.deviation, self.dihedral_NABB, self.occupancy, self.altloc)

  def as_kinemage(self):
    key = "cb %3s%2s%4s%1s  %.3f %.2f" % (self.resname.lower(),
      self.chain_id, self.resseq, self.icode, self.deviation,
      self.dihedral_NABB)
    return "{%s} r=%.3f %s %.3f, %.3f, %.3f" % (key,
      self.deviation, '', #The blank char is a placeholder for optional color
      self.ideal_xyz[0], self.ideal_xyz[1], self.ideal_xyz[2])

  def as_bullseye_point(self):
    #print a point in kinemage format for the "bulleye" plot of cbdev distribution
    import numpy as np
    #original point id format: {4hum  trp A 427; dev=0.090}
    key = "%s%3s%2s%4s%1s;  dev=%.3f" % (self.altloc.strip(), self.resname.lower(),
      self.chain_id, self.resseq, self.icode, self.deviation)
    #convert polar position to cartesian
    angle = np.radians(self.dihedral_NABB)
    x = self.deviation * np.cos(angle)
    y = self.deviation * np.sin(angle)
    return "{%s} %.3f %.3f 0" % (key,x,y)

  def as_bullseye_label(self):
    #label a point in the cbdev bulleye, intended for outliers
    import numpy as np
    #expected label format: { 4hum  trp A 427}
    #shorter than full point id, leading space to offset from point
    key = "  %s%3s%2s%4s%1s" % (self.altloc.strip(), self.resname.lower(),
      self.chain_id, self.resseq, self.icode)
    #convert polar position to cartesian
    angle = np.radians(self.dihedral_NABB)
    x = self.deviation * np.cos(angle)
    y = self.deviation * np.sin(angle)
    return "{%s} %.3f %.3f 0" % (key,x,y)

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    return json.dumps(slots_as_dict, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_table_row_phenix(self):
    return [ self.chain_id, "%1s%s %s" % (self.altloc,self.resname, self.resid),
             self.deviation, self.dihedral_NABB ]

class cbetadev(validation):
  __slots__ = validation.__slots__ + ["beta_ideal","_outlier_i_seqs","stats",
                                      'percent_outliers', 'new_outliers', 'outliers_removed']
  program_description = "Analyze protein sidechain C-beta deviation"
  output_header = "pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:"
  gui_list_headers = ["Chain", "Residue","Deviation","Angle"]
  gui_formats = ["%s", "%s", "%.3f", "%.2f"]
  wx_column_widths = [75, 125, 100, 100]

  def get_result_class(self) : return cbeta

  def __init__(self, pdb_hierarchy,
      outliers_only=False,
      out=sys.stdout,
      collect_ideal=False,
      apply_phi_psi_correction=False,
      display_phi_psi_correction=False,
      exclude_d_peptides=False,
      quiet=False):
    validation.__init__(self)
    self._outlier_i_seqs = flex.size_t()
    self.beta_ideal = {}
    output_list = []
    self.stats = group_args(n_results=0,
                            n_weighted_results = 0,
                            n_weighted_outliers = 0)
    total_residues = 0
    new_outliers = None
    outliers_removed = None
    if apply_phi_psi_correction:
      phi_psi_angles = get_phi_psi_dict(pdb_hierarchy)
      new_outliers = 0
      outliers_removed = 0
    from mmtbx.validation import utils
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)
    for model in pdb_hierarchy.models():
      if model.id not in self.n_total_by_model:
        self.n_total_by_model[model.id] = 0
        self.n_outliers_by_model[model.id] = 0
      for chain in model.chains():
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain)
        else:
          chain_id = chain.id
        for rg in chain.residue_groups():
          for i_cf,cf in enumerate(rg.conformers()):
            for i_residue,residue in enumerate(cf.residues()):
              if (residue.resname == "GLY"):
                continue
              is_first = (i_cf == 0)
              is_alt_conf = False
              relevant_atoms = {}
              for atom in residue.atoms():
                if (atom.name in relevant_atom_names):
                  relevant_atoms[atom.name] = atom
                  if (len(atom.parent().altloc) != 0):
                    is_alt_conf = True
              if ((is_first or is_alt_conf) and len(relevant_atoms) == 4):
                result = calculate_ideal_and_deviation(
                  relevant_atoms=relevant_atoms,
                  resname=residue.resname)
                dev = result.deviation
                dihedralNABB = result.dihedral
                betaxyz = result.ideal
                if (dev is None) : continue
                resCB = relevant_atoms[" CB "]
                self.stats.n_results += 1
                self.n_total_by_model[model.id] += 1
                self.stats.n_weighted_results += resCB.occ
                if (is_alt_conf):
                  altchar = cf.altloc
                else:
                  altchar = " "
                total_residues+=1
                if apply_phi_psi_correction:
                  id_str = '|%s:%s|' % (residue.id_str(), altchar)
                  phi_psi = phi_psi_angles.get(id_str, None)
                  if phi_psi:
                    rc = cbd_utils.get_phi_psi_correction(
                      result,
                      residue,
                      phi_psi,
                      display_phi_psi_correction=display_phi_psi_correction,
                      )
                    if rc:
                      dev, dihedralNABB, start, finish = rc
                      # if start or finish: print(id_str,dev,dihedralNABB,start,finish)
                      if start and not finish:
                        outliers_removed += 1
                      elif not start and finish:
                        new_outliers += 1
                if(exclude_d_peptides and dev>=2.):
                  pass
                elif(dev >=0.25 or outliers_only==False):
                  if(dev >=0.25):
                    self.n_outliers+=1
                    self.n_outliers_by_model[model.id]+=1
                    self.stats.n_weighted_outliers += resCB.occ
                    self._outlier_i_seqs.append(atom.i_seq)
                  res=residue.resname.lower()
                  sub=chain.id
                  if(len(sub)==1):
                    sub=" "+sub
                  result = cbeta(
                    model_id=model.id,
                    chain_id=chain_id,
                    resname=residue.resname,
                    resseq=residue.resseq,
                    icode=residue.icode,
                    altloc=altchar,
                    xyz=resCB.xyz,
                    occupancy=resCB.occ,
                    deviation=dev,
                    dihedral_NABB=dihedralNABB,
                    ideal_xyz=betaxyz,
                    outlier=(dev >= 0.25))
                  self.results.append(result)
                  key = result.id_str()
                  if (collect_ideal):
                    self.beta_ideal[key] = betaxyz
      if apply_phi_psi_correction:
        print('''
  Outliers removed : %5d
  New outliers     : %5d
  Num. of outliers : %5d
  Num. of residues : %5d
  ''' % (outliers_removed,
         new_outliers,
         self.n_outliers,
         total_residues,
        ))
    self.new_outliers=new_outliers
    self.outliers_removed=outliers_removed
    if total_residues:
      self.percent_outliers=self.n_outliers/total_residues*100
    else:
      self.percent_outliers = None

  def show_old_output(self, out, verbose=False, prefix="pdb"):
    if (verbose):
      print(self.output_header, file=out)
    for result in self.results :
      print(prefix + " :" + result.format_old(), file=out)
    if (verbose):
      self.show_summary(out)

  def show_summary(self, out, prefix=""):
    print(prefix + \
      'SUMMARY: %d C-beta deviations >= 0.25 Angstrom (Goal: 0)' % \
      self.n_outliers, file=out)

  #functions for internal access of summary statistics
  def get_outlier_count(self):
    return self.n_outliers

  def get_weighted_outlier_count(self):
    return self.stats.n_weighted_outliers

  def get_result_count(self):
    return self.stats.n_results

  def get_weighted_result_count(self):
    return self.stats.n_weighted_results

  def get_outlier_percent(self):
    if self.stats.n_results == 0:
      return 0
    return self.n_outliers/self.stats.n_results*100

  def get_weighted_outlier_percent(self):
    weighted_result_count = self.get_weighted_result_count()
    if weighted_result_count == 0:
      return 0
    return self.get_weighted_outlier_count()/weighted_result_count*100

  def get_expected_count(self):
    return 0

  def get_beta_ideal(self):
    return self.beta_ideal

  def as_kinemage(self, chain_id=None):
    cbeta_out = "@subgroup {CB dev} dominant\n"
    cbeta_out += "@balllist {CB dev Ball} color= magenta master= {Cbeta dev}\n"
    for result in self.results :
      if result.is_outlier():
        if (chain_id is None) or (chain_id == result.chain_id):
          cbeta_out += result.as_kinemage() + "\n"
    return cbeta_out

  def as_bullseye_kinemage(self, pdbid=""):
    from mmtbx.validation.molprobity import kinemage_templates
    header = []
    header.append(kinemage_templates.cbetadev_bullseye(pdbid=pdbid))
    header.append("@group {Cbeta dev}")
    cbeta_main = ["@dotlist {Cb scatter} color= white"]
    cbeta_main_labels = ["@labellist {outlier labels} color= white"]
    cbeta_alt = ["@dotlist {alt conf Cb scatter} color= pink"]
    cbeta_alt_labels = ["@labellist {alt conf outlier labels} color= pinktint"]

    #cbeta_main.append("@dotlist {Cb scatter} color= white")
    for result in self.results:
      if result.altloc in ['',' ','A']:
        cbeta_main.append(result.as_bullseye_point())
        if result.is_outlier():
          cbeta_main_labels.append(result.as_bullseye_label())
      else:
        cbeta_alt.append(result.as_bullseye_point())
        if result.is_outlier():
          cbeta_alt_labels.append(result.as_bullseye_label())
    #if no residues were added to any of the extra lists, also empty the header for that list
    if len(cbeta_main_labels) == 1: cbeta_main_labels = []
    if len(cbeta_alt) == 1: cbeta_alt = []
    if len(cbeta_alt_labels) == 1: cbeta_alt_labels = []
    return "\n".join(header + cbeta_main + cbeta_main_labels + cbeta_alt + cbeta_alt_labels)

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "cbetadev"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results
    data['summary_results'] = summary_results
    for model_id in self.n_total_by_model.keys():
      if self.n_total_by_model[model_id]:
        summary_results[model_id] = {
          "num_outliers" : self.n_outliers_by_model[model_id],
          "num_cbeta_residues" : self.n_total_by_model[model_id],
          "outlier_percentage" : self.n_outliers_by_model[model_id]/self.n_total_by_model[model_id]*100,
          "outlier_goal" : 0,
        }
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def as_coot_data(self):
    data = []
    for result in self.results :
      if result.is_outlier():
        data.append((result.chain_id, result.resid, result.resname,
          result.altloc, result.deviation, result.xyz))
    return data

class calculate_ideal_and_deviation(object):
  __slots__ = ["deviation", "ideal", "dihedral"]
  def __init__(self, relevant_atoms, resname):
    assert (resname != "GLY")
    from cctbx.geometry_restraints import chirality
    self.deviation = None
    self.ideal = None
    self.dihedral = None
    resCA = relevant_atoms[" CA "]
    resN  = relevant_atoms[" N  "]
    resC  = relevant_atoms[" C  "]
    resCB = relevant_atoms[" CB "]
    if None not in [resCA, resN, resCB, resC]:
      chiral_volume = chirality([resCA.xyz, resN.xyz, resCB.xyz, resC.xyz],
                               volume_ideal=0.,
                               both_signs=True,
                               weight=1.,
                               ).volume_model
    else:
      chiral_volume = None
    dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal= \
      idealized_calpha_angles(resname, chiral_volume)
    betaNCAB = construct_fourth(resN,
                                resCA,
                                resC,
                                dist,
                                angleCAB,
                                dihedralNCAB,
                                method="NCAB")
    betaCNAB = construct_fourth(resN,
                                resCA,
                                resC,
                                dist,
                                angleNAB,
                                dihedralCNAB,
                                method="CNAB")
    if (not None in [betaNCAB, betaCNAB]):
      betaxyz = (col(betaNCAB) + col(betaCNAB)) / 2
      betadist = abs(col(resCA.xyz) - betaxyz)
      if betadist != 0:
        if(betadist != dist):
          distTemp = betaxyz - col(resCA.xyz)
          betaxyz = col(resCA.xyz) + distTemp * dist/betadist
        self.deviation = abs(col(resCB.xyz) - betaxyz)
        self.dihedral = dihedral_angle(
          sites=[resN.xyz,resCA.xyz,betaxyz.elems,resCB.xyz], deg=True)
        self.ideal = betaxyz.elems

def idealized_calpha_angles(resname, chiral_volume=None):
  from iotbx.pdb import common_residue_names_get_class
  if (resname in ["ALA", "DAL"]):
    #target values are for L-aminao acids.
    #D-amino acids will require the signs of the dihedrals to be flipped
    #This will be performed just before the return
    dist = 1.536
    angleCAB = 110.1
    dihedralNCAB = 122.9
    angleNAB = 110.6
    dihedralCNAB = -122.6
    angleideal = 111.2
  elif (resname in ["PRO", "DPR"]):
    dist = 1.530
    angleCAB = 112.2
    dihedralNCAB = 115.1
    angleNAB = 103.0
    dihedralCNAB = -120.7
    angleideal = 111.8
  elif (resname in ["VAL", "THR", "ILE",
                    "DVA", "DTH", "DIL"]):
    dist = 1.540
    angleCAB = 109.1
    dihedralNCAB = 123.4
    angleNAB = 111.5
    dihedralCNAB = -122.0
    angleideal = 111.2
  elif (resname == "GLY"):
    dist = 1.10
    angleCAB = 109.3
    dihedralNCAB = 121.6
    angleNAB = 109.3
    dihedralCNAB = -121.6
    angleideal = 112.5
  else:
    dist = 1.530
    angleCAB = 110.1
    dihedralNCAB = 122.8
    angleNAB = 110.5
    dihedralCNAB = -122.6
    angleideal = 111.2
  #check if dihedral signs should be flipped becasue residue has D chirality
  #rely on residue names first, so that mis-named residues gat marked as outliers
  #if the residue name is nonstandard, use the chiral volume for best guess about target
  res_class = common_residue_names_get_class(resname)
  if res_class == "common_amino_acid" or chiral_volume is None:
    pass
  elif res_class == "d_amino_acid" or chiral_volume > 0:
    dihedralNCAB = dihedralNCAB * -1
    dihedralCNAB = dihedralCNAB * -1
  return dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal

def construct_fourth(resN,resCA,resC,dist,angle,dihedral,method="NCAB"):
  if (not None in [resN, resCA, resC]):
    if (method == "NCAB"):
      res0 = resN
      res1 = resC
      res2 = resCA
    elif (method == "CNAB"):
      res0 = resC
      res1 = resN
      res2 = resCA
    a = col(res2.xyz) - col(res1.xyz)
    b = col(res0.xyz) - col(res1.xyz)
    c = a.cross(b)
    cmag = abs(c)
    if(cmag > 0.000001):
      c *= dist/cmag
    c += col(res2.xyz)
    d = c
    angledhdrl = dihedral - 90
    a = col(res1.xyz)
    b = col(res2.xyz)
    # XXX is there an equivalent method for 'col'?
    newD = col(rotate_point_around_axis(
      axis_point_1=res1.xyz,
      axis_point_2=res2.xyz,
      point=d.elems,
      angle=angledhdrl,
      deg=True))
    a = newD - col(res2.xyz)
    b = col(res1.xyz) - col(res2.xyz)
    c = a.cross(b)
    cmag = abs(c)
    if(cmag > 0.000001):
      c *= dist/cmag
    angledhdrl = 90 - angle;
    a = col(res2.xyz)
    c += a
    b = c
    if a == b: return newD.elems
    return rotate_point_around_axis(
      axis_point_1=a.elems,
      axis_point_2=b.elems,
      point=newD.elems,
      angle=angledhdrl,
      deg=True)
  return None

def extract_atoms_from_residue_group(residue_group):
  """
  Given a residue_group object, which may or may not have multiple
  conformations, extract the relevant atoms for each conformer, taking into
  account any atoms shared between conformers.  This is implemented
  separately from the main validation routine, which accesses the hierarchy
  object via the chain->conformer->residue API.  Returns a list of hashes,
  each suitable for calling calculate_ideal_and_deviation.
  """
  atom_groups = residue_group.atom_groups()
  if (len(atom_groups) == 1):
    relevant_atoms = {}
    for atom in atom_groups[0].atoms():
      relevant_atoms[atom.name] = atom
    return [ relevant_atoms ]
  else :
    all_relevant_atoms = []
    expected_names = [" CA ", " N  ", " CB ", " C  "]
    main_conf = {}
    for atom_group in atom_groups :
      if (atom_group.altloc.strip() == ''):
        for atom in atom_group.atoms():
          if (atom.name in expected_names):
            main_conf[atom.name] = atom
      else :
        relevant_atoms = {}
        for atom in atom_group.atoms():
          if (atom.name in expected_names):
            relevant_atoms[atom.name] = atom
        if (len(relevant_atoms) == 0) : continue
        for atom_name in main_conf.keys():
          if (not atom_name in relevant_atoms):
            relevant_atoms[atom_name] = main_conf[atom_name]
        if (len(relevant_atoms) != 0):
          all_relevant_atoms.append(relevant_atoms)
    if (len(main_conf) == 4):
      all_relevant_atoms.insert(0, main_conf)
    return all_relevant_atoms
