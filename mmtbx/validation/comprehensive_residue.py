from __future__ import absolute_import, division, print_function
import iotbx.phil
from mmtbx.validation import ramalyze
from mmtbx.validation import rotalyze
import boost_adaptbx.boost.python as bp
import sys
from six.moves import range

ext = bp.import_ext("mmtbx_ramachandran_restraints_ext")

def get_master_phil():
  return iotbx.phil.parse(input_string="""
    model = None
      .type = path
      .help = "Model file (PDB or mmCIF)"
    verbose = True
      .type = bool
      .help = "Be verbose"
    help = False
      .type = bool
      .help = "Prints this help message if true"
""", process_includes=True)

usage_string = """
To be determined.
"""

def construct_complete_sidechain(residue,altloc):
  if residue is None : return {}
  complete_dict = {}
  atom_dict = {}
  #print residue.resname
  if not rotalyze.has_heavy_atoms(residue.atoms()) : return {}
  for atom in residue.atoms():
    #handle hydrogen/deuterium swaps
    #print atom.name
    if atom_dict.get(atom.name) == None:
      if atom_dict.get(atom.name.replace("H","D",1)) != None:
        del(atom_dict[atom.name.replace("H","D",1)])
      elif atom_dict.get(atom.name.replace("D","H",1)) != None:
        del(atom_dict[atom.name.replace("D","H",1)])
      atom_dict[atom.name] = atom
  clone_dict = {}
  clone_dict.update(atom_dict)
  complete_dict[altloc] = clone_dict
  if len(complete_dict) > 0:
    return complete_dict
  return {}

class ValidationResidue(object):

  def __init__(self,three,rama_eval,rota_eval,rotamer_id,index=1):
    self.three = three
    self.rama_eval = rama_eval
    self.rota_eval = rota_eval
    self.rotamer_id = rotamer_id
    self.index = index
    self.chain_id = self.three[index].parent().parent().id
    self.resseq = self.three[index].resseq
    self.icode = self.three[index].icode
    self.resname = self.three[index].resname.upper()
    if self.three[index].is_pure_main_conf : self.altloc = ' '
    else : self.altloc = self.three[index].parent().altloc
    # when index is 1 we can do rama if not we can't
    # when index is 0 we can't do omega
    assert index in range(3)

    #print self.three[index].resseq
    phi_psi_atoms = self.three.get_phi_psi_atoms()
    if phi_psi_atoms :
      # FIXME you are adding non alts twice because threes method.
      self.set_rama_result(phi_psi_atoms)
      print(self.rama_result)
    self.set_rota_result()
    if self.resname not in ['GLY','ALA'] : print(self.rota_result)
    if self.resname != "GLY" :
      self.set_cbeta_result()
      print(self.cbeta_result)
    if self.index >= 1 :
      self.set_omega_result()
      print(self.omega_result.as_string())

  def get_start_kwargs(self):
    return {'chain_id':self.chain_id,
            'resseq':self.resseq,
            'icode':self.icode,
            'resname':self.resname,
            'altloc':self.altloc}

  def set_omega_result(self):
    from mmtbx.validation import omegalyze
    kwargs = self.get_start_kwargs()
    if self.index == 1 : omega_return = 'middle'
    elif self.index == 2 : omega_return = 'last'
    else : raise RuntimeError('unexpected index')
    omega = self.three.get_omega_value(omega_return=omega_return)
    if self.resname == "PRO" : res_type = omegalyze.OMEGA_PRO
    else : res_type = omegalyze.OMEGA_GENERAL
    omega_type = omegalyze.find_omega_type(omega)
    is_nontrans = False
    if omega_type in [omegalyze.OMEGALYZE_CIS, omegalyze.OMEGALYZE_TWISTED] :
      is_nontrans = True
    highest_mc_b = 0
    bb = [' N  ',' CA ',' CB ',' C  ',' O  ']
    atms = [a for a in self.three[self.index].atoms() if a.name not in bb]
    atms+= [a for a in self.three[self.index-1].atoms() if a.name not in bb]
    for a in atms :
      if a.b > highest_mc_b : highest_mc_b = a.b
    kwargs.update({'prev_resseq':self.three[self.index-1].resseq,
                   'prev_icode':self.three[self.index-1].icode,
                   'prev_resname':self.three[self.index-1].resname,
                   'prev_altloc':self.altloc,
                   'omega':omega,
                   'omega_type':omega_type,
                   'res_type':res_type,
                   'is_nontrans':is_nontrans,
                   'highest_mc_b':highest_mc_b})
    self.omega_result = omegalyze.omega_result(**kwargs)

  def set_cbeta_result(self):
    from mmtbx.validation import cbetadev
    kwargs = self.get_start_kwargs()
    # get relevant_atoms
    relevant_atoms = {}
    for atom in self.three[self.index].atoms():
      if (atom.name in [" CA ", " N  ", " C  ", " CB "]):
        relevant_atoms[atom.name] = atom
    result = cbetadev.calculate_ideal_and_deviation(
      relevant_atoms=relevant_atoms,
      resname=self.resname)
    dev = result.deviation
    dihedralNABB = result.dihedral
    betaxyz = result.ideal
    if (dev is None) : return
    res=self.resname.lower()
    resCB = relevant_atoms[" CB "]
    kwargs.update({'xyz':resCB.xyz,
                  'occupancy':resCB.occ,
                  'deviation':dev,
                  'dihedral_NABB':dihedralNABB,
                  'ideal_xyz':betaxyz,
                  'outlier':(dev >= 0.25)})
    self.cbeta_result = cbetadev.cbeta(**kwargs)

  def set_rota_result(self):
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    sidechain_angles = SidechainAngles(show_errs=True)
    if self.resname == "MSE" : curres = "MET"
    else : curres = self.resname
    all_dict = construct_complete_sidechain(self.three[self.index],self.altloc)
    #print all_dict
    kwargs = self.get_start_kwargs()
    try :
      chis = sidechain_angles.measureChiAngles(
                         self.three[self.index],
                         all_dict.get(self.altloc))
    except AttributeError :
      kwargs['incomplete'] = True
      result = rotamer(**kwargs)
      print('%s is missing some sidechain atoms'%result.id_str(), file=sys.stderr)
      self.rota_result = rotalyze.rotamer(**kwargs)
      return
    #print all_dict.get(self.altloc).keys()
    if None in chis :
      result = rotamer(**kwargs)
      print('%s may be missing some sidechain atoms' %\
                               result.id_str(), file=sys.stderr)
      self.rota_result = rotalyze.rotamer(**kwargs)
      return
    cur_res = self.resname.lower().strip()
    if cur_res == 'mse':
      cur_res = 'met'
    value = self.rota_eval.evaluate(cur_res, chis)
    if value == None :
      self.rota_result = rotalyze.rotamer(**kwargs)
      return
    kwargs['score'] = value * 100
    wrap_chis = self.rotamer_id.wrap_chis(self.resname.strip(),
                          chis,symmetry=False)
    if value >= rotalyze.ALLOWED_THRESHOLD   : evaluation = "Favored"
    elif value >= rotalyze.OUTLIER_THRESHOLD : evaluation = "Allowed"
    else                                     :
      evaluation = "OUTLIER"
      kwargs['outlier'] = True
      kwargs['rotamer_name'] = evaluation
    kwargs['evaluation'] = evaluation
    if evaluation != "OUTLIER" :
      kwargs['outlier'] = False
      kwargs['rotamer_name'] = self.rotamer_id.identify(self.resname, wrap_chis)
      # deal with unclassified rotamers
      if kwargs['rotamer_name'] == '' : kwargs['rotamer_name'] = "UNCLASSIFIED"
    # fill out wrap_chis win None til len == 4
    while (len(wrap_chis) < 4) : wrap_chis.append(None)
    kwargs['chi_angles'] = wrap_chis
    self.rota_result = rotalyze.rotamer(**kwargs)

  def set_rama_result(self,phi_psi_atoms):
    import mmtbx.rotamer
    res_type = self.three.get_ramalyze_key()
    r_name = self.three.get_resnames()[1]
    assert res_type in range(6)
    phi_atoms, psi_atoms = phi_psi_atoms
    res1C,res2N,res2CA,res2C = phi_atoms
    res2N,res2CA,res2C,res3N = psi_atoms
    phi = mmtbx.rotamer.phi_from_atoms(res1C, res2N, res2CA, res2C)
    psi = mmtbx.rotamer.psi_from_atoms(res2N, res2CA, res2C, res3N)
    value = self.rama_eval.evaluate(r_name, [phi, psi])
    ramaType = ramalyze.ramalyze.evalScore(res_type, value)
    is_outlier = ramaType == ramalyze.RAMALYZE_OUTLIER
    c_alphas = None
    if is_outlier :
      c_alphas = []
      for a in self.three.atoms():
        if (a.name.strip() == "CA"):
          a_ = ramalyze.atom(pdb_atom=a)
          c_alphas.append(ramalyze.c_alpha(
            id_str=a_.atom_group_id_str(),
            xyz=a_.xyz))
      assert (len(c_alphas) == 3)
    kwargs = self.get_start_kwargs()
    kwargs.update({'phi':phi,
                   'psi':psi,
                   'rama_type':ramaType,
                   'res_type':res_type,
                   'score':value*100,
                   'outlier':is_outlier,
                   'xyz':ramalyze.get_center(self.three[1]),
                   'c_alphas':c_alphas})
    self.rama_result = ramalyze.ramachandran(**kwargs)

class ComprehensiveResidueValidation(object):

  def __init__(self, pdb_hierarchy):
    self.pdb_hierarchy = pdb_hierarchy
    self.residues = []
    self.validate_residues()

  def validate_residues(self):
    from mmtbx.conformation_dependent_library import generate_protein_threes
    from mmtbx.rotamer import ramachandran_eval,rotamer_eval
    # this is so we generate rama_eval only once
    rama_eval = ramachandran_eval.RamachandranEval()
    rota_eval = rotamer_eval.RotamerEval()
    rotamer_id = rotamer_eval.RotamerID() # loads in the rotamer names
    threes = generate_protein_threes(
        hierarchy = self.pdb_hierarchy,
        include_non_linked=True,
        backbone_only=False,
        geometry=None)
    for i,three in enumerate(threes):
      if i == 0 :
        self.residues.append(ValidationResidue(three,rama_eval,
                                               rota_eval,rotamer_id,index=0))
      self.residues.append(ValidationResidue(three,rama_eval,
                                               rota_eval,rotamer_id))
      if three.end :
        self.residues.append(ValidationResidue(three,rama_eval,
                                               rota_eval,rotamer_id,index=2))

def run(args, out=sys.stdout, quiet=False):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()

  import time
  start_time = time.time()
  pdb_in = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.hierarchy
  hierarchy.atoms().reset_i_seq()
  result = ComprehensiveResidueValidation(pdb_hierarchy = hierarchy)
  print("Elapsed time = %.5f" % (time.time() - start_time))

if (__name__ == "__main__"):
  run(sys.argv[1:])
