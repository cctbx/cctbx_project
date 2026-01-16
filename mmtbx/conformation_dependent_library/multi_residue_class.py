from __future__ import absolute_import, division, print_function
import copy

from mmtbx.conformation_dependent_library.cdl_utils import \
  get_c_ca_n
from mmtbx.conformation_dependent_library.cdl_utils import \
  distance2, get_omega_value, get_phi_psi_angles
from mmtbx.conformation_dependent_library.cdl_utils import \
  get_ca_dihedrals
from mmtbx.conformation_dependent_library.LinkedResidues import LinkedResidues
from six.moves import range

class RestraintsRegistry(dict):
  def __init__(self):
    self.n = {}

  def __repr__(self):
    outl = "RestraintsRegistry"
    outl += "\n  %s(%d)" % (list(self.keys()), len(self))
    outl += "\n  %s" % self.n
    return outl

  def __setitem__(self, key, item):
    if key in self:
      if self[key]!=item:
        self.n.setdefault(key,1)
        self.n[key]+=1
        dict.__setitem__(self, key, (self[key]+item))
    else:
      dict.__setitem__(self, key, item)

class ProteinResidues(LinkedResidues):
  def __init__(self,
               geometry,
               length=3, # CDL & other psi/phi apps
               allow_poly_ca=False,
               registry=None,
               include_non_linked=False,
              ):
    LinkedResidues.__init__(self,
                            geometry,
                            length=length,
                            allow_poly_ca=allow_poly_ca,
                            registry=registry,
                            include_non_linked=include_non_linked,
                            )

  def show(self):
    outl = "%sProteinResidues" % self.length
    for residue in self:
      if residue is not None: outl += " %s(%s)" % (residue.resname, residue.resseq)
      else: outl += ' "%s"' % residue
    outl += " %s" % self.are_linked(return_value=True)
    if self.start is not None: outl += " start=T"
    if self.end is not None: outl += " end=T"
    return outl

  def show_detailed(self):
    outl = "%sProteinResidues" % self.length
    outl += "\nREMARK"
    for residue in self:
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def get_omega_value(self): assert 0

  def _define_omega_a_la_duke_using_limit(self,
                                          omega,
                                          limit=45.,
                                          ):
    if omega is None: return None
    if abs(omega)<limit: return 'cis'
    elif 180-abs(omega)<limit: return 'trans'
    else: return 'twisted'

  def cis_group(self,
                limit=45.,
#                omega_cdl=False, # need last not middle
                verbose=False):
    # is any omega a cis angle?
    # assert not omega_cdl
    #cis_peptide_bond = False
    #omega = self.get_omega_value(omega_cdl=omega_cdl)
    #if omega is None: return None
    omegas = self.get_omega_values()
    assert omegas
    def _is_cis(angle):
      return self._define_omega_a_la_duke_using_limit(angle, limit=limit)=='cis'
    if list(filter(_is_cis, omegas)): return True
    return False
    #if self._define_omega_a_la_duke_using_limit(omega, limit=limit)=='cis':
    #  cis_peptide_bond = True

  def trans_group(self, limit=45.):
    return not self.cis_group(limit=limit)

  def cis_trans_twisted_list(self, limit=45.):
    omegas = self.get_omega_values()
    def _is_cis_trans_twisted(angle):
      return self._define_omega_a_la_duke_using_limit(angle, limit=limit)
    return [_is_cis_trans_twisted(o) for o in omegas]

  def enol_group(self, validate=True):
    from libtbx.utils import Sorry
    assert len(self) in [2,3], 'Enol-peptide only coded for 2,3 peptides'
    for i, residue in enumerate(self):
      if i!=len(self)-2: continue
      if rc:=residue.find_atom_by(name=' HNO'): break
    if rc and validate:
      for i, residue in enumerate(self):
        if i!=len(self)-1: continue
        h_atom=residue.find_atom_by(name=' H  ')
        if h_atom:
          raise Sorry('Enol-peptide should not have a "H" hydrogen on following residue : %s' % h_atom.quote())
    return rc

  def are_linked(self,
                 return_value=False,
                 return_atoms=False,
                 use_distance_always=False,
                 bond_cut_off=2.,
                 allow_poly_ca=False,
                 poly_ca_cut_off=4.,
                 verbose=True):
    '''
    Need to add poly-Calpha chains
      CA-CA 4.5 is use in CaBLAM, maybe shorter
    '''
    allow_poly_ca = allow_poly_ca or self.allow_poly_ca
    d2 = None
    bond_cut_off *= bond_cut_off
    poly_ca_cut_off *= poly_ca_cut_off
    for i, residue in enumerate(self):
      if i==0: continue
      ccn1, outl1 = get_c_ca_n(residue, return_subset=True)
      if self[i-1] is None: # place holder for omega CDL
        return False
      ccn2, outl2 = get_c_ca_n(self[i-1], return_subset=True)
      if ccn1 is None:
        for line in outl1:
          if line not in self.errors:
            self.errors.append(line)
        break
      if ccn2 is None:
        for line in outl2:
          if line not in self.errors:
            self.errors.append(line)
        break
      n = ccn1[2]
      c = ccn2[0]
      if n is None or c is None:
        if not allow_poly_ca: return False
        #poly ca "bonding" is checked only if peptide bond is missing
        #  and if poly ca chains are allowed
        ca1 = ccn1[1]
        ca2 = ccn2[1]
        if ca1 is None or ca2 is None: return False
        d2 = distance2(ca1,ca2)
        if d2<poly_ca_cut_off:
          bond=True
          continue
        else:
          bond=False
          break
      if self.bond_params_table is None:
        d2 = distance2(n,c)
        if d2<bond_cut_off: bond=True
        else: bond=False
      else:
        bond=self.bond_params_table.lookup(c.i_seq, n.i_seq)
        if not bond and use_distance_always:
          # needed for situations where atoms are added and the i_seq is updated
          if distance2(n,c)<bond_cut_off: bond=True
      if not bond:
        break
    else:
      if return_atoms: return c,n
      return True
    if return_value: return d2
    return False

  def get_phi_psi_angles(self): assert 0

  def get_omega_values(self, verbose=False):
    rc=[]
    for i, residue in enumerate(self):
      if i==0: continue
      omega = get_omega_value(residue, self[i-1], verbose=verbose)
      rc.append(omega)
    return rc

  def get_residue_group_from_hierarchy(self, hierarchy, index):
    atom = self[index].atoms()[0]
    for i in range(atom.i_seq, len(hierarchy.atoms())):
      tmp = hierarchy.atoms()[i]
      if tmp.id_str()==atom.id_str(): break
    atom = hierarchy.atoms()[i]
    return atom.parent().parent()

  def trace(self, include_side_chain=False):
    atoms={
      'N':None,
      'CA':None,
      'C':None,
      }
    bonds={
      ('N', 'CA'): [],
      ('CA', 'C'): [],
      # ('C', 'N'): [],
      }
    for atom in self.atoms():
      if atom.name.strip() in atoms:
        atoms[atom.name.strip()]=atom
    # print(atoms)
    for key in bonds:
      n1, n2 = key
      bonds[key]=[atoms[n1], atoms[n2]]
      yield bonds[key]

class TwoProteinResidues(ProteinResidues):
  def get_omega_value(self):
    return get_omega_value(self[1], self[0])

class ThreeProteinResidues(ProteinResidues):
  def get_omega_values(self,
                       #omega_cdl=None,
                       verbose=False,
                       ):
    #assert omega_cdl is None, 'can not use omega_cdl for %sProteinResidues' % self.length
    return ProteinResidues.get_omega_values(self, verbose=verbose)

  def get_phi_psi_atoms(self,
                        only_psi_phi_pairs=True,
                        force_plus_one=False,
                        verbose=False,
                        ):
    if len(self)!=self.length: return None, None
    if force_plus_one: only_psi_phi_pairs=False
    if self[0] is None:
      backbone_i_minus_1 = None
    else:
      backbone_i_minus_1, junk = get_c_ca_n(self[0], return_subset=True)
      assert len(backbone_i_minus_1)==self.length
    backbone_i, junk = get_c_ca_n(self[1], return_subset=True)
    if verbose: print(backbone_i)
    if None in backbone_i: return None
    backbone_i_plus_1, junk = get_c_ca_n(self[2], return_subset=True)
    if verbose: print(backbone_i_plus_1, junk)
    if None in backbone_i_plus_1: return None
    assert len(backbone_i)==self.length
    assert len(backbone_i_plus_1)==self.length
    phi_atoms = [
      backbone_i_minus_1[0],
      backbone_i[2],
      backbone_i[1],
      backbone_i[0],
      ]
    psi_atoms = [
      backbone_i[2],
      backbone_i[1],
      backbone_i[0],
      backbone_i_plus_1[2],
      ]
    atoms = [phi_atoms, psi_atoms]
    if verbose: print(atoms)
    if not only_psi_phi_pairs:
      if self.start:
        psi_atoms = [
          backbone_i_minus_1[2],
          backbone_i_minus_1[1],
          backbone_i_minus_1[0],
          backbone_i[2],
          ]
        atoms.insert(0, psi_atoms)
      if self.end or force_plus_one:
        phi_atoms = [
          backbone_i[0],
          backbone_i_plus_1[2],
          backbone_i_plus_1[1],
          backbone_i_plus_1[0],
          ]
        atoms.append(phi_atoms)
    if verbose:
      for dihedral in atoms:
        print('-'*80)
        for atom in dihedral:
          print(atom.quote())
    return atoms

  def get_phi_psi_angles(self, verbose=False):
    if verbose:
      for residue in self:
        print(residue.id_str())
    return get_phi_psi_angles(self, verbose=verbose)

  def get_ramalyze_key(self,
                       limit=30.,
                       verbose=False,
                       ):
    from mmtbx.validation import ramalyze
    # defined in mmtbx.validation.ramalyze:
    # res_types = ["general", "glycine", "cis-proline", "trans-proline",
    #              "pre-proline", "isoleucine or valine"]
    #
    # This should be consistent with mmtbx/validation/ramalyze.py,
    # lines 179-195. Particularly, prepro comes before ile/val
    if self[1].resname == "PRO":
      if self.cis_group(limit=limit): return ramalyze.RAMA_CISPRO
      else: return ramalyze.RAMA_TRANSPRO
    elif self[1].resname == "GLY": return ramalyze.RAMA_GLYCINE
    elif self[2].resname == "PRO": return ramalyze.RAMA_PREPRO
    elif self[1].resname in ["ILE", "VAL"]: return ramalyze.RAMA_ILE_VAL
    else: return ramalyze.RAMA_GENERAL

  def provide_second_sub_unit_if_unlinked(self):
    # used if residue is appended using superclass method
    if not self.are_linked():
      sub_unit = copy.copy(self) # calls append to delete first sub unit
      while not self.are_linked():
        del self[-1]
      return sub_unit
    return None

  def get_dummy_dihedral_proxies(self, only_psi_phi_pairs=True):
    #
    # Needs testing. One of the candidates is 3j0d, chain I, the first
    # residue is missing CA atom.
    #
    from cctbx.geometry_restraints import dihedral_proxy
    atoms = self.get_phi_psi_atoms(only_psi_phi_pairs=only_psi_phi_pairs)
    proxies = []
    if atoms is None: return proxies
    for dihedral in atoms:
      if None not in dihedral:
        proxy = dihedral_proxy(
            i_seqs=[atom.i_seq for atom in dihedral],
            angle_ideal=0,
            weight=1)
        proxies.append(proxy)
    return proxies

def id_str_for_phi_psi_2(residue, ignore_altloc=False):
  id_str = residue.id_str(suppress_segid=True).replace('pdbres=','').replace('"','')
  if ignore_altloc:
    return ' %s' % id_str
  altloc=''
  for atom in residue.atoms():
    if atom.name not in [' N  ', ' CA ', ' C  ']: continue
    atom_group = atom.parent()
    if atom_group.altloc:
      altloc = atom_group.altloc
      break
  return '%1s%s' % (altloc, id_str)

class FourProteinResidues(ThreeProteinResidues):
  def get_ca_dihedrals(self, verbose=False):
    if verbose:
      for residue in self:
        print(residue.id_str())
    return get_ca_dihedrals(self)

  def id_str_for_phi_psi_2(self, ignore_altloc=False):
    return '%s ~> %s' % ( id_str_for_phi_psi_2(self[1], ignore_altloc),
                          id_str_for_phi_psi_2(self[2], ignore_altloc))

class FiveProteinResidues(FourProteinResidues):
  def get_cablam_info(self):
    assert 0

if __name__=="__main__":
  import sys
  from iotbx import pdb
  from mmtbx.conformation_dependent_library.tst_rdl import get_geometry_restraints_manager
  filename=sys.argv[1]
  pdb_inp = pdb.input(filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(filename)
  pdb_hierarchy.reset_i_seq_if_necessary()
  from mmtbx.conformation_dependent_library import generate_protein_fragments
  for i in range(2,11):
    for threes in generate_protein_fragments(pdb_hierarchy,
                                             geometry_restraints_manager,
                                             length=i,
                                             #verbose=verbose,
                                             ):
      print(threes)
      try: print('  omega   %5.1f' % threes.get_omega_value())
      except: print('  omega is not valid') # intentional
      print('  omegas  %s' % threes.get_omega_values())
      try: print("  cis?    %-5s %s" % (threes.cis_group(), threes.cis_group(limit=30)))
      except: print('  cis? is not valid') # intentional
      try: print("  trans?  %-5s %s" % (threes.trans_group(), threes.trans_group(limit=30)))
      except: print('  trans? is not valid') # intentional
      print(threes.cis_trans_twisted_list())
      # print('  cis/trans/twisted? %s' % ' '.join(threes.cis_trans_twisted_list()))
      try: print("  rama    %s" % threes.get_ramalyze_key())
      except: print('  rama not specified') # intentional
      print('  conf    %s' % threes.is_pure_main_conf())
      try: print('  phi/psi %s' % threes.get_phi_psi_angles())
      except: print('  phi/psi not specified') # intentional
      try: print('  CA dihedrals %s' % threes.get_ca_dihedrals())
      except: print('  CA dihedrals not specified') # intentional
    print("OK",i+2)
