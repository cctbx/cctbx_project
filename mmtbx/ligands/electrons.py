from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
import sys
import time

from functools import cmp_to_key

from libtbx.utils import Sorry

#
# TODO:
#  test on radicals - MTN
#

get_class = iotbx.pdb.common_residue_names_get_class

base_amino_acid_charges = { # all the rest are zero
  'ARG' :  1,
  'ASP' : -1,
  'CYS' :  0, # just to be clear
  'GLU' : -1,
  'HIS' :  1,
  'LYS' :  1,
  }
other_charges = {
  'HOH' : 0,
}
disallowed_element_charges = {
  'N' :  1,
  'O' : -1,
  'C' :  1,
}
default_metal_charges = {
  'Li' : 1,
  'Na' : 1,
  'Mg' : 2,
  'K'  : 1,
  'Ca' : 2,
  'Cu' : 1,
  'Zn' : 2,
}

class atoms(dict):
  def __repr__(self):
    from mmtbx.ligands.chemistry import elements
    outl = ''
    for element in elements:
      outl += '  %2s : %s\n' % (element, self.get(element, None))
    return outl

def get_atom_database():
  from mmtbx.ligands.chemistry import elements, valences, lone_pairs
  from mmtbx.ligands.chemistry import non_metal_indices
  atom_database = atoms()
  for i, element in enumerate(elements):
    atom_database[element.upper()]={'number':i}
  for element, valence in zip(elements, valences):
    if valence>-1:
      atom_database[element.upper()]['valence'] = valence
  for element, lone_pair in zip(elements, lone_pairs):
    if lone_pair>0:
      atom_database[element.upper()]['lone pairs'] = lone_pair
  for i, element in enumerate(elements):
    if i not in non_metal_indices:
      atom_database[element.upper()]['metal']=True
      if element in default_metal_charges:
        # atom_database['valence'] = default_metal_charges[element]*-1
        atom_database[element.upper()]['charge'] = default_metal_charges[element]
      else:
        atom_database[element.upper()]['charge'] = None
  atom_database['D']=atom_database['H']
  return atom_database

class atom_property(dict):
  def __init__(self):
    atom_database = get_atom_database()
    for element, data in atom_database.items():
      self[element] = data

  def __repr__(self):
    outl = 'atom properties\n'
    for element, data in self.items():
      outl += '  %-2s : %s\n' % (element, data)
    return outl

  def get_valence(self, element, effective=True):
    assert effective
    return self.get(element.strip(), {}).get('valence', None)

  def get_lone_pairs(self, element):
    return self.get(element.strip(), {}).get('lone pairs', 0)

  def get_atomic_number(self, element):
    return self.get(element.strip(), {}).get('number', 0)

  def is_metal(self, element):
    return self.get(element.strip(), {}).get('metal', False)

  def get_charge(self):
    return self.get(element.strip(), {}).get('charge', None)

class electron_distribution(dict):
  def __init__(self,
               hierarchy,
               grm,
               specific_atom_charges=None, # a list of selections and charges
               specific_atom_multiplicities=None,
               alternative_location_id=None,
               alternative_location_index=None,
               log=None,
               verbose=False,
               ):
    alternative_location_id='A'
    self.properties = atom_property()
    self.hierarchy = hierarchy
    self.atoms = self.hierarchy.atoms()
    self.grm = grm
    #
    self.xrs = self.hierarchy.extract_xray_structure(
      crystal_symmetry=self.grm.crystal_symmetry)
    self.simple, self.asu = self.grm.get_all_bond_proxies(
      sites_cart=self.xrs.sites_cart())
    #
    self.specific_atom_charges = specific_atom_charges
    self.specific_atom_multiplicities = specific_atom_multiplicities
    self.atoms_with_charges_set = []
    if log is None:
      self.logger = sys.stdout
    else:
      self.logger = log
    self.verbose=verbose
    if [_f for _f in hierarchy.get_conformer_indices().conformer_indices if _f]:
      assert (alternative_location_id is not None or
              alternative_location_index is not None)
    for atom in self.atoms:
      e = self.properties.get_valence(atom.element)
      assert e is not None, ' element %s not found' % atom.element
      metal = self.properties.is_metal(atom.element)
      if metal:
        if atom.element.capitalize() in default_metal_charges:
          self[atom.i_seq]=default_metal_charges[atom.element.capitalize()]*-1
          # self[atom.i_seq]=None
        else:
          print(atom.quote())
          # assert 0, ' charge not found for %s' % atom.quote()
      else:
        self[atom.i_seq] = e
    self.set_charges() # place for selections
    self.validate_metals()
    self.form_bonds()
    self.adjust_for_multiplicity()

  def __repr__(self):
    return self.show()

  def _repr_(self,
             show_all = False,
             show_unpaired = True,
             show_empty_bonds = True,
             ):
    atoms = self.hierarchy.atoms()
    header = 'elec. dist.\n'
    outl = ''
    for key, electrons in self.items():
      if type(key)==type(tuple([])):
        if(show_empty_bonds and electrons==0) or show_all:
          outl += '  %s-%s : %d electrons\n' % (atoms[key[0]].quote(),
                                      atoms[key[1]].quote(),
                                      electrons,
          )
      else:
        assert abs(electrons)<20
        if(show_unpaired and electrons) or show_all:
          outl += '  %s  : %3de\n' % (atoms[key].quote(), electrons)
    if not outl:
      outl = '  molecule neutral'
    outl = '%s%s' % (header, outl)
    return outl

  def show_detailed(self):
    return self._repr_(show_all=True)

  def show(self):
    return self._repr_()

  def show_bonds(self):
    assert 0

  def show_residues(self):
    outl = 'elec. dist.\n'
    for residue_group in self.hierarchy.residue_groups():
      sum_e = 0
      for atom in residue_group.atoms():
        sum_e += self[atom.i_seq]
      outl += '  %s : %2d\n' % (atom.parent().id_str(), sum_e)
    return outl

  def show_dot(self, filename='molecule.png'):
    import graphviz
    f = graphviz.Digraph(filename=filename)
    indx=[]
    element=[]
    edges=[]
    atoms=self.hierarchy.atoms()
    for key, electrons in self.items():
      if type(key)==type(tuple([])):
        if atoms[key[0]].element in ['H', 'D']: continue
        if atoms[key[1]].element in ['H', 'D']: continue
        for e in range(electrons):
          edges.append([str(key[0]), str(key[1])])
      else:
        if atoms[key].element in ['H', 'D']: continue
        indx.append(str(key))
        element.append('%s' % (atoms[key].name))
        if electrons: element[-1]+=' (%s)' % electrons
    for name, position in zip(indx, element):
      f.node(name, position)
    for e1, e2 in edges:
      f.edge(e1,e2)
    print(f.source)
    f.render(view=True)

  def _generate_atoms(self):
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      yield key

  def _generate_bonds(self):
    for key, electrons in self.items():
      if type(key)==type(tuple([])):
        yield key

  def __setitem__(self, i_seq, electrons):
    if electrons<-1:
      if self.properties.get_lone_pairs(self.atoms[i_seq].element):
        electrons+=2
    dict.__setitem__(self, i_seq, electrons)

  def _add_electron_to_bond(self, i_seqs, verbose=False):
    if verbose:
      atoms = self.hierarchy.atoms()
      print('_add_electron_to_bond')
      print(i_seqs, atoms[i_seqs[0]].quote(), atoms[i_seqs[1]].quote())
      print(self)
    if i_seqs not in self:
      tmp = (i_seqs[1], i_seqs[0])
      i_seqs=tmp
    self[i_seqs]+=1
    self[i_seqs[0]]-=1
    self[i_seqs[1]]-=1

  def _subtract_electron_from_bond(self, i_seqs, verbose=False):
    if verbose:
      atoms = self.hierarchy.atoms()
      print('_subtract_electron_from_bond')
      print(i_seqs, atoms[i_seqs[0]].quote(), atoms[i_seqs[1]].quote())
      print(self)
    if i_seqs not in self:
      tmp = (i_seqs[1], i_seqs[0])
      i_seqs=tmp
    self[i_seqs]-=1
    self[i_seqs[0]]+=1
    self[i_seqs[1]]+=1

  def set_charges(self):
    atoms = self.hierarchy.atoms()
    for key, electrons in self.items():
      element = atoms[key].element.strip()
      if element.capitalize() in default_metal_charges:
        self[key]=default_metal_charges[element.capitalize()]*-1

    for atom in atoms:
      element = atom.element.strip()
      if self.properties.is_metal(element):
        if atom.charge:
          self[atom.i_seq]=atom.charge_as_int()*-1

    if self.specific_atom_charges:
      for i, sac in enumerate(self.specific_atom_charges):
        metal_asc = self.hierarchy.atom_selection_cache()
        metal_sel = metal_asc.selection(sac.atom_selection)
        metal_hierarchy = self.hierarchy.select(metal_sel)
        for i, atom in enumerate(metal_hierarchy.atoms()):
          self[atom.i_seq]=sac.charge*-1
          self.atoms_with_charges_set.append(atom.i_seq)
          assert i<1

  def adjust_for_multiplicity(self):
    if self.specific_atom_multiplicities:
      for i, mac in enumerate(self.specific_atom_multiplicities):
        radical_asc = self.hierarchy.atom_selection_cache()
        radical_sel = radical_asc.selection(mac.atom_selection)
        radical_hierarchy = self.hierarchy.select(radical_sel)
        for i, atom in enumerate(radical_hierarchy.atoms()):
          if self[atom.i_seq] and not atom.i_seq in self.atoms_with_charges_set:
            if mac.multiplicity==2:
              self[atom.i_seq]=0
              print('\nCharge on %s changed to zero because of multiplicity. CHECK!' % atom.quote(),
                    file=self.logger)
              break

  def _has_metal(self, atom1, atom2):
    is_metal_count = [0,1][self.properties.is_metal(atom1.element)]
    is_metal_count+= [0,1][self.properties.is_metal(atom2.element)]
    if is_metal_count==2:
      print('\nMore than one metal in a bond can lead to issues. CHECK!', file=self.logger)
    return is_metal_count

  def is_metal_bond(self, key):
    assert key in self
    atoms = self.hierarchy.atoms()
    return self._has_metal(atoms[key[0]], atoms[key[1]])

  def validate_metals(self):
    for atom in self.hierarchy.atoms():
      element = atom.element.strip()
      if self.properties.is_metal(element):
        if atom.i_seq not in self:
          raise Sorry('''Charge error:
  Atom %s does not have a charge specified. Use PHIL parameter
  specific_atom_charges or specify in the input model.
  ''' % (
    atom.quote(),
    )
  )

  def _is_max_bond_valence(self, i_seq):
    max_valence = self.properties.get_valence(self.atoms[i_seq].element)
    lp = self.properties.get_lone_pairs(self.atoms[i_seq].element)
    if lp: max_valence += lp*2
    for bond in self._generate_bonds():
      if i_seq in bond:
        max_valence -= self[bond]
        if max_valence==0: break
    return max_valence==0

  def _is_c_c_bond(self, i_seq, j_seq, tetra=['C']):
    # not checked for bonding
    if (self.atoms[i_seq].element.strip() in tetra and
        self.atoms[j_seq].element.strip() in tetra
      ): return True
    return False

  def _can_denote_electron_to_covalent_bond(self,
                                            i_seq,
                                            j_seq,
                                            dangling=False,
                                            verbose=False):
    if verbose:
      print('processing %s %s' % (self.atoms[i_seq].quote(), self.atoms[j_seq].quote()))
    if self[i_seq]>0 and self[j_seq]>0:
      if verbose:
        print('bonding %s %s' % (self.atoms[i_seq].quote(), self.atoms[j_seq].quote()))
      return True
    elif self[i_seq]==0 and self[j_seq]==0:
      return False
    atom1 = self.atoms[i_seq]
    if atom1.element_is_hydrogen() and self[i_seq]==0: return False
    atom2 = self.atoms[j_seq]
    if atom2.element_is_hydrogen() and self[j_seq]==0: return False
    if self._is_max_bond_valence(i_seq) or self._is_max_bond_valence(j_seq):
      return False
    assert i_seq==atom1.i_seq
    assert j_seq==atom2.i_seq
    if atom1.element_is_hydrogen():
      hydrogen = atom1
      other = atom2
    elif atom2.element_is_hydrogen():
      hydrogen = atom2
      other = atom1
    elif dangling and self.properties.get_lone_pairs(atom1.element)==0:
      if verbose: print('atom has no lone pairs',atom1.quote())
      return False
    elif dangling and self.properties.get_lone_pairs(atom2.element)==0:
      if verbose: print('atom has no lone pairs',atom2.quote())
      return False
    else:
      an1 = self.properties.get_atomic_number(atom1.element)
      an2 = self.properties.get_atomic_number(atom2.element)
      if an2>an1:
        dummy = atom1
        atom1 = atom2
        atom2 = dummy
      if self.properties.get_lone_pairs(atom1.element):
        lone_pair = atom1
        other = atom2
      elif self.properties.get_lone_pairs(atom2.element):
        lone_pair = atom2
        other = atom1
      else:
        return False
      if verbose:
        print('other-lp   %s-%s' % (other.quote(), lone_pair.quote()))
        print(self.properties.get_lone_pairs(atom1.element))
        print(self.properties.get_lone_pairs(atom2.element))
        print(self[lone_pair.i_seq], self[other.i_seq])
      if self[other.i_seq]>0:
        return True
      return False
    if self.properties.get_lone_pairs(other.element):
      #self[other.i_seq]+=2
      if verbose: print('hydrogen-X lone pair TRUE')
      return True
    return None

  def form_bonds_using_simple(self, extend_based_on_proximity=False, verbose=False):
    if self.verbose or verbose: verbose=1
    # verbose=1
    def _get_sum_lone_pairs(bp):
      i_seq, j_seq = bp.i_seqs
      lp1 = self.properties.get_lone_pairs(atoms[i_seq].element)
      lp2 = self.properties.get_lone_pairs(atoms[j_seq].element)
      return lp1+lp2
    def _sort_lone_pairs(bp1, bp2):
      slp1 = _get_sum_lone_pairs(bp1)
      slp2 = _get_sum_lone_pairs(bp2)
      if slp2>slp1: return -1
      return 1
    def generate_bonds_from_simple(simple, sort_on_lone_pairs=False):
      if sort_on_lone_pairs:
        l = []
        for bp in simple:
          l.append(bp)
        l.sort(key=cmp_to_key(_sort_lone_pairs))
        l.reverse()
        later = []
        for bp in l:
          atom1 = atoms[bp.i_seqs[0]]
          atom2 = atoms[bp.i_seqs[1]]
          if atom1.parent().parent().resseq != atom2.parent().parent().resseq:
            later.append(bp)
            continue
          yield bp
        for bp in later:
          yield bp
      else:
        assert 0
    def generate_bonds(simple, asu, sort_on_lone_pairs=False):
      for rc in generate_bonds_from_simple(simple,
                                           sort_on_lone_pairs=sort_on_lone_pairs,
                                           ):
        yield rc
    def generate_atoms_from_simple(simple):
      for bp in simple:
        # assert bp.origin_id in [0,3], ' origin_id "%s"' % bp.origin_id
        i_seq, j_seq = bp.i_seqs
        assert i_seq in self
        assert j_seq in self
        atom1 = atoms[i_seq]
        atom2 = atoms[j_seq]
        yield bp, i_seq, j_seq, atom1, atom2
    def generate_atoms(simple, asu):
      for rc in generate_atoms_from_simple(simple):
        yield rc
      assert not len(asu)
    ###
    xrs = self.hierarchy.extract_xray_structure(
      crystal_symmetry=self.grm.crystal_symmetry)
    simple, asu = self.grm.get_all_bond_proxies(sites_cart=xrs.sites_cart())
    # need to check asu...
    # need to filter out H-bonds
    # look for metal coordination
    # maybe use origin_id
    metal_coordination = []
    for bp, i_seq, j_seq, atom1, atom2 in generate_atoms(simple, asu):
      assert bp.i_seqs not in self
      if is_metal(atom1, atom2):
        self[bp.i_seqs]=0
        metal_coordination.append(i_seq)
        metal_coordination.append(j_seq)
        self[bp.i_seqs]+=1
        if verbose: print('metal',self)
    # look for single (non-metal) bonds
    for bp, i_seq, j_seq, atom1, atom2 in generate_atoms(simple, asu):
      if is_metal(atom1, atom2): continue
      mc = None
      if i_seq in metal_coordination:
        mc = atoms[i_seq]
        other = atoms[j_seq]
      elif j_seq in metal_coordination:
        mc = atoms[j_seq]
        other = atoms[i_seq]
      if mc:
        if other.element_is_hydrogen():
          continue
      self[bp.i_seqs]=0
      if _can_denote_electron_to_covalent_bond(i_seq, j_seq):
        self._add_electron_to_bond(bp.i_seqs)
        if verbose: print('single: %s-%s\n%s' % (atoms[i_seq].quote(),
                                                 atoms[j_seq].quote(),
                                                 self))
    # look for double bonds
    for bp in generate_bonds_from_simple(simple,
                                         sort_on_lone_pairs=True,
                                         ):
      if bp.i_seqs not in self: continue
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      while self[i_seq]>0 and self[j_seq]>0:
        self._add_electron_to_bond(bp.i_seqs)
        if verbose: print('double',self)
        if verbose: print('bonding 2',atoms[i_seq].quote(), atoms[j_seq].quote())
    hypers = []
    for bp in simple:
      if bp.i_seqs not in self: continue
      if verbose: print('hyper',self)
      i_seq, j_seq = bp.i_seqs
      assert i_seq in self
      assert j_seq in self
      while _can_denote_electron_to_covalent_bond(i_seq,
                                                  j_seq,
                                                  verbose=verbose):
        self._add_electron_to_bond(bp.i_seqs)
    # remove HG on sulfur bridge
    # self.check_sulfur_bridge()

  def generate_bond_i_seqs(self, verbose=False):
    for bpc in [self.simple, self.asu]:
      for bp in bpc:
        if hasattr(bp, 'j_seq'):
          i_seqs = [bp.i_seq, bp.j_seq]
          i_seqs.sort()
          i_seqs=tuple(i_seqs)
          # i_seq, j_seq = i_seqs
        else:
          i_seqs = bp.i_seqs
          # i_seq, j_seq = bp.i_seqs
        if verbose:
          print(i_seqs,
                self.atoms[i_seqs[0]].quote(),
                self.atoms[i_seqs[1]].quote(),
                )
        yield i_seqs

  def get_bonds_containing_i_seq(self, i_seq):
    rc = []
    for i_seqs in self.generate_bond_i_seqs():
      if i_seq in i_seqs:
        rc.append(i_seqs)
    return rc

  def get_cycle_charge_count(self, cycle):
    tmp=[]
    for c in cycle:
      tmp.append(c[0])
      tmp.append(c[1])
    rc=0
    for key, item in self.items():
      if key in tmp:
        if item: rc+=1
    return rc

  def get_cycle_charge(self, cycle):
    tmp=[]
    for c in cycle:
      tmp.append(c[0])
      tmp.append(c[1])
    rc=0
    for key, item in self.items():
      if key in tmp:
        if item: rc-=self[key]
    return rc

  def process_dangling_heavy_atoms(self, verbose=False):
    for key, electrons in self.items():
      bonds = self.get_bonds_containing_i_seq(key)
      if len(bonds)==1:
        i_seq, j_seq = bonds[0]
        if self._can_denote_electron_to_covalent_bond(i_seq,
                                                      j_seq,
                                                      dangling=True,
                                                      verbose=verbose):
          self._add_electron_to_bond((i_seq, j_seq))
          if verbose: print('dangling: %s-%s\n' % (self.atoms[i_seq].quote(),
                                                   self.atoms[j_seq].quote(),
                                                  ))

  def does_using_hyper_remove_electrons(self, cycle):
    def _generate_ij(cycle):
      for i_seq, j_seq in cycle:
        yield i_seq, j_seq
        yield j_seq, i_seq
    if self.get_cycle_charge(cycle)!=-1: return
    for i_seq, j_seq in _generate_ij(cycle):
      if self[i_seq]==1:
        bonds = self.get_bonds_containing_i_seq(i_seq)
        for b_i_seq, b_j_seq in bonds:
          rc = self._can_denote_electron_to_covalent_bond(b_i_seq, b_j_seq)
          if rc:
            self._add_electron_to_bond((i_seq, j_seq))
          if self.get_cycle_charge(cycle)==1:
            break
      if self.get_cycle_charge(cycle)==1:
        break

  def form_bonds_using_networkx(self, verbose=False):
    import networkx as nx
    g = nx.DiGraph()
    #
    def generate_atom_nodes():
      for atom in self.atoms:
        yield (atom.i_seq,
               {'element':atom.element.strip(),
                'i_seq': atom.i_seq,
               })
    def generate_bond_edges(extend_based_on_proximity=False, verbose=False):
      for i_seqs in self.generate_bond_i_seqs():
        if i_seqs in self: continue
        self[i_seqs]=0
        i_seq, j_seq = i_seqs
        if self._can_denote_electron_to_covalent_bond(i_seq, j_seq):
          self._add_electron_to_bond(i_seqs)
          if verbose: print('single: %s-%s\n%s' % (self.atoms[i_seq].quote(),
                                                   self.atoms[j_seq].quote(),
                                                   self))
        yield i_seqs
    #
    t0=time.time()
    g.add_nodes_from(generate_atom_nodes())
    g.add_edges_from(generate_bond_edges(verbose=verbose))
    h = g.to_undirected()
    if verbose: print('  Created graphs of molecule : %0.1fs' % (time.time()-t0))
    self.process_dangling_heavy_atoms()
    cycle_bases = nx.cycle_basis(h)
    done_cycles = []
    t0=time.time()
    for i_seq, (node, attrs) in enumerate(g.nodes(data=True)):
      if attrs['element'] in ['H', 'D']: continue
      assert i_seq==node, '%s %s' % (i_seq, node)
      # =O
      if len(h.adj[i_seq])==1:
        j_seq=list(h.adj[i_seq].keys())[0]
        if self._can_denote_electron_to_covalent_bond(i_seq, j_seq):
          self._add_electron_to_bond((i_seq, j_seq))
          if verbose: print('double: %s-%s\n' % (self.atoms[i_seq].quote(),
                                                 self.atoms[j_seq].quote(),
                                                ))
      # rings
      cycle=[]
      for cb in cycle_bases:
        if i_seq in cb:
          for e in g.edges:
            if e[0]in cb and e[1] in cb:
              cycle.append(e)

      # try:
      #   cycle = nx.find_cycle(g, i_seq, orientation='ignore')
      # except Exception:
      #   pass
      if not cycle: continue
      tmp = []
      for bond in list(cycle): tmp.append(bond[0])
      tmp.sort()
      if tmp in done_cycles: continue
      done_cycles.append(tmp)
      tries=10
      cycle_charge_count=self.get_cycle_charge_count(cycle)
      subtract=[]
      while cycle_charge_count and tries:
        tries-=1
        if not tries:
          self.does_using_hyper_remove_electrons(cycle)
        while subtract:
          i_seqs = subtract.pop()
          self._subtract_electron_from_bond(i_seqs)
        for filter_non_tetra_coordinate in range(2,-1,-1):
          import random
          # cycle=_sort_on_element(cycle, self.atoms)
          for i_seq, j_seq in cycle:
            if filter_non_tetra_coordinate:
              if not self._is_c_c_bond(i_seq, j_seq):
                if verbose: print('skipping C-C bond')
                continue
            if self[i_seq]>0 and self[j_seq]>0:
              self._add_electron_to_bond((i_seq, j_seq))
              subtract.append((i_seq, j_seq))
              if verbose: print('rings : %s-%s\n%s' % (self.atoms[i_seq].quote(),
                                                       self.atoms[j_seq].quote(),
                                                       self))
        cycle_charge_count=self.get_cycle_charge_count(cycle)
        random.shuffle(cycle)
    if verbose: print('  Double & rings : %0.1fs' % (time.time()-t0))
    #
    # hyper and triple
    #
    t0=time.time()
    for i_seqs in self.generate_bond_i_seqs():
    # for bp in self.simple:
      if i_seqs not in self: continue
      if verbose: print('hyper',self)
      i_seq, j_seq = i_seqs
      assert i_seq in self
      assert j_seq in self
      while self._can_denote_electron_to_covalent_bond(i_seq,
                                                       j_seq,
                                                       verbose=verbose):
        self._add_electron_to_bond(i_seqs)
        if verbose: print('hyper : %s-%s\n' % (self.atoms[i_seq].quote(),
                                               self.atoms[j_seq].quote(),
                                              ))
    if verbose: print('  Hyper & triple : %0.1fs' % (time.time()-t0))

  def form_bonds(self, extend_based_on_proximity=False, verbose=False):
    if self.verbose or verbose: verbose=1
    # verbose=1
    self.form_bonds_using_networkx(verbose=verbose)

  def check_sulfur_bridge(self, verbose=False):
    assert 0
    atoms = self.hierarchy.atoms()
    for i_seq in self._generate_atoms():
      for j_seq in self._generate_atoms():
        if j_seq==i_seq: break
        atom1 = atoms[i_seq]
        atom2 = atoms[j_seq]
        if self[i_seq]<0 and self[j_seq]<0:
          bond0 = bond1 = bond2 = None
          for key in self:
            if type(key)==type(tuple([])):
              if i_seq in key and j_seq in key:
                bond0 = key
              elif i_seq in key and not bond1:
                other1=list(key)
                other1.remove(i_seq)
                if atoms[other1[0]].element_is_hydrogen():
                  bond1 = key
              elif j_seq in key and not bond2:
                other2=list(key)
                other2.remove(j_seq)
                if atoms[other2[0]].element_is_hydrogen():
                  bond2 = key
          if bond0 and bond1 and bond2:
            if verbose:
              print('-'*80)
              print(bond0, bond1, bond2)
              print(atoms[bond0[0]].quote())
              print(atoms[bond0[1]].quote())
              print(atoms[bond1[0]].quote())
              print(atoms[bond1[1]].quote())
              print(atoms[bond2[0]].quote())
              print(atoms[bond2[1]].quote())
            self[bond1]-=1
            self[bond2]-=1
            self[bond1[0]]+=1
            self[bond1[1]]+=1
            self[bond2[0]]+=1
            self[bond2[1]]+=1
            assert 0

  def extend_based_on_proximity(self):
    # use available electrons and proximity
    # needs more care or does not need bond proxies
    rc = self.get_possible_covalent_bonds()
    atoms = self.hierarchy.atoms()
    for i_seq, j_seq in rc:
      if verbose:
        print('  forming bond between %s %s' % (atoms[i_seq].quote(),
                                                atoms[j_seq].quote()))
      assert (i_seq, j_seq) not in self
      self[(i_seq, j_seq)] = 1
      self[i_seq]-=1
      self[j_seq]-=1

  def get_possible_covalent_bonds(self):
    def distance2(xyz1, xyz2):
      sum = 0
      for i in range(3): sum+=(xyz2[i]-xyz1[i])**2
      return sum
    rc = []
    atoms = self.hierarchy.atoms()
    for i_seq in self._generate_atoms():
      if self[i_seq]<1: continue
      for j_seq in self._generate_atoms():
        if j_seq==i_seq: break
        if self[j_seq]<1: continue
        atom1 = atoms[i_seq]
        atom2 = atoms[j_seq]
        # exclude H-H
        if atom1.element_is_hydrogen() and atom2.element_is_hydrogen(): continue
        # terminal atoms on a single amino acid C..N
        if not (atom1.element_is_hydrogen() or atom2.element_is_hydrogen()):
          continue
        d2 = distance2(atoms[i_seq].xyz, atoms[j_seq].xyz)
        if atom1.element_is_hydrogen() or atom2.element_is_hydrogen():
          if d2<1.5:
            rc.append([i_seq, j_seq])
          continue
        assert d2>9, ' %s-%s is %0.1f' % (atoms[i_seq].quote(),
                                          atoms[j_seq].quote(),
                                          d2,
                                          )
    return rc

  def validate_atomic_formal_charges(self, verbose=False):
    data = {'*'   : {'N'  : [-1,0,1],
                     'OXT': [1,0],
                     },
            'LYS' : {'NZ' : [-1]},
            'GLU' : {'OE2': [1]},
            'ASP' : {'OD2': [1]},
            }
    rc = []
    for i_seq in self._generate_atoms():
      atom = self.atoms[i_seq]
      residue_data = data.get(atom.parent().resname, {})
      residue_data.update(data['*'])
      if not residue_data:
        if self[i_seq]:
          rc.append(i_seq)
        assert 0
      else:
        if self[i_seq] in residue_data.get(atom.name.strip(), [0]): continue
      residue_data = data['*'] # needs to be only AA?
      if self[i_seq] in residue_data.get(atom.name.strip(), [0]): continue
      rc.append(i_seq)
    if verbose:
      for i_seq in rc:
        print(i_seq, self.atoms[i_seq].quote())
    return rc

  def get_total_charge(self):
    total=0
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      total+=electrons
    return total*-1

  def get_charged_atoms(self):
    rc = []
    atoms = self.hierarchy.atoms()
    for key, electrons in self.items():
      if type(key)==type(tuple([])): continue
      if electrons:
        rc.append([ atoms[key],electrons])
    return rc

  def validate(self, ignore_water=False, raise_if_error=True):
    charged_atoms = self.get_charged_atoms()
    charged_residues = {}
    rc = {}

    atoms = self.hierarchy.atoms()
    for key, electrons in self.items():
      if type(key)==type(tuple([])):
        if self.is_metal_bond(key): pass
        elif electrons==0:
          outl = 'No electrons allocated to bond: %s-%s' % (
            atoms[key[0]].quote(),
            atoms[key[1]].quote(),
          )
          if raise_if_error: raise Sorry(outl)
          rc.setdefault(outl, [])
          rc[outl].append([ atoms[key[0]].quote(),
                            atoms[key[1]].quote(),
                            key])
      else:
        assert abs(electrons)<10
        disallowed = disallowed_element_charges.get(atoms[key].element, None)
        outl = 'Element has strange number of electrons  %s  : %d' % (
          atoms[key].element,
          electrons)
        if electrons!=0 and disallowed is not None:
          def _comp_disallowed(actual, disallowed):
            if disallowed<0: return actual<=disallowed
            elif disallowed>0: return actual>=disallowed
            assert 0
          if _comp_disallowed(electrons, disallowed):
            if raise_if_error: raise Sorry(outl)
            rc.setdefault(outl, [])
            rc[outl].append([atoms[key].quote(), key])

    terminals = {}
    for atom, charge in charged_atoms:
      if atom.name in [' OXT']: terminals[atom.parent().id_str()]=charge
      ag = atom.parent()
      if get_class(ag.resname) in ['common_amino_acid']:
        base = base_amino_acid_charges.get(ag.resname, 0)
        tmp = charged_residues.setdefault(ag.id_str(), base)
        tmp += charge
        charged_residues[ag.id_str()] = tmp

      if ag.resname in other_charges:
        if ignore_water and ag.resname in ['HOH']: continue
        if charge!=other_charges[ag.resname]:
          outl = '  Residue %s has a problem with the charge : %s!=%s' % (
            ag.resname,
            charge,
            other_charges[ag.resname]
            )
        if raise_if_error: raise Sorry(outl)
        rc.setdefault(outl, [])
        rc[outl].append(atom.quote())

    for ag in self.hierarchy.atom_groups():
      delta = 1
      if ag.resname in ['HIS']: delta=2
      terminal_adjust = ag.id_str() in terminals
      charge = charged_residues.get(ag.id_str(), 0)
      outl = 'Unlikely charge for %s of %s' % (ag.resname, charge)
      if abs(charge-base_amino_acid_charges.get(ag.resname, 0)-int(terminal_adjust)) > delta:
        if raise_if_error: raise Sorry(outl)
        rc.setdefault(outl, [])
        rc[outl].append('"%s"' % ag.id_str())
    return rc

  def report(self, ignore_water=False, show_detailed=False):
    answers = {
      'Residue HOH has a problem with the charge : 2!=0' : \
        'Hydrogen atoms not added to water',
      'Element has strange number of electrons  N  : 1' : \
        'N terminal (or break) missing hydrogen atoms',
      'Element has strange number of electrons  O  : -1' : \
        'C terminal (or break) missing oxygen atoms',
    }
    report = self.validate(ignore_water=ignore_water,
                           raise_if_error=False)
    outl=''
    for key, item in sorted(report.items()):
      outl += '\n  %s\n' % key.strip()
      for instance in item:
        i=instance
        if type(instance)==type([]):
          i=instance[0]
        outl += '    %s\n' % i
      if show_detailed:
        answer = answers.get(key.strip(), None)
        if answer:
          outl += '\n     HINT: %s\n' % answer
        else:
          if key.find('Unlikely charge for')>-1 and int(key.split()[-1])>1:
            outl += '\n     HINT: %s\n' % 'Missing side chain atoms'
          elif key.find('No electrons allocated to bond:')>-1:
            outl += '\n     HINT: %s\n' % 'Too many hydrogen atoms'
          else:
            pass
    if outl:
      outl = 'Validation report\n%s' % outl
      print(outl)
    return report

from libtbx.program_template import ProgramTemplate
from libtbx.utils import null_out
from libtbx import group_args

master_phil_str = '''
input
{
  selection = None
    .type = atom_selection
  ignore_water = False
    .type = bool
}
output
  .style = menu_item auto_align
{
  file_name_prefix = None
    .type = path
    .short_caption = Prefix for file name
    .help = Prefix for file name
    .input_size = 400
}
'''

class Program(ProgramTemplate):
  description = '''
Count electrons

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
'''
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    model = self.data_manager.get_model()
    model.set_log(null_out())
    model.process(make_restraints=True)
    if self.params.input.selection:
      new_model = model.selection(self.params.input.selection)
      new_model = model.select(new_model)
      model = new_model
    t0=time.time()
    self.atom_valences = electron_distribution(
      model.get_hierarchy(), # needs to be altloc free
      model.get_restraints_manager().geometry,
      verbose=False,
    )
    print('Distribution time : %01.fs' % (time.time()-t0))
    print('='*80)
    print(self.atom_valences)
    self.report = self.atom_valences.report(
      ignore_water=self.params.input.ignore_water,
      show_detailed=True,
      )
    self.total_charge = self.atom_valences.get_total_charge()

  def get_results(self):
    return group_args(atom_valences = self.atom_valences,
                      validation = self.report,
                      total_charge = self.total_charge,
                      )

def run(pdb_filename=None,
        raw_records=None,
        return_formal_charges=False,
        verbose=False,
        cif_objects=None,
        ):
  # legacy from Q|R...
  if pdb_filename:
    # Read file into pdb_input class
    inp = iotbx.pdb.input(file_name=pdb_filename)
  elif raw_records:
    inp = iotbx.pdb.input(lines=raw_records, source_info='lines from PDB')
  else:
    assert 0

  # create a model manager
  from io import StringIO
  log = StringIO()
  default_scope = mmtbx.model.manager.get_default_pdb_interpretation_scope()
  working_params = default_scope.extract()
  # optional???
  working_params.pdb_interpretation.automatic_linking.link_metals=True
  model = mmtbx.model.manager(
    model_input = inp,
    restraint_objects=cif_objects,
    log = log,
  )
  model.process(make_restraints=True,
    pdb_interpretation_params = working_params)
  # get xray structure
  xrs = model.get_xray_structure()
  grm = model.get_restraints_manager()
  t0=time.time()
  atom_valences = electron_distribution(
    model.get_hierarchy(), # needs to be altloc free
    model.get_restraints_manager().geometry,
    verbose=verbose,
  )
  if verbose: print(atom_valences)
  total_charge = atom_valences.get_total_charge()
  rc = atom_valences.validate_atomic_formal_charges()
  if return_formal_charges: return atom_valences
  return total_charge
