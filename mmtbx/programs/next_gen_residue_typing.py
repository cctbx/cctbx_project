"""Report on Ramachandran values for a model"""
from __future__ import absolute_import, division, print_function

import os
import iotbx.phil
from libtbx.program_template import ProgramTemplate
# from six.moves import range
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
from libtbx.utils import Sorry

from cctbx.geometry_restraints.linking_class import linking_class
origin_ids = linking_class()

master_phil_str = """
display_peptide_types = True
  .type = bool
display_rmsz_outliers = True
  .type = bool
rmsz_limit = 3.
  .type = float
legend = True
  .type = bool
quiet = False
  .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_phil_str)

def print_list(ll, l=3, sep=', '):
  tmp=[]
  # for t in ll: tmp.append(t.strip())
  tmp=ll
  if len(ll)>l:
    return '%s ...' % sep.join(tmp[:l])
  return '%s' % sep.join(tmp)

class lookup_list(list):
  # def __init__(self):

  def __repr__(self):
    outl='MATTS key'
    attrs = [ 'resn.',
              'start',
              'end',
              'side',
      ]
    for i, attr in enumerate(attrs):
      outl += ', %s : %s' % (attr, self[i])
      if i+1>=len(self): break
    return outl

  def as_tuple(self):
    return tuple(self)

  def append1(self, item):
    print(item)
    if len(self)==1 and item is None: assert 0
    list.append(self, item)

def get_atom_names(rg):
  atoms=rg.atoms()
  tmp=set()
  for atom in atoms:
    tmp.add(atom.name.strip())
  return tmp

def check_atom_names(atom_names, loop, verbose=False):
  rc=[]
  if verbose: print(atom_names)
  for names in loop:
    if set(names).intersection(atom_names) == set(names):
      rc.append(names)
      break
  else:
    rc=None
  return rc

def Boedeker(rc):
  info=[]
  resn, start, end, side = rc
  if start:
    if start[0]==['H1', 'H2', 'H3']:
      info.append('NH3+ termini')
    elif start[0]==['H2', 'H3']:
      if resn in ['PRO']:
        info.append('NH2+ PRO termini')
    elif start[0]:
      info.append(start)
  if end:
    if end[0]==['OXT', 'HXT']:
      info.append('COOH termini')
    elif end[0]==['OXT']:
      info.append('COO- termini')
    elif end[0]:
      info.append(end)
  if side:
    if side[0]==['HD1', 'HE2']:
      assert resn in ['HIS']
      info.append('Doubly protonated %s' % resn)
    elif side[0] in[['HD1'], ['HE2']]:
      assert resn in ['HIS']
      info.append('Singly protonated %s (%s)' % (resn, side[0][0]))
    elif side[0]==['HG']:
      assert resn in ['CYS']
      info.append('Protonated %s' % (resn))
    else:
      print(side)
      assert 0
  return info

def process_residue(one, dna_rna_residues=False, index=0, verbose=False):
  #
  # resname
  #
  key=lookup_list()
  key.append(one.get_resnames()[index].strip())
  rg=one[index]
  #
  # start of chain
  #
  if one.start:
    if dna_rna_residues:
      loop=[["HO3'"]]
    else:
      loop=[['H1', 'H2', 'H3'],
            ['H1', 'H2'],
            ['H2', 'H3'],
            # ['H'],
           ]
    atom_names=get_atom_names(rg)
    adding = check_atom_names(atom_names, loop, verbose=verbose)
    if adding is None:
      if dna_rna_residues:
        adding = "Unknown 3'"
      else:
        if 'H' in atom_names:
          adding = 'One H N terminal'
        else:
          adding = 'Broken N terminal'
    key.append(adding)
  else:
    if verbose: print('start else')
    key.append(None)
  #
  # end of chain
  #
  if one.end:
    atom_names=get_atom_names(rg)
    if dna_rna_residues:
      loop=[]
    else:
      loop=[['OXT', 'HXT'],
            ['OXT'],
            ]
    atom_names=get_atom_names(rg)
    adding = check_atom_names(atom_names, loop, verbose=verbose)
    if adding is None:
      if dna_rna_residues:
        adding = "Unknown 5'"
      else:
        adding = 'Broken C terminal'
    key.append(adding)
  else:
    if verbose: print('end else')
    key.append(None)
  #
  # middle
  #
  atom_names=get_atom_names(rg)
  if key[0]=='HIS':
    loop=[['HD1', 'HE2'],
          ['HD1'],
          ['HE2']]
    adding = check_atom_names(atom_names, loop, verbose=verbose)
  elif key[0]=='CYS':
    loop=[['HG']]
    adding = check_atom_names(atom_names, loop, verbose=verbose)
  elif key[0]=='ASP':
    loop=[['HD2']]
    adding = check_atom_names(atom_names, loop, verbose=verbose)
  elif key[0]=='GLU':
    loop=[['HE2']]
    adding = check_atom_names(atom_names, loop, verbose=verbose)
  # elif key[0]=='LYS':
  #   loop=[['HZ1']]
  else:
    # key.append(None)
    adding=None
    if verbose: print('middle else')
  key.append(adding)
  return key

def compute(hierarchy, grm, verbose=False):
  if verbose: hierarchy.show()
  def _process_rg_for_compute(one, two, index):
    rc = process_residue(one, dna_rna_residues=dna_rna_residues)
    interesting = Boedeker(rc)
    if interesting:
      rg=two.get_residue_group_from_hierarchy(hierarchy, index)
      for atom in rg.atoms(): break
      keys.setdefault(tuple(interesting), [])
      id_str = atom.parent().id_str()
      if two.altloc():
        id_str=two.altloc()+id_str[1:]
      if two.conformer:
        id_str='"%s"' % two.conformer+id_str
      if id_str not in keys[tuple(interesting)]:
        keys[tuple(interesting)].append(id_str)
      if verbose: print('  %s : "%s" |%s|' %(one, rc, interesting))
  #
  cache={}
  cache['start']={}
  cache['middle']={}
  cache['end']={}
  atoms=hierarchy.atoms()
  # geometry=grm.geometry
  # from mmtbx.conformation_dependent_library import generate_protein_threes
  import copy
  from mmtbx.conformation_dependent_library import generate_residue_tuples
  keys={}
  for dna_rna_residues in [False, True]:
    for two in generate_residue_tuples(hierarchy=hierarchy,
                                       geometry=None,
                                       length=2,
                                       backbone_only=False,
                                       dna_rna_residues=dna_rna_residues,
                                       include_non_linked=True,
                                       ):
      if verbose: print(two)
      are_linked = two.are_linked(return_value=True)
      one=copy.copy(two)
      del one[1]
      post_process=False
      if are_linked not in [True]:
        # print(two)
        post_process=True
        one.end=True
      elif one.start and one.end:
        assert 0
        if are_linked in [True]:
          one.end=None
        else:
          print(two)
          assert 0
        post_process=True
        assert 0
      elif one.end:
        one.end=None
        post_process=True

      _process_rg_for_compute(one, two, index=0)
      if post_process:
        one=copy.copy(two)
        if two.start and two.end:
          assert 0
          if are_linked in [True]:
            one.start=None
          else:
            print(two)
            assert 0
        if are_linked not in [True]:
          two.end=True
        del one[0]
        _process_rg_for_compute(one, two, index=1)
  # think of an object to send back
  if verbose:
    for part, item in cache.items():
      print(part)
      for an in item:
        print(an)
  return keys

def print_voids(rc):
  print()
  print('='*80)
  print('\n  Missing atoms analysis (comparison of model atom names and restraints)')
  if rc:
    errors=rc.get('errors', {})
    warnings=rc.get('warnings', {})
    waters=rc.get('water',[])
    comments=rc.get('comments', {})
    incomplete=rc.get('incomplete', {})
    remove=[]
    for key, item in warnings.items():
      if not item:
        remove.append(key)
    for key in remove: del warnings[key]
    if incomplete:
      print('\n    Gross errors')
      for key, item in incomplete.items():
        print('      %d %-25s %s' % (len(item), key, print_list(item)))
    if errors:
      print('\n    ERRORS')
      for resname, item in errors.items():
        print('    %3d %s residue has missing atoms'% (len(item), resname))
        for i, (resi, missing_heavies) in enumerate(item):
          print('         %s : %s' % (resi,','.join(missing_heavies)))
    if warnings:
      print('\n    Warnings')
      for resname, item in warnings.items():
        if item:
          print('    %3d %s residue has missing H' % (len(item), resname))
          for i, (resi, missing_hs) in enumerate(item):
            print('         %s : %s' % (resi,','.join(missing_hs)))
    if waters:
      print('\n    Missing water H : %d ~> %s' % (len(waters), print_list(waters)))
    if comments:
      print('\n    Comments')
      for resname, item in comments.items():
        print('    %3d %s residue has missing atoms'% (len(item), resname))
        for i, (resi, missing_heavies) in enumerate(item):
          print('         %s : %s' % (resi,','.join(missing_heavies)))
  else:
    print('\n  No missing atoms found.')

  print('  ...\n')

def print_types(rc):
  print('='*80)
  print('\n  Residue types analysis')
  if rc:
    for key, item in rc.items():
      outl = '    %3d ' % (len(item))
      key = '\|'.join(key)
      outl += '%-22s ~> %s' % (key, print_list(item))
      print(outl)

  print('  ...\n')

def Nietzsche(void):
  # looking into the void
  errors={}
  warnings={}
  comments={}
  incomplete={}
  waters=[]
  for resi, item in void.items():
    resc = item.get('class', None)
    resname = item.get('resname')
    if resc == 'water':
      h=item.get('missing',{})
      h=h.get('hydrogens',{})
      if list(h.keys())==['H1', 'H2']:
        waters.append(resi)
    elif resc == 'peptide':
      if 'incomplete' in item:
        incomplete.setdefault(item['incomplete'], [])
        incomplete[item['incomplete']].append(resi)
        continue
      h=item.get('missing',{})
      c=h.get('heavy',{})
      h=h.get('hydrogens',{})
      if c:
        errors.setdefault(resname, [])
        errors[resname].append([resi, list(c.keys())])
      if resname == 'CYS':
        warnings.setdefault(resname, [])
        if list(h.keys())==['HG']:
          warnings[resname].append([resi, list(h.keys())])
        else:
          errors.setdefault(resname, [])
          errors[resname].append([resi, list(h.keys())])
      elif resname == 'HIS':
        warnings.setdefault(resname, [])
        if list(h.keys())==['HD1', 'HE2']:
          warnings[resname].append([resi, list(h.keys())])
      else:
        if h:
          warnings.setdefault(resname, [])
          warnings[resname].append([resi, list(h.keys())])
    else:
      h=item.get('missing',{})
      c=h.get('heavy',{})
      h=h.get('hydrogens',{})
      comments.setdefault(resname, [])
      comments[resname].append([resi, list(c.keys())+list(h.keys())])
  rc={}
  rc['water']=waters
  rc['warnings']=warnings
  rc['errors']=errors
  rc['comments']=comments
  rc['incomplete']=incomplete
  return rc

def _generate_bonds_with_origin_ids_in_list(bond_proxies, specific_origin_ids=None):
  assert specific_origin_ids
  for specific_origin_id in specific_origin_ids:
    for p in bond_proxies.get_proxies_with_origin_id(specific_origin_id):
      yield p

def process_voids(v, model):
  grm=model.get_restraints_manager()
  sg_bonds=[]
  specific_origin_ids = [origin_ids.get_origin_id('SS BOND'),
                         # origin_ids.get_origin_id('metal coordination'),
                         # origin_ids.get_origin_id('Misc. bond'),
    ]
  for bond in grm.geometry.get_all_bond_proxies():
    if not hasattr(bond, 'get_proxies_with_origin_id'): continue
    for p in _generate_bonds_with_origin_ids_in_list(bond, specific_origin_ids):
      sg_bonds.append(p.i_seqs)
  atoms=model.get_hierarchy().atoms()
  bound_sg=[]
  for i, sgb in enumerate(sg_bonds):
    bound_sg.append(atoms[sgb[0]].parent().id_str().strip())
    bound_sg.append(atoms[sgb[1]].parent().id_str().strip())
  for key, item in v.items():
    if key=='warnings':
      for resname, issues in item.items():
        if resname=='CYS':
          remove=[]
          for i, iss in enumerate(issues):
            if iss[1]==['HG']:
              if iss[0].find('"')>-1:
                tmp=iss[0].split('"')[-1]
              else:
                tmp=iss[0]
              if tmp.strip() in bound_sg:
                remove.append(i)
          if remove:
            remove.reverse()
            for r in remove:
              del issues[r]
  return v

def rmsz(internal):
  if len(internal)==12:
    values=6
  elif len(internal)==7:
    values=3
  delta=abs(internal[values])
  return delta/internal[values+1]

def display_internals(internal):
  if len(internal)==12:
    atoms=2
  elif len(internal)==7:
    atoms=0
  outl=''
  tmp=[]
  for atom in internal[atoms]:
    tmp.append(' %s' % atom.split('"')[1])
  outl+=' -'.join(tmp)
  outl+=' ::  %7.3f' % (internal[atoms+1])
  outl+=' %7.3f' % (internal[atoms+2])
  outl+=' rmsZ: %4.1f' % (rmsz(internal))
  return outl

def process_restraints(model, limit=3):
  grm=model.get_restraints_manager()
  ph=model.get_hierarchy()
  atoms=ph.atoms()
  sites_cart=atoms.extract_xyz()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  bonds=[]
  for bp in grm.geometry.get_all_bond_proxies():
    if not hasattr(bp, 'get_proxies_with_origin_id'): continue
    for bl in bp.get_sorted(by_value="residual",
                            sites_cart=sites_cart,
                            site_labels = site_labels):
      if type(bl)==type(1): continue
      for bond in bl:
        if rmsz(bond)>limit:
          bonds.append(bond)
        else:
          break
  angles=[]
  for ap in grm.geometry.angle_proxies.get_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=site_labels,
    ):
    if type(ap)==type(1): continue
    for angle in ap:
      if rmsz(angle)>limit:
        angles.append(angle)
      else:
        break
  return bonds, angles

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = master_phil_str

  datatypes = ['model','phil']
  # data_manager_options = ['model_skip_expand_with_mtrix']
  # known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    model=self.data_manager.get_model()
    if not model.has_hd():
      raise Sorry('Model must have Hydrogen atoms')

  def run(self):
    self.results={}
    model = self.data_manager.get_model()
    # p = m.get_default_pdb_interpretation_params()
    # p.pdb_interpretation.nonbonded_distance_cutoff = nonbonded_distance_cutoff
    model.process(make_restraints=True) #pdb_interpretation_params=p

    missing_atoms_rc = model._missing_atoms
    v=Nietzsche(missing_atoms_rc)
    process_voids(v, model) #.get_restraints_manager())
    if not self.params.quiet: print_voids(v)
    self.results['missing_atoms']=v

    if self.params.display_peptide_types:
      hierarchy = model.get_hierarchy()
      chains={}
      for chain in hierarchy.chains(): chains[chain.id]=1
      grm = model.get_restraints_manager()
      t = compute(hierarchy, grm)
      if not self.params.quiet: print_types(t)
      self.results['types']=t

    if self.params.display_rmsz_outliers:
      limit=self.params.rmsz_limit
      show=5
      print('='*80)
      print('\n  Internal outliers')
      bonds, angles = process_restraints(model, limit=limit)
      if bonds:
        print('    Bond rmsZ over %4.1f' % limit)
        for i, bond in enumerate(bonds):
          print('    %s' % display_internals(bond))
          if i+1==show:
            if len(bonds)-show:
              print('      Plus %d more outliers' % (len(bonds)-show))
            break
      if angles:
        print('    Angle rmsZ over %4.1f' % limit)
        for i, angle in enumerate(angles):
          print('    %s' % display_internals(angle))
          if i+1==show:
            if len(angles)-show:
              print('      Plus %d more outliers' % (len(angles)-show))
            break
      print('  ...')
      self.results['rmsz']={'bonds': bonds, 'angles': angles}

    if self.params.legend and not self.params.quiet:
      legend = [
        '',
        '='*80,
        '',
        '  Residue legend',
        '',
        '    Eacn residue is identified by the internal conformer. For each altloc',
        '    letter used in the model, Phenix will generate an internal copy of the',
        '    entire chain. The conformer is added in double quotes while the altloc',
        '    is before the residue code. So:',
        '',
        '      "A"ASER A   1',
        '',
        '    is a SER which happens to have an altloc A in the conformer "A". Each',
        '    is added when detected.'
        ]
      print('\n'.join(legend))

  def get_results(self):
    return self.results
