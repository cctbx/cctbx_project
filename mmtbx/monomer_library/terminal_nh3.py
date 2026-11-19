"""
N-terminal NH3 for peptide monomers whose dictionary does not fit mod_NH3:
no atom H (free-form dictionaries: HN1/HN2, HN/HNA, ...), or an atom named
H1, H2 or H3 on another atom than N (2KZ H2 on CB, NPH H3 on C3). Such a
terminus keeps the dictionary's N-H names and gets the missing H under names
computed from the dictionary alone (extra_h_names()); pdb_interpretation
restrains it with a generated modification (nh3_mod()), on every read.
Shared by pdb_interpretation and reduce2 so that both use the same names.
"""
from __future__ import absolute_import, division, print_function
import re
from mmtbx.monomer_library import cif_types

FALLBACK = ('HN3', 'HNT', 'HT3', 'H3N', 'HNB', 'HNC')
STANDARD = ('H1', 'H2', 'H3')

def _types(monomer):
  return dict((a.atom_id, (a.type_symbol or '').strip().upper())
    for a in monomer.atom_list)

def _neighbours(monomer, name):
  result = []
  for b in monomer.bond_list:
    if b.atom_id_1 == name: result.append(b.atom_id_2)
    elif b.atom_id_2 == name: result.append(b.atom_id_1)
  return result

def n_h_names(monomer):
  """The dictionary's H on N, sorted by length and name."""
  t = _types(monomer)
  return sorted((x for x in _neighbours(monomer, 'N')
    if t.get(x) in ('H', 'D')),
    key=lambda x: (len(x), x))

def standard_names_collide(monomer):
  """An atom named H1, H2 or H3 that the dictionary bonds to another atom
  than N: mod_NH3 would take its name."""
  for name in STANDARD:
    if monomer.i_atom_by_id(name) is None: continue
    if 'N' not in _neighbours(monomer, name): return True
  return False

def mod_takes_name(monomer, mod_mod_id):
  """A modification adds an atom whose name the monomer has on another atom
  than N (an H already on N is replaced, as for CCD-style NH2)."""
  for a in mod_mod_id.atom_list:
    if a.function != 'add': continue
    if monomer.i_atom_by_id(a.new_atom_id) is None: continue
    if 'N' not in _neighbours(monomer, a.new_atom_id): return True
  return False

def _bump(n):
  # digits only (HN1 HN2 -> HN3, H3 H4 -> H5), or the letter after HN
  # (HN HNA -> HNB); never a bare letter (HN -> HO)
  m = re.match(r'^(\D*)(\d+)(\D*)$', n)
  if m: return '%s%d%s' % (m.group(1), int(m.group(2)) + 1, m.group(3))
  m = re.match(r'^HN([A-Y])$', n)
  if m: return 'HN' + chr(ord(m.group(1)) + 1)
  return None

# N-CA at least this long (A) is single where the dictionary gives no bond
# order (mon_lib "coval"): imines 1.27-1.30, 2KZ's enamine N-CA 1.402
N_CA_SINGLE_MIN = 1.35

def n_ca_single(monomer):
  """Whether the dictionary's N-CA is a single bond: by its type, or, where
  the type gives no order (coval, missing), by the ideal length
  (N_CA_SINGLE_MIN); without either, not single. An imine N (HHI, ITR: N=CA 1.265) is not an amine."""
  for b in monomer.bond_list:
    if set([b.atom_id_1, b.atom_id_2]) != set(['N', 'CA']): continue
    t = (b.type or '').strip().lower()
    if t in ('single', 'sing'): return True
    if t in ('double', 'doub', 'triple', 'trip', 'aromatic', 'arom',
        'deloc', 'metal'):
      return False
    if b.value_dist is None: return False
    return b.value_dist >= N_CA_SINGLE_MIN
  return False

def applies(monomer):
  """Whether a terminus of this monomer is built and restrained here:
  a peptide whose N has one heavy neighbour, CA, by a single bond
  (n_ca_single()), and either no atom H or an H1/H2/H3 not on N."""
  if not monomer.is_peptide(): return False
  t = _types(monomer)
  heavy = [x for x in _neighbours(monomer, 'N') if t.get(x) not in ('H', 'D')]
  if heavy != ['CA']: return False
  if not n_ca_single(monomer): return False
  return monomer.i_atom_by_id('H') is None or standard_names_collide(monomer)

def extra_h_names(monomer):
  """Names of the H that make the dictionary's N-H a terminal NH3: the
  dictionary's own pattern for one missing H; after a lone H the free ones
  of H2, H3 (the series of the CCD and old reduce; NPH: H2, HN3); then the
  first free fallback names. [] if none is missing; None if no free name
  is left."""
  hs = n_h_names(monomer)
  need = 3 - len(hs)
  if need <= 0: return []
  taken = set(a.atom_id for a in monomer.atom_list)
  if need == 1 and hs:
    b = _bump(hs[-1])
    if b and len(b) <= 4 and b.startswith('H') and b not in taken:
      return [b]
  got = []
  if hs == ['H']:
    # never H1: H with H1 fits no convention
    got = [x for x in ('H2', 'H3') if x not in taken]
  got += [f for f in FALLBACK if f not in taken][:need - len(got)]
  if len(got) < need: return None
  return got

def ca_sp3(monomer):
  """CA with four dictionary neighbours (2KZ, CA=CB: three). By neighbour
  count, not bond type: mon_lib-format entries type every bond coval."""
  return len(_neighbours(monomer, 'CA')) == 4

# mod_NH3's N-H (x-ray, esd, neutron), where no mod_NH3 is given
NH3_N_H = (0.89, 0.02, 1.04)

# an N-H this short (A) is an x-ray length (GeoStd 0.86-0.89; nuclear ~1.02)
XRAY_N_H_MAX = 0.95

def n_h_values(monomer, mod_nh3=None):
  """The one N-H value set (x-ray, esd, neutron) of all three terminal N-H,
  and whether it is the dictionary's own: the dictionary's N-H if every N-H
  on N has the same ideal, esd and neutron value (a consistent custom CIF
  keeps its values); otherwise mod_NH3's (TSY 0.86/0.89; NPH without a
  neutron value: its neutron column repeats the x-ray 0.86, shorter than
  XRAY_N_H_MAX)."""
  hs = n_h_names(monomer)
  own = set((b.value_dist, b.value_dist_esd, b.value_dist_neutron)
    for b in monomer.bond_list
    if 'N' in (b.atom_id_1, b.atom_id_2) and
      set([b.atom_id_1, b.atom_id_2]) & set(hs))
  if len(own) == 1:
    x = list(own)[0]
    # no neutron value: missing, or an x-ray length repeated (NPH 0.86);
    # one nuclear length in both columns (1.03/1.03) is kept
    if None not in x and not (x[2] == x[0] and x[0] < XRAY_N_H_MAX):
      return x, True
  for b in (mod_nh3.bond_list if mod_nh3 is not None else []):
    if (b.function == 'add' and 'N' in (b.atom_id_1, b.atom_id_2) and
        None not in (b.new_value_dist, b.new_value_dist_esd,
          b.new_value_dist_neutron)):
      return (b.new_value_dist, b.new_value_dist_esd,
        b.new_value_dist_neutron), False
  return NH3_N_H, False

def nh3_mod(monomer, extra, mod_nh3=None):
  """mod_NH3's geometry over the dictionary's N-H names plus extra: N type
  NT3, each extra H bonded to N; the three N-H with one value set
  (n_h_values(): the dictionary's own N-H if they agree, else mod_NH3's),
  every H-N-CA and H-N-H angle 109.47; N-CA as mod_NH3 (1.491) where CA is
  sp3 (ca_sp3()), else the dictionary's value (an sp2 CA, 2KZ: open)."""
  comp = monomer.chem_comp.id
  hs = n_h_names(monomer)
  (dist, esd, neutron), own = n_h_values(monomer, mod_nh3)
  mod = cif_types.mod_mod_id(source_info='generated by terminal_nh3',
    chem_mod=cif_types.chem_mod(id='NH3-%s' % comp,
      name='NH3 amino terminus, %s dictionary names' % comp, comp_id='.',
      group_id='.'))
  mod.atom_list.append(cif_types.chem_mod_atom(function='change',
    atom_id='N', new_atom_id='N', new_type_symbol='',
    new_type_energy='NT3', new_partial_charge=-0.2))
  for h in extra:
    mod.atom_list.append(cif_types.chem_mod_atom(function='add', atom_id='',
      new_atom_id=h, new_type_symbol='H', new_type_energy='HNT3',
      new_partial_charge=0.1))
    mod.bond_list.append(cif_types.chem_mod_bond(function='add',
      atom_id_1='N', atom_id_2=h, new_type='single', new_value_dist=dist,
      new_value_dist_esd=esd, new_value_dist_neutron=neutron))
  for h in ([] if own else hs):
    mod.bond_list.append(cif_types.chem_mod_bond(function='change',
      atom_id_1='N', atom_id_2=h, new_type='single', new_value_dist=dist,
      new_value_dist_esd=esd, new_value_dist_neutron=neutron))
  for b in (mod_nh3.bond_list if mod_nh3 is not None and ca_sp3(monomer)
      else []):
    if b.function == 'change' and set([b.atom_id_1, b.atom_id_2]) == set(
        ['N', 'CA']):
      mod.bond_list.append(cif_types.chem_mod_bond(function='change',
        atom_id_1='N', atom_id_2='CA', new_type=b.new_type,
        new_value_dist=b.new_value_dist,
        new_value_dist_esd=b.new_value_dist_esd,
        new_value_dist_neutron=b.new_value_dist_neutron))
  every = hs + list(extra)
  existing = set()
  for a in monomer.angle_list:
    if a.atom_id_2 == 'N':
      existing.add(frozenset([a.atom_id_1, a.atom_id_3]))
  pairs = [(h, 'CA') for h in every] + [(every[i], every[j])
    for i in range(len(every)) for j in range(i + 1, len(every))]
  for x, y in pairs:
    mod.angle_list.append(cif_types.chem_mod_angle(
      function='change' if frozenset([x, y]) in existing else 'add',
      atom_id_1=x, atom_id_2='N', atom_id_3=y, new_value_angle=109.47,
      new_value_angle_esd=3.0))
  return mod

def _h_forms(names):
  # a D conformer (neutron H/D pairs) names them DN3 for HN3
  result = set(names)
  result.update('H' + x[1:] for x in names if x[:1] == 'D')
  return result

def terminal_mod(monomer, unexpected, mod_nh3=None):
  """The generated modification for a terminus of monomer, if any extra
  name is among the unexpected atoms; else None. Built per call from the
  monomer as it is (about 0.05 ms), so a replacement CIF counts.
  mod_nh3: the library's mod_NH3 (its N-CA term)."""
  if not applies(monomer): return None
  extra = extra_h_names(monomer)
  # any of them, as mod_NH3 takes any of H1-H3 (a partial terminus)
  if not extra or not set(extra) & _h_forms(unexpected): return None
  return nh3_mod(monomer, extra, mod_nh3)

# the names pdb_interpretation's N-terminal rules test (mod_NH3, NH2, NH1)
RULE_NAMES = ('H1', 'H2', 'H3', '1H', '2H', '3H', 'HN1', 'HN2', '1HN',
  '2HN', 'HN')

def computed_names_mod(monomer, unexpected, mod_nh3=None):
  """terminal_mod() where every unexpected name that pdb_interpretation's
  N-terminal rules test (RULE_NAMES) is a computed name, so that the rules
  do not take it (0AA, CTE: computed HN2 -> mod_NH2); else None."""
  mod = terminal_mod(monomer, unexpected, mod_nh3)
  if mod is None: return None
  if (_h_forms(unexpected) & set(RULE_NAMES)) - set(extra_h_names(monomer)):
    return None
  return mod
