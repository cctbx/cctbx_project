from __future__ import absolute_import, division, print_function
from cctbx.geometry_restraints.auto_linking_types import origin_ids
from cctbx.geometry_restraints.auto_linking_types import covalent_headers

class linking_class(dict):
  def __init__(self):
    self.data = {}
    origin_id = 0
    for oi in origin_ids:
      for i, item in oi.items():
        if item[0] in self: continue
        self[item[0]] = i #origin_id
        self.data[item[0]] = item
        origin_id+=1

  def __repr__(self):
    outl = 'links\n'
    for i, item in sorted(self.items()):
      if type(i)==type(0) and 0:
        outl += '  %3d\n' % (i)
        for j in item:
          outl += '      %s\n' % (j)
      else:
        outl += '  %-20s : %s\n' % (i, item)
    return outl

  def __getitem__(self, key):
    try:
      return dict.__getitem__(self, key)
    except KeyError as e:
      print('''
Look for a key in the list below

%s
      ''' % self)
      raise e

  def get_origin_id(self,
                    key,
                    add_if_absent=False,
                    return_none_if_absent=False,
                    ):
    rc = self.get(key, None)
    if not rc:
      if return_none_if_absent:
        pass
      elif add_if_absent:
        rc = max(self.values())+1 # = origins(key, internals=[0,1,2,3,4,5])
        self[key] = rc
        self.data[key] = ['key', 'user supplied link']
      else:
        assert rc is not None, 'linking origin id not found for "%s"' % key
    return rc

  def get_origin_key(self, origin_id):
    for key, item in self.items():
      if item==origin_id:
        return key
    return None

  def _get_origin_id_labels(self, internals=None):
    keys = list(self.keys())
    def _sort_on_values(k1, k2):
      if self[k1]<self[k2]: return -1
      return 1
    def _filter_on_internals(k1):
      ptr = {'bonds':0,
             'angles':1,
             'dihedrals':2,
             'planes':3,
             'chirals':4,
             'parallelity':5,
             }[internals]
      if internals is None: return True
      if ptr in self.data[k1].internals: return True
      return False
    from functools import cmp_to_key
    keys.sort(key = cmp_to_key(_sort_on_values))
    keys = filter(_filter_on_internals, keys)
    return keys

  def get_bond_origin_id_labels(self):
    return self._get_origin_id_labels(internals='bonds')

  def get_angle_origin_id_labels(self):
    return self._get_origin_id_labels(internals='angles')

  def get_dihedral_origin_id_labels(self):
    return self._get_origin_id_labels(internals='dihedrals')

  def get_chiral_origin_id_labels(self):
    return self._get_origin_id_labels(internals='chirals')

  def get_plane_origin_id_labels(self):
    return self._get_origin_id_labels(internals='planes')

  def get_parallelity_origin_id_labels(self):
    return self._get_origin_id_labels(internals='parallelity')

  def get_geo_file_header(self, origin_id_label, internals=None):
    info = self.data.get(origin_id_label, None)
    assert info
    if len(info)>=4:
      rc = info[3]
      assert type(rc)==type([])
      if internals in [None, 'bonds']: return rc[0]
      elif internals in ['angles']: return rc[1]
      elif internals in ['dihedrals']: return rc[2]
      elif internals in ['chirals']: return rc[3]
      elif internals in ['planes']: return rc[4]
      elif internals in ['parallelities']: return rc[5]
      else: assert 0
    else: return info[0]

  def parse_geo_file_header(self, origin_id_label, subheader=None, internals=None):
    if origin_id_label in ['Nonbonded']:
      # special case for Nonbonded
      return 0, origin_id_label
    elif not origin_id_label in covalent_headers:
      assert 0, 'origin_id_label "%s" not in %s' % (origin_id_label, covalent_headers)
    info = self.data.get(origin_id_label, None)
    if info:
      assert 0
    else:
      for i, (origin_label, info) in enumerate(self.data.items()):
        if len(info)==2:
          if subheader in info:
            return i, '%s %s' % (origin_id_label, subheader)
        elif len(info)>=4:
          header_info = info[3]
          if isinstance(header_info, list):
            if subheader in header_info:
              return i, '%s %s' % (origin_id_label, subheader)

    print(origin_id_label, subheader)
    assert 0

  def get_origin_label_and_internal(self, query_header, verbose=False):
    '''
    Input:
      Line from a .geo file
    Returns:
      None if not a header line
      If header:
      tuple of:
        Origin_id :
        Bond type : Bond, Bond angle, ...
        Subtype   : A string related mostly to the cif_link used but also
                    origin_id
        Number of restraints: int - number from the header.
    '''
    if verbose:
      for origin_label, info in self.data.items():
        print('origin_label, info',origin_label,info)
    if query_header.find('|')==-1: return None
    tmp = query_header.split('|')
    header=tmp[0].strip()
    subheader=tmp[1].strip()
    oi, rc = self.parse_geo_file_header(header, subheader=subheader)
    rc = rc.replace(header, '').strip()
    tmp = query_header.split(':')
    num = int(tmp[-1])
    return oi, header, rc, num

if __name__=='__main__':
  lc = linking_class()
  print(lc)
  for line in [ 'Bond | covalent geometry | restraints: -1',
                'Bond | Misc. | restraints: -1',
                'Bond | link_BETA1-4 | restraints: -1',
                'Bond | link_TRANS | restraints: -1',
                'Bond angle | covalent geometry | restraints: -1',
                'Bond angle | link_BETA1-4 | restraints: -1',
                'Bond angle | link_TRANS | restraints: -1',
                'Dihedral angle | covalent geometry | restraints: -1',
                'Dihedral angle | C-Beta improper | restraints: -1',
                'Dihedral angle | link_TRANS | restraints: -1',
                'Chirality | covalent geometry | restraints: -1',
                'Chirality | link_BETA1-4 | restraints: -1',
                'Planarity | covalent geometry | restraints: -1',
                'Planarity | link_TRANS | restraints: -1',

                "Bond | Bond-like | restraints: -1",
                "Bond angle | Secondary Structure restraints around h-bond | restraints: -1",
                "Parallelity | Stacking parallelity | restraints: -1",
                "Parallelity | Basepair parallelity | restraints: -1",
                'random line',

                #148L
                'Bond | covalent geometry | restraints: 1390',
                'Bond | Misc. | restraints: 4',
                'Bond | link_BETA1-4 | restraints: 1',
                'Bond | link_TRANS | restraints: 1',
                'Bond angle | covalent geometry | restraints: 1868',
                'Bond angle | link_BETA1-4 | restraints: 3',
                'Bond angle | link_TRANS | restraints: 3',
                'Dihedral angle | covalent geometry | restraints: 563',
                'Dihedral angle | C-Beta improper | restraints: 308',
                'Dihedral angle | link_TRANS | restraints: 3',
                'Chirality | covalent geometry | restraints: 210',
                'Chirality | link_BETA1-4 | restraints: 1',
                'Planarity | covalent geometry | restraints: 238',
                'Planarity | link_TRANS | restraints: 1',
                'Nonbonded | unspecified | interactions: 15086',
    ]:
    print('.........',line, lc.get_origin_label_and_internal(line))
