from __future__ import division, print_function
from cctbx.geometry_restraints.auto_linking_types import origin_ids

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

  def get_origin_id(self, key):
    rc = self.get(key, None)
    assert rc is not None, 'linking origin id not found for "%s"' % key
    return rc

  def _get_origin_id_labels(self, internals=None):
    keys = self.keys()
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
    keys.sort(_sort_on_values)
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

if __name__=='__main__':
  lc = linking_class()
  print(lc)
