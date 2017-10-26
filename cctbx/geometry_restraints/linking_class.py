from cctbx.geometry_restraints.auto_linking_types import origin_ids

class linking_class(dict):
  def __init__(self):
    for i, item in origin_ids[0].items(): # bond
      self[i]=item
      self[item[0]] = i
      
  def __repr__(self):
    outl = 'links\n'
    for i, item in self.items():
      if type(i)==type(0) and 0:
        outl += '  %3d\n' % (i)
        for j in item:
          outl += '      %s\n' % (j)
      else:
        outl += '  %-20s : %s\n' % (i, item)
    return outl

  def get_origin_id(self, key):
    rc = self.get(key, None)
    assert rc is not None, 'linking origin id not found for "%s"' % key
    return rc

if __name__=='__main__':
  lc = linking_class()
  print lc
