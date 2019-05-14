from __future__ import absolute_import, division, print_function
class bond_angle_registry(dict):
  def __init__(self):
    pass

  def __repr__(self):
    outl = "CDL Bond Angle Registry"
    for i, residue in enumerate(sorted(self)):
      outl += "\n  %s : %7.2f" % (str(residue),
                                  self[residue].angle_ideal,
                                  )
    return outl

if __name__=="__main__":
  pass
