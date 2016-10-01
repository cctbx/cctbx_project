from __future__ import division
from mmtbx.conformation_dependent_library import cdl_database

errors = {'Pro_nonxpro': [
  (-90, 60),
  (-90, 70),
  (-90, 80),
  (-80, 60),
  (-80, 70),
  (-80, 80),
  ],
  }

def run():
  bond_esds = []
  angle_esds = []
  for res_group_type in cdl_database:
    print res_group_type, len(cdl_database[res_group_type])
    e = errors.get(res_group_type, [])
    cd = cdl_database[res_group_type]
    for phi in range(-180,180,10):
      for psi in range(-180,180,10):
        if (phi, psi) in e: continue
        restraints = cd.get((phi, psi), None)
        if restraints is None: continue
        for i in range(2,len(restraints),2):
          if restraints[i]==-1: continue # GLY has no CB
          print res_group_type, phi, psi, i, restraints[i],restraints[i+1]
          if restraints[i]>2: # angles
            angle_esds.append(restraints[i+1])
            print 'angle',min(angle_esds), max(angle_esds)
            assert restraints[i+1]<5, "angle esd too large"
            assert restraints[i+1]>.9, "angle esd too small"
          else:
            bond_esds.append(restraints[i+1])
            print 'bond',min(bond_esds),max(bond_esds)
            assert restraints[i+1]<0.05, "bond esd too large"
            assert restraints[i+1]>0.001, "bond esd too small"

if __name__=="__main__":
  run()#sys.argv[1])
