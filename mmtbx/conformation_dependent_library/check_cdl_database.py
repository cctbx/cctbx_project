from __future__ import absolute_import, division, print_function
import os, sys
from mmtbx.conformation_dependent_library import cdl_database
from mmtbx.conformation_dependent_library import cdl_setup
from six.moves import range

errors = {'Pro_nonxpro': [
  (-90, 60),
  (-90, 70),
  (-90, 80),
  (-80, 60),
  (-80, 70),
  (-80, 80),
  ],
  }

print_to_disk = " -hdevice JPEG -hardcopy"

def generate_restraints():
  for res_group_type in cdl_database:
    print(res_group_type, len(cdl_database[res_group_type]))
    e = errors.get(res_group_type, [])
    cd = cdl_database[res_group_type]
    for phi in range(-180,180,10):
      for psi in range(-180,180,10):
        #if (phi, psi) in e: continue
        restraints = cd.get((phi, psi), None)
        if restraints is None: continue
        yield res_group_type, phi, psi, restraints

def find_esd_extrema():
  bond_esds = []
  angle_esds = []
  for res_group_type, phi, psi, restraints in generate_restraints():
    for i in range(2,len(restraints),2):
      if restraints[i]==-1: continue # GLY has no CB
      print(res_group_type, phi, psi, i, restraints[i],restraints[i+1])
      if restraints[i]>2: # angles
        angle_esds.append(restraints[i+1])
        print('angle',min(angle_esds), max(angle_esds))
        assert restraints[i+1]<5, "angle esd too large"
        assert restraints[i+1]>.9, "angle esd too small"
      else:
        bond_esds.append(restraints[i+1])
        print('bond',min(bond_esds),max(bond_esds))
        assert restraints[i+1]<0.05, "bond esd too large"
        assert restraints[i+1]>0.001, "bond esd too small"

def analysis_esd(number_of_observations_max=20,
                 maximum_z_score=3,
                 minimum_z_score=0.333,
                 ):

  from elbow.utilities import rmsd_utils

  outl = ""
  py = ""

  defaults = {}
  data = {}
  for res_group_type, phi, psi, restraints in generate_restraints():
    defaults.setdefault(res_group_type, None)
    if restraints[0]=='I':
      defaults[res_group_type]=restraints
      continue

  for res_group_type, phi, psi, restraints in generate_restraints():
    if restraints[0]=='I': continue
    print(res_group_type, phi, psi, restraints)
    data.setdefault(res_group_type, {})
    data[res_group_type].setdefault(restraints[1], {})
    for i in range(2, len(restraints), 2):
      if restraints[i]==-1: continue # GLY has no CB
      print(res_group_type, phi, psi, i, restraints[i],restraints[i+1])
      if restraints[i]>2: # angles
        ptr = 'angles'
      else:
        ptr = 'bonds'
      data[res_group_type][restraints[1]].setdefault(ptr,[])
      d = defaults[res_group_type][i+1]
      z = restraints[i+1]/d
      print(i,z,d,restraints[i+1])
      if z>=maximum_z_score or z<=minimum_z_score:
        for j in range(i,i+2):
          py += "%s[(%d, %d)][%d]=%f # %s\n" % (
            res_group_type,
            phi,
            psi,
            j,
            defaults[res_group_type][j],
            cdl_setup.columns[j],
            )
        outl += " %-18s %4d %4d %2d %2d %-4s %9.4f %7.4f %9.4f %7.4f %5.2f\n" % (
          res_group_type,
          phi,
          psi,
          i,
          restraints[1],
          cdl_setup.columns[i],
          restraints[i],
          restraints[i+1],
          defaults[res_group_type][i],
          defaults[res_group_type][i+1],
          z,
          )
      data[res_group_type][restraints[1]][ptr].append(z)
    if 'Pro_nonxpro' in data:
      if 3 in data['Pro_nonxpro']:
        print('3'*80)
        print(data['Pro_nonxpro'][3])
        #assert 0
    #if len(data)>1: break
  print(data)
  print(outl)
  print(py)
  if 0: xmgace_data(data)

def xmgace_data(data):
  for res_group_type in data:
    for i in range(1,20):
      print(i)
      outl = ""
      d = data[res_group_type].get(i, None)
      if d is None: continue
      print(d)
      for ptr in ['bonds', 'angles']:
        h = rmsd_utils.histogram(d[ptr])
        print(h)
        outl += "@type xy\n"
        for z,c in sorted(h.items()):
          print(z,c)
          outl += " %f %f\n" % (z,c)
        outl += "&\n"
      print(outl)
      df = "%s_%02d.dat" % (res_group_type, i)
      f=open(df, "w")
      f.write(outl)
      f.close()
      cmd  = 'xmgrace -geometry 1100x900 -param z.par %s' % df
      cmd += """ -pexec 'title "%s"'""" % res_group_type
      cmd += """ -pexec 'subtitle "%s"'""" % 'Sample size: %d' % i
      os.system(cmd)
      cmd += print_to_disk
      print(cmd)
      os.system(cmd)
      if os.path.exists("Untitled.jpg"):
        os.rename("Untitled.jpg", "%s" % df.replace(".dat",".jpg"))

def run():
  if 0: find_esd_extrema()
  if 1: analysis_esd()

if __name__=="__main__":
  run()#sys.argv[1])
