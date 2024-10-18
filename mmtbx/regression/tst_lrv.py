from __future__ import division
import os
import libtbx.load_env
from libtbx import easy_run

repository_dir = libtbx.env.dist_path("chem_data")

def main(code='NWM'):
  print(repository_dir)
  rf = os.path.join(repository_dir,
                    'geostd',
                    code.lower()[0],
                    'data_%s.cif' % code)
  assert os.path.exists(rf), 'path %s not found' % rf

  for args in ['',
               'write_pdb=lrv_test.pdb',
               'write_geo=lrv_test.geo',
               'write_pdb=lrv_test.pdb write_geo=lrv_test.geo',
               'use_hydrogens=False',
               ]:
    cmd = 'mmtbx.development.ligand_restraints_validation %s' % rf
    cmd += ' %s' % args
    print(cmd)
    rc=easy_run.go(cmd)
    for line in rc.stdout_lines:
      if line.find('RMSD')>-1: break
    else:
      rc.show_stdout()
      assert 0

if __name__ == '__main__':
  main()
