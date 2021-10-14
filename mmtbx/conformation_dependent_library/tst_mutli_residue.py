from __future__ import absolute_import, division, print_function
from iotbx import pdb
from mmtbx.conformation_dependent_library.testing_utils import get_geometry_restraints_manager

from mmtbx.regression import model_1yjp
from six.moves import range

from libtbx.test_utils import approx_equal


answers = {'omegas' : [ 171.5156431406436,
                       -176.23617714135813,
                        166.21212041473837,
                        173.18755254778898,
                       -176.72682825711064,
                       -171.77764536433386],
           'calphas' : [-141.35004030187673,
                         169.2753355533696,
                         179.25726892530153,
                        -169.00952867101745],
           'phi_psi' : [ -60.58200175013997,
                         141.18718965243215,
                        -119.25059933081734,
                         125.14691314895126,
                        -126.15928469653304,
                         112.80786803182154,
                        -114.97595140659884,
                         126.75642562595709,
                        -116.41993898593708,
                          97.69496718589357],
}

def main():
  pdb_inp = pdb.input(lines=model_1yjp, source_info='model_1yjp')
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(raw_records=model_1yjp)
  pdb_hierarchy.reset_i_seq_if_necessary()
  from mmtbx.conformation_dependent_library import generate_protein_fragments
  for k in range(2,6):
    for j, threes in enumerate(generate_protein_fragments(
      pdb_hierarchy,
      geometry_restraints_manager,
      length=k,
      #verbose=verbose,
      )):
      i=k-2
      print(i, j, k, threes)
      rc = None
      try: rc = threes.get_omega_value()
      except: print('  omega is not valid') # intentional
      if i>0: assert rc==None
      else:
        print('  omega   %5.1f' % rc)
        assert approx_equal(rc, answers['omegas'][j], 0.1), '%f != %f' % (rc, answers['omegas'][j])
      rc = threes.get_omega_values()
      print('  omegas  %s' % rc)
      assert approx_equal(rc, answers['omegas'][j:j+i+1]), '%s != %s' % (rc,
                                                             answers['omegas'][j:j+i+1]
                                                             )
      rc = None
      try: rc = threes.cis_group()
      except: pass # intentional
      try: print("  cis?    %-5s %s" % (rc, threes.cis_group(limit=30)))
      except: print('  cis? is not valid') # intentional
      if i>=2: assert (rc==None or rc==False), '%s!=%s' % (rc, None)
      else: assert rc == False
      try: print("  trans?  %-5s %s" % (threes.trans_group(), threes.trans_group(limit=30)))
      except: print('  tran? is not valid') # intentional
      print('  cis/trans/twisted? %s' % ' '.join(threes.cis_trans_twisted_list()))
      try: print("  rama    %s" % threes.get_ramalyze_key())
      except: print('  rama not specified') # intentional
      print('  conf    %s' % threes.is_pure_main_conf())
      rc = None
      try: rc = threes.get_phi_psi_angles()
      except: print('  phi/psi not specified') # intentional
      print('  phi/psi %s' % rc)
      if i<1: assert rc==None
      else:
        test = answers['phi_psi'][j*2:(j+i)*2]
        assert approx_equal(rc[0], test[0]), '%s!=%s' % (rc, test)
        assert approx_equal(rc[1], test[1]), '%s!=%s' % (rc, test)
      rc = None
      try: rc = threes.get_ca_dihedrals()
      except: print('  CA dihedrals not specified') # intentional
      print('  CA dihedrals %s' % rc)
      if i<=1: assert rc == None
      else:
        test = answers['calphas'][j:j+i+1-2]
        assert approx_equal(rc[0], test[0]), '%s != %s' % (rc, test)

    print("OK",i+2)


if __name__ == '__main__':
  main()
