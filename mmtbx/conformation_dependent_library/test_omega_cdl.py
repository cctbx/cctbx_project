import os, sys
import time
from StringIO import StringIO

from mmtbx import conformation_dependent_library as cdl
from mmtbx.conformation_dependent_library import cdl_utils
from mmtbx.conformation_dependent_library.omega_database import omega_database
from mmtbx.conformation_dependent_library import omega

from test_cdl import output_filenames, filenames, get_managers

filenames["cdl_test_1.pdb"][1] = [(130, -100, 110, -70),
                                  (-70, 150, -80),
                                  (-80, 160, -60),
                                  ]
filenames["cdl_test_1.pdb"][2] = [['I', 3528, 178.64, 6.31],
                                  ['I', 20342, 178.94, 6.38],
                                  ['I', 2599, 180.43, 5.63],
                                  ]
filenames["cdl_test_2.pdb"][1] = [(80, 170, -40, -100),
                                  (110, -110, 160, -60),
                                  (80, 170, -40, -100),
                                  (110, -110, 160, -60),
                                  ]
filenames["cdl_test_2.pdb"][2] = [['I', 20342, 178.94, 6.38],
                                  ['I', 20342, 178.94, 6.38],
                                  ['I', 20342, 178.94, 6.38],
                                  ['I', 20342, 178.94, 6.38],
                                  ]
filenames["cdl_test_3.pdb"][1] = [(0, -130, 20, -120),
                                  (-120, 90, -50),
                                  (-50, 130, 90),
                                  (90, 0, -80),
                                  ]
filenames["cdl_test_3.pdb"][2] = [['I', 20342, 178.94, 6.38],
                                  ['I', 928, 180.09, 6.71],
                                  ['I', 1214, 178.49, 6.3],
                                  ['I', 2599, 180.43, 5.63],
                                  ]
filenames["cdl_test_4.pdb"][1] = [(50, -150, 150, -70),
                                  (-70, 140, -60),
                                  (50, -150, 150, -70),
                                  (-70, 140, -60),
                                  ]
filenames["cdl_test_4.pdb"][2] = [['I', 177, 180.02, 5.36],
                                  ['I', 1214, 178.49, 6.3],
                                  ['I', 177, 180.02, 5.36],
                                  ['I', 1214, 178.49, 6.3],
                                  ]
filenames["cdl_test_5.pdb"][1] = [(80, 170, 70, -30),
                                  (-30, 160, -30),
                                  (-30, 160, -180),
                                  (-180, -160, 100),
                                  ]
filenames["cdl_test_5.pdb"][2] = [['I', 100, 183.21, 5.51],
                                  ['I', 29, 178.55, 6.54],
                                  ['I', 1214, 178.49, 6.3],
                                  ['I', 2599, 180.43, 5.63],
                                  ]
filenames["cdl_test_7.pdb"][1] = [(-10, -110, 100, -60),
                                  (-60, -30, -60),
                                  (-60, -30, -130),
                                  (-10, -110, 100, -50),
                                  (-50, -40, -70),
                                  (-70, -20, -110),
                                  ]
filenames["cdl_test_7.pdb"][2] = [['I', 20342, 178.94, 6.38],
                                  ['I', 20342, 178.94, 6.38],
                                  ['I', 20342, 178.94, 6.38],
                                  ['B', 3, 179.78, 7.53],
                                  ['I', 20342, 178.94, 6.38],
                                  ['I', 20342, 178.94, 6.38],
                                  ]
filenames["cdl_test_7.pdb"][3] = [(10, 11, 27, 28)]
filenames["cdl_test_8.pdb"] = [
    [
      "NonPGIV_xpro",
      ],
    [
        (30, -70, 130, -80),
      ],
    [
      ['I', 928, 180.09, 6.71],
      ],
    [
      True,
      ],
    ]

def test_phi_psi_key(hierarchy,
                     filename,
                     restraints_manager,
                     ):
  for i, threes in enumerate(cdl.generate_protein_threes(
    hierarchy,
    #restraints_manager=restraints_manager
    geometry=restraints_manager.geometry,
    )
                             ):
    key = threes.get_cdl_key(force_plus_one=True)
    print key, filenames[filename][1]
    assert key == filenames[filename][1][i]

def test_cdl_lookup(hierarchy,
                    filename,
                    restraints_manager,
                    ):
  for i, threes in enumerate(cdl.generate_protein_threes(
    hierarchy,
    #restraints_manager=restraints_manager
    geometry=restraints_manager.geometry,
    )
                             ):
    res_type_group = cdl_utils.get_res_type_group(
      threes[1].resname,
      threes[2].resname,
      )
    key = threes.get_cdl_key(force_plus_one=True)
    key = key[-2:]
    restraint_values = omega_database[res_type_group][key]
    print i, key, restraint_values[:4], filenames[filename][2]
    del threes.registry.n
    threes.registry.n = {}
    assert restraint_values[:4] == filenames[filename][2][i]

def test_average(hierarchy,
                 filename,
                 restraints_manager,
                 ):
  for i, threes in enumerate(cdl.generate_protein_threes(
    hierarchy,
    geometry=restraints_manager.geometry,
    )
                             ):
    if threes.registry.n: 
      atoms = hierarchy.atoms()
      for key in threes.registry.n:
        print key
        for atom in key:
          print atoms[atom].quote()
      #assert threes.registry.n.keys() == filenames[filename][3]
      assert filenames[filename][3][0] in threes.registry.n

def run_apply(filename, testing=False, verbose=False):
  processed_pdb_file, restraints_manager = get_managers(filename)
  if testing:
    #test_res_type(processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    #              filename,
    #              restraints_manager,
    #              )
    test_phi_psi_key(processed_pdb_file.all_chain_proxies.pdb_hierarchy,
                     filename,
                     restraints_manager,
                     )
    test_cdl_lookup(processed_pdb_file.all_chain_proxies.pdb_hierarchy,
                    filename,
                    restraints_manager,
                    )
    #test_cis_trans_peptide(processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    #                       filename,
    #                       restraints_manager,
    #                       )
    print "OK"

  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart_exact()
  if True:
    lines = StringIO()
    restraints_manager.geometry.show_sorted(sites_cart=sites_cart,
                                                  f=lines)
    f=file("%s_pre.geo" % os.path.basename(filename)[:-4], "wb")
    f.write(lines.getvalue())
    f.close()

  t0=time.time()
  cdl_proxies = omega.setup_restraints(restraints_manager.geometry,
                                       verbose=verbose,
                                       )
  omega.update_restraints(processed_pdb_file.all_chain_proxies.pdb_hierarchy,
                          geometry=restraints_manager.geometry,
                          cdl_proxies=cdl_proxies,
                          verbose=verbose,
                          )
  if testing:
    test_average(processed_pdb_file.all_chain_proxies.pdb_hierarchy,
                 filename,
                 restraints_manager,
                 )
    print 'OK'

  print 'Time to update restraints %0.5fs' % (time.time()-t0)
  if True:
    lines = StringIO()
    restraints_manager.geometry.show_sorted(sites_cart=sites_cart,
                                                  f=lines)
    f=file("%s_post.geo" % os.path.basename(filename)[:-4], "wb")
    f.write(lines.getvalue())
    f.close()

def run(filename=None):
  if filename:
    print 'run oCDL on',filename
    run_apply(filename,
              #verbose=True,
              )
  else:
    print 'Running tests'
    for filename in sorted(output_filenames):
      print filename
      f=file(filename, "wb")
      f.write(output_filenames[filename])
      f.close()
    for filename in sorted(filenames):
      #if filename.find("_7")==-1: continue
      print ('%s ' % filename)*5
      run_apply(filename, testing=True)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
