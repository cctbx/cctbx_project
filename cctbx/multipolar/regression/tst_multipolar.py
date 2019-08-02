from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import os, sys
import libtbx
#import libtbx.load_env
#from libtbx.test_utils import approx_equal
from cctbx import multipolar
from six.moves import range
from six.moves import zip

pdbs = [
  """\
CRYST1   16.372   15.451   19.829  90.00  90.00  90.00 P 1
HETATM    1  N   ALA A   1       5.669   9.679  12.421  1.00 20.00      A    N+1
HETATM    2  H   ALA A   1       5.000   9.066  12.578  1.00 20.00      A    H
HETATM    3  H2  ALA A   1       5.697   9.879  11.525  1.00 20.00      A    H
HETATM    4  H3  ALA A   1       5.501  10.451  12.904  1.00 20.00      A    H
HETATM    5  CA  ALA A   1       6.928   9.127  12.826  1.00 20.00      A    C
HETATM    6  HA  ALA A   1       7.607   9.822  12.759  1.00 20.00      A    H
HETATM    7  CB  ALA A   1       6.843   8.674  14.268  1.00 20.00      A    C
HETATM    8  HB1 ALA A   1       6.188   7.951  14.345  1.00 20.00      A    H
HETATM    9  HB2 ALA A   1       7.720   8.353  14.561  1.00 20.00      A    H
HETATM   10  HB3 ALA A   1       6.567   9.426  14.829  1.00 20.00      A    H
HETATM   11  C   ALA A   1       7.311   7.978  11.888  1.00 20.00      A    C
HETATM   12  O   ALA A   1       6.496   7.106  11.628  1.00 20.00      A    O
HETATM   13  N  AALA A   2       8.599   7.974  11.244  1.00 20.00      A    N
HETATM   14  H  AALA A   2       9.195   8.638  11.416  1.00 20.00      A    H
HETATM   15  CA AALA A   2       8.954   6.905  10.360  1.00 20.00      A    C
HETATM   16  HA AALA A   2       8.104   6.605   9.983  1.00 20.00      A    H
HETATM   17  CB AALA A   2       9.512   5.685  11.065  1.00 20.00      A    C
HETATM   18  HB1AALA A   2       8.788   5.050  11.242  1.00 20.00      A    H
HETATM   19  HB2AALA A   2       9.924   5.957  11.913  1.00 20.00      A    H
HETATM   20  HB3AALA A   2      10.186   5.259  10.496  1.00 20.00      A    H
HETATM   21  C  AALA A   2       9.780   7.353   9.143  1.00 20.00      A    C
HETATM   22  O  AALA A   2      10.743   8.078   9.296  1.00 20.00      A    O
HETATM   13  N  BGLY A   2      -0.366  -0.669   0.828  1.00 20.00      A    N
HETATM   14  H  BGLY A   2      -0.084  -1.525   0.575  1.00 20.00      A    H
HETATM   15  CA BGLY A   2       0.151   0.485   0.139  1.00 20.00      A    C
HETATM   16  HA BGLY A   2      -0.586   1.189   0.119  1.00 20.00      A    H
HETATM   17  HB3BGLY A   2       0.863   0.830   0.654  1.00 20.00      A    H
HETATM   18  C  BGLY A   2       0.668   0.298  -1.316  1.00 20.00      A    C
HETATM   19  O  BGLY A   2       0.556  -0.754  -1.853  1.00 20.00      A    O
HETATM   23  N   ALA A   3       9.411   6.919   7.825  1.00 20.00      A    N
HETATM   24  H   ALA A   3       8.563   6.640   7.660  1.00 20.00      A    H
HETATM   25  CA  ALA A   3      10.371   6.978   6.751  1.00 20.00      A    C
HETATM   26  HA  ALA A   3      11.265   6.846   7.125  1.00 20.00      A    H
HETATM   27  CB  ALA A   3      10.094   5.867   5.758  1.00 20.00      A    C
HETATM   28  HB1 ALA A   3      10.225   5.000   6.194  1.00 20.00      A    H
HETATM   29  HB2 ALA A   3       9.170   5.937   5.439  1.00 20.00      A    H
HETATM   30  HB3 ALA A   3      10.708   5.948   5.000  1.00 20.00      A    H
HETATM   31  C   ALA A   3      10.320   8.345   6.099  1.00 20.00      A    C
HETATM   32  O   ALA A   3      11.372   9.029   5.987  1.00 20.00      A    O
HETATM   33  OXT ALA A   3       9.238   8.772   5.617  1.00 20.00      A    O-1
TER
END
""",
]

def should_have_imported_proper_methods():
  assert hasattr(multipolar, 'multipolar_test' ) == True, 'multipolar has no method multipolar_test'
  assert hasattr(multipolar, 'assign_atom_types' ) == True, 'multipolar has no method assign_atom_types'

def should_return_atom_types():
  atom_number = flex.int()
  coordinates = flex.vec3_double()
  for i in range (3):
    atom_number.append(i)
    coordinates.append( (i,i,i) )
  atom_types = multipolar.assign_atom_types( atom_number, coordinates )
  #print atom_types
  for atom_type in atom_types:
    print(atom_type)

def exercise(file_name=None):
  import libtbx.load_env
  if (not libtbx.env.has_module(name="mmtbx")):
    print("Skipping exercise():", \
      "mmtbx.monomer_library.pdb_interpretation not available")
    return
  if (libtbx.env.find_in_repositories(relative_path="chem_data") is None):
    print("Skipping exercise(): chem_data directory not available")
    return
  from mmtbx import monomer_library
  import mmtbx.monomer_library.server
  import mmtbx.monomer_library.pdb_interpretation
  import mmtbx.restraints
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  #pdb_file = libtbx.env.find_in_repositories(
  #                 relative_path="phenix_regression/pdb/enk.pdb",
  #                test=os.path.isfile)
  if file_name is None:
    pdb_file = "tst_multipoloar.pdb"
    f=open(pdb_file, "w")
    f.write(pdbs[0])
    f.close()

  from iotbx import pdb
  pdb_inp = pdb.input(pdb_file)
  hierarchy = pdb_inp.construct_hierarchy()

  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)

  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  atoms = pdb_hierarchy.atoms()
  atom_altlocs = flex.std_string()
  for atom in atoms:
    atom_altlocs.append(atom.parent().altloc)

  xrs = processed_pdb_file.xray_structure()

  atom_numbers = flex.int()
  atom_symbols = flex.std_string()
  atom_charges = flex.std_string()
  for scatterer in xrs.scatterers():
    e,c = scatterer.element_and_charge_symbols()
    atom_numbers.append(scatterer.electron_count())
    atom_symbols.append(e)
    atom_charges.append(c)

  for i, (n, e, c, a) in enumerate(zip(atom_numbers,
                                       atom_symbols,
                                       atom_charges,
                                       atom_altlocs)):
    print(' "%2d" "%2s" "%-2s" "%1s" %s' % (n, e, c, a, atoms[i].charge))
  print()

  bond_proxies_simple, asu = geometry.get_covalent_bond_proxies(sites_cart =
    xrs.sites_cart())

  sites_cart = xrs.sites_cart()

  #bond_deltas = flex.double()

  #loop over bonds
  i_seqs = flex.int()
  j_seqs = flex.int()
  for proxy in bond_proxies_simple:
    #print proxy.i_seqs, proxy.distance_ideal
    i,j = proxy.i_seqs
    i_seqs.append(i)
    j_seqs.append(j)

    #site_1 = sites_cart[i]
    #site_2 = sites_cart[j]
    #atom1 = atoms[i]
    #atom2 = atoms[j]
    #print i,j,atom1.quote(), atom2.quote(),
    #print dir(atom1)
    #print " '%s' '%s' " % (atom1.parent().altloc, atom2.parent().altloc)
    #assert 0
    #dist_model = math.sqrt((site_1[0]-site_2[0])**2+(site_1[1]-site_2[1])**2+(site_1[2]-site_2[2])**2)
    #bond_deltas.append(proxy.distance_ideal-dist_model)
  #print list(flex.abs(bond_deltas))
  print('Bond i_seqs j_seqs')
  for i,j in zip(i_seqs, j_seqs):
    print(i,j)

def run(filename=None):
  should_have_imported_proper_methods()
  print("OK")
  should_return_atom_types()
  print("OK")
  exercise(filename)
  print("OK")

if (__name__ == "__main__"):
  run(*tuple(sys.argv[1:]))
