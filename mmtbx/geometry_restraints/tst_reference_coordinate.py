from __future__ import division

from libtbx.test_utils import approx_equal
from mmtbx.geometry_restraints import reference_coordinate
import iotbx.pdb
from cctbx.array_family import flex
from cctbx import adp_restraints # import dependency
from mmtbx.monomer_library import server, pdb_interpretation
from cStringIO import StringIO

def simple_pdb () :
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
CRYST1   36.670   40.710   66.290  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   THR A   2      -3.791  -8.769  29.092  1.00 24.15           N
ATOM      2  CA  THR A   2      -3.627  -7.675  28.090  1.00 25.97           C
ATOM      3  C   THR A   2      -2.202  -7.127  28.152  1.00 24.18           C
ATOM      4  O   THR A   2      -1.633  -6.984  29.233  1.00 24.71           O
ATOM      5  CB  THR A   2      -4.627  -6.527  28.357  1.00 26.50           C
ATOM      6  OG1 THR A   2      -5.961  -7.056  28.404  1.00 28.79           O
ATOM      7  CG2 THR A   2      -4.548  -5.486  27.255  1.00 27.05           C
ATOM      8  N   LYS A   3      -1.629  -6.832  26.988  1.00 24.44           N
ATOM      9  CA  LYS A   3      -0.266  -6.307  26.901  1.00 25.16           C
ATOM     10  C   LYS A   3      -0.196  -4.896  27.485  1.00 23.66           C
ATOM     11  O   LYS A   3      -1.094  -4.084  27.265  1.00 23.75           O
ATOM     12  CB  LYS A   3       0.199  -6.262  25.438  1.00 26.61           C
ATOM     13  CG ALYS A   3       0.312  -7.619  24.754  0.50 27.88           C
ATOM     14  CG BLYS A   3       0.201  -7.603  24.718  0.50 27.66           C
ATOM     15  CD ALYS A   3       1.436  -8.454  25.347  0.50 27.58           C
ATOM     16  CD BLYS A   3       1.205  -8.570  25.325  0.50 27.30           C
ATOM     17  CE ALYS A   3       1.585  -9.783  24.621  0.50 28.69           C
ATOM     18  CE BLYS A   3       1.213  -9.893  24.575  0.50 28.17           C
ATOM     19  NZ ALYS A   3       0.362 -10.624  24.732  0.50 28.63           N
ATOM     20  NZ BLYS A   3       2.149 -10.873  25.188  0.50 27.40           N
ATOM     21  N   LYS A   4       0.873  -4.612  28.225  1.00 22.24           N
ATOM     22  CA  LYS A   4       1.068  -3.295  28.826  1.00 21.81           C
ATOM     23  C   LYS A   4       2.337  -2.642  28.295  1.00 19.26           C
ATOM     24  O   LYS A   4       3.417  -3.243  28.310  1.00 18.66           O
ATOM     25  CB  LYS A   4       1.156  -3.398  30.354  1.00 23.29           C
ATOM     26  CG  LYS A   4      -0.170  -3.685  31.031  1.00 27.60           C
ATOM     27  CD  LYS A   4      -0.049  -3.681  32.551  1.00 32.16           C
ATOM     28  CE  LYS A   4       0.797  -4.842  33.052  1.00 33.04           C
ATOM     29  NZ  LYS A   4       0.827  -4.892  34.541  1.00 36.05           N
""")
  return pdb_in

def exercise_simple () :
  pdb_in = simple_pdb()
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  params = pdb_interpretation.master_params.extract()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    pdb_inp=pdb_in,
    log=StringIO())
  grm = processed_pdb_file.geometry_restraints_manager()
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  #pdb_hierarchy.atoms().reset_i_seq()
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  reference_coordinate_proxies = \
    reference_coordinate.build_proxies(
      sites_cart=sites_cart,
      pdb_hierarchy=pdb_hierarchy,
      c_alpha_only=True).reference_coordinate_proxies
  grads = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  residual = reference_coordinate.target_and_gradients(
               proxies=reference_coordinate_proxies,
               sites_cart=sites_cart,
               gradient_array=grads)
  assert approx_equal(residual, 0.0)

  assert grm.generic_restraints_manager.reference_coordinate_proxies is None

  grm.generic_restraints_manager.add_reference_restraints(
    sites_cart=sites_cart,
    pdb_hierarchy=pdb_hierarchy)

  assert grm.generic_restraints_manager.reference_coordinate_proxies \
    is not None
  assert len(grm.generic_restraints_manager.reference_coordinate_proxies) \
    == 29

  #test selection
  tst_iselection = flex.size_t()
  for atom in pdb_hierarchy.atoms():
    if atom.name == " CA " or atom.name == " N  ":
      tst_iselection.append(atom.i_seq)
  selection = flex.bool(len(sites_cart), tst_iselection)
  grm.generic_restraints_manager.add_reference_restraints(
    sites_cart=sites_cart,
    pdb_hierarchy=pdb_hierarchy,
    selection=selection)
  assert len(grm.generic_restraints_manager.reference_coordinate_proxies) \
    == 6

if (__name__ == "__main__") :
  exercise_simple()
  print "OK"
