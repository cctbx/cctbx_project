from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.utils
from scitbx.matrix import rotate_point_around_axis

pdb_1 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  1.00 10.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  1.00 10.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  1.00 10.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  1.00 10.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  1.00 10.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  1.00 10.00           N
TER
END
"""

pdb_2 = """\
CRYST1   12.298   14.179   12.287  90.00  90.00  90.00 P 1
ATOM    128  N   THR A  17       7.298   6.978   5.045  1.00 14.30           N
ATOM    129  CA  THR A  17       5.897   7.356   5.000  1.00 15.23           C
ATOM    130  C   THR A  17       5.803   8.761   5.487  1.00 12.29           C
ATOM    131  O   THR A  17       6.614   9.179   6.304  1.00 14.53           O
ATOM    132  CB  THR A  17       5.026   6.450   5.930  1.00 17.46           C
ATOM    133  OG1 THR A  17       5.485   6.544   7.287  1.00 19.58           O
ATOM    134  CG2 THR A  17       5.000   5.000   5.469  1.00 20.72           C
TER
END
"""

pdb_3="""\
ATOM    477  N   ARG A 160       0.631  -7.517 -35.759  1.00  0.00           N  
ATOM    478  CA  ARG A 160      -0.057  -7.335 -34.456  1.00  0.00           C  
ATOM    479  C   ARG A 160       0.314  -8.359 -33.357  1.00  0.00           C  
ATOM    480  O   ARG A 160       0.562  -7.977 -32.214  1.00  0.00           O  
ATOM    481  CB  ARG A 160      -1.600  -7.314 -34.675  1.00  0.00           C  
ATOM    482  CG  ARG A 160      -2.278  -5.934 -34.545  1.00  0.00           C  
ATOM    483  CD  ARG A 160      -3.243  -5.543 -35.683  1.00  0.00           C  
ATOM    484  NE  ARG A 160      -4.534  -6.286 -35.680  1.00  0.00           N  
ATOM    485  CZ  ARG A 160      -5.611  -5.955 -34.966  1.00  0.00           C  
ATOM    486  NH1 ARG A 160      -5.632  -5.023 -34.038  1.00  0.00           N  
ATOM    487  NH2 ARG A 160      -6.718  -6.601 -35.202  1.00  0.00           N  
"""

pdb_4 = """
CRYST1    5.573    3.411    5.179  90.00  90.00  90.00 P 1
ATOM    520  N   PHE A 162       0.000   1.817   3.119  1.00  0.01           N
ATOM    521  CA  PHE A 162       1.436   2.203   3.251  1.00  0.01           C
ATOM    522  C   PHE A 162       1.724   3.411   4.222  1.00  0.01           C
ATOM    523  O   PHE A 162       2.490   3.281   5.179  1.00  0.01           O
ATOM    524  CB  PHE A 162       2.039   2.358   1.819  1.00  0.01           C
ATOM    525  CG  PHE A 162       3.328   1.565   1.499  1.00  0.01           C
ATOM    526  CD1 PHE A 162       4.400   1.486   2.396  1.00  0.01           C
ATOM    527  CD2 PHE A 162       3.403   0.853   0.293  1.00  0.01           C
ATOM    528  CE1 PHE A 162       5.516   0.708   2.099  1.00  0.01           C
ATOM    529  CE2 PHE A 162       4.517   0.075   0.000  1.00  0.01           C
ATOM    530  CZ  PHE A 162       5.573   0.000   0.902  1.00  0.01           C
TER
END
"""


  
def torsion_search_nested(residue, clusters, rotamer_eval, states):
  n_angles = len(clusters)
  print n_angles
  rid = rotamer_eval.evaluate_residue(residue = residue)
  if(n_angles == 2):
    r1 = [0,0]
    r2 = [360,360]
  elif(n_angles == 3):
    r1 = [0,0,0]
    r2 = [360,360,360]
  elif(n_angles == 4):
    r1 = [0,0,0,0]
    r2 = [360,360,360,360]
  else: return
  nested_loop = flex.nested_loop(begin=r1, end=r2, open_range=False)
  xyz_moved_dc = residue.atoms().extract_xyz().deep_copy()
 
  nested_loop = []
  # 2
  #for a1 in range(-60, 61, 1):
  #  #for a2 in range(-90, 95, 5):
  #  for a2 in range(-90, 91, 1):
  #    nested_loop.append([a1,a2])
  # 4
  for a1 in range(-90, 91, 10):
    for a2 in range(-90, 91, 10):
      for a3 in range(-90, 91, 10):
        for a4 in range(-90, 91, 10):
          nested_loop.append([a1,a2,a3,a4])
  print "nested_loop", len(nested_loop)
  good = []
  good_sites = []
  for ia, angles in enumerate(nested_loop):
    #print ia
    xyz_moved = xyz_moved_dc.deep_copy()
    for i, angle in enumerate(angles):
      cl = clusters[i]
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = xyz_moved[cl.axis[0]],
          axis_point_2 = xyz_moved[cl.axis[1]],
          point        = xyz_moved[atom],
          angle        = angle, deg=True)
        xyz_moved[atom] = new_xyz
    residue.atoms().set_xyz(xyz_moved)
    if(rotamer_eval.evaluate_residue(residue = residue) ==rid):#!= "OUTLIER"):
      states.add(sites_cart=xyz_moved)
      #print "%4d %4d"%(angles[0], angles[1])
      good.append(angles)
  #
  x1 = flex.double()
  x2 = flex.double()
  x3 = flex.double()
  x4 = flex.double()
  for a in good:
    x1.append(a[0])
    x2.append(a[1])
    x3.append(a[2])
    x4.append(a[3])
  print flex.min(x1), flex.max(x1), "<>", flex.min(x2), flex.max(x2),\
        flex.min(x3), flex.max(x3), "<>", flex.min(x4), flex.max(x4), x1.size()
  #
  return states                                                   

def exercise():
  # XXX
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_1)
  pdb_inp.write_pdb_file(file_name = "answer.pdb")
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  mon_lib_srv = monomer_library.server.server()
  residue = pdb_hierarchy.only_residue()
  clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = residue,
    mon_lib_srv     = mon_lib_srv,
    backbone_sample = False,
    torsion_search_sidechain_start = 0,
    torsion_search_sidechain_stop  = 364,
    torsion_search_sidechain_step  = 4).clusters
  rotamer_eval = RotamerEval()
  ri = mmtbx.refinement.real_space.fit_residue.get_rotamer_iterator(
    mon_lib_srv = mon_lib_srv,
    residue     = residue)
  cntr = 0
  for rotamer, rotamer_sites_cart in ri:
    if 1:#(cntr == 0):
      residue.atoms().set_xyz(rotamer_sites_cart)
      xrs= xrs.replace_sites_cart(rotamer_sites_cart)
      states = mmtbx.utils.states(xray_structure=xrs, pdb_hierarchy=pdb_hierarchy)
    #print rotamer.id 
    #residue.atoms().set_xyz(rotamer_sites_cart)
    #states.add(sites_cart=residue.atoms().extract_xyz().deep_copy())
    states = torsion_search_nested( #torsion_search(
    #states = torsion_search(
      residue      = residue, 
      clusters     = clusters, 
      rotamer_eval = rotamer_eval, 
      states       = states)
    cntr+=1
    states.write(file_name="all_%d.pdb"%cntr)


if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "Time: %6.4f"%(time.time()-t0)
