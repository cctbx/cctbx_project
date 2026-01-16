from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from cctbx.maptbx.bcr import qmap
from cctbx import maptbx
import mmtbx.model
import iotbx.map_manager

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

pdb_str = """
CRYST1   30.410   14.706   21.507  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       5.043   9.706  12.193  1.00  3.00           N
ATOM      2  CA  GLY A   1       5.000   9.301  10.742  1.00  3.00           C
ATOM      3  C   GLY A   1       6.037   8.234  10.510  1.00  3.00           C
ATOM      4  O   GLY A   1       6.529   7.615  11.472  1.00  3.00           O
ATOM      5  N   ASN A   2       6.396   8.017   9.246  1.00  3.00           N
ATOM      6  CA  ASN A   2       7.530   7.132   8.922  1.00  3.00           C
ATOM      7  C   ASN A   2       8.811   7.631   9.518  1.00  3.00           C
ATOM      8  O   ASN A   2       9.074   8.836   9.517  1.00  3.00           O
ATOM      9  CB  ASN A   2       7.706   6.975   7.432  1.00  3.00           C
ATOM     10  CG  ASN A   2       6.468   6.436   6.783  1.00  3.00           C
ATOM     11  OD1 ASN A   2       6.027   5.321   7.107  1.00  3.00           O
ATOM     12  ND2 ASN A   2       5.848   7.249   5.922  1.00  3.00           N
ATOM     13  N   ASN A   3       9.614   6.684   9.996  1.00  3.00           N
ATOM     14  CA  ASN A   3      10.859   6.998  10.680  1.00  3.00           C
ATOM     15  C   ASN A   3      12.097   6.426   9.986  1.00  3.00           C
ATOM     16  O   ASN A   3      12.180   5.213   9.739  1.00  3.00           O
ATOM     17  CB  ASN A   3      10.793   6.472  12.133  1.00  3.00           C
ATOM     18  CG  ASN A   3      12.046   6.833  12.952  1.00  3.00           C
ATOM     19  OD1 ASN A   3      12.350   8.019  13.163  1.00  3.00           O
ATOM     20  ND2 ASN A   3      12.781   5.809  13.397  1.00  3.00           N
ATOM     21  N   GLN A   4      13.047   7.322   9.689  1.00 50.00           N
ATOM     22  CA  GLN A   4      14.436   6.982   9.290  1.00 50.00           C
ATOM     23  C   GLN A   4      15.487   7.700  10.179  1.00 50.00           C
ATOM     24  O   GLN A   4      15.599   8.937  10.206  1.00 50.00           O
ATOM     25  CB  GLN A   4      14.708   7.242   7.802  1.00 50.00           C
ATOM     26  CG  GLN A   4      15.996   6.552   7.304  1.00 50.00           C
ATOM     27  CD  GLN A   4      16.556   7.138   6.002  1.00 50.00           C
ATOM     28  OE1 GLN A   4      16.796   8.362   5.901  1.00 50.00           O
ATOM     29  NE2 GLN A   4      16.802   6.255   5.000  1.00 50.00           N
ATOM     30  N   GLN A   5      16.206   6.915  10.962  1.00 50.00           N
ATOM     31  CA  GLN A   5      17.322   7.455  11.731  1.00 50.00           C
ATOM     32  C   GLN A   5      18.646   6.862  11.263  1.00 50.00           C
ATOM     33  O   GLN A   5      18.820   5.640  11.145  1.00 50.00           O
ATOM     34  CB  GLN A   5      17.108   7.277  13.238  1.00 50.00           C
ATOM     35  CG  GLN A   5      15.881   8.044  13.738  1.00 50.00           C
ATOM     36  CD  GLN A   5      15.396   7.508  15.045  1.00 50.00           C
ATOM     37  OE1 GLN A   5      14.826   6.419  15.093  1.00 50.00           O
ATOM     38  NE2 GLN A   5      15.601   8.281  16.130  1.00 50.00           N
ATOM     39  N   ASN A   6      19.566   7.758  10.947  1.00 50.00           N
ATOM     40  CA  ASN A   6      20.883   7.404  10.409  1.00 50.00           C
ATOM     41  C   ASN A   6      21.906   7.855  11.415  1.00 50.00           C
ATOM     42  O   ASN A   6      22.271   9.037  11.465  1.00 50.00           O
ATOM     43  CB  ASN A   6      21.117   8.110   9.084  1.00 50.00           C
ATOM     44  CG  ASN A   6      20.013   7.829   8.094  1.00 50.00           C
ATOM     45  OD1 ASN A   6      19.850   6.698   7.642  1.00 50.00           O
ATOM     46  ND2 ASN A   6      19.247   8.841   7.770  1.00 50.00           N
ATOM     47  N   TYR A   7      22.344   6.911  12.238  1.00 50.00           N
ATOM     48  CA  TYR A   7      23.211   7.238  13.390  1.00 50.00           C
ATOM     49  C   TYR A   7      24.655   7.425  12.976  1.00 50.00           C
ATOM     50  O   TYR A   7      25.093   6.905  11.946  1.00 50.00           O
ATOM     51  CB  TYR A   7      23.113   6.159  14.460  1.00 50.00           C
ATOM     52  CG  TYR A   7      21.717   6.023  14.993  1.00 50.00           C
ATOM     53  CD1 TYR A   7      20.823   5.115  14.418  1.00 50.00           C
ATOM     54  CD2 TYR A   7      21.262   6.850  16.011  1.00 50.00           C
ATOM     55  CE1 TYR A   7      19.532   5.000  14.887  1.00 50.00           C
ATOM     56  CE2 TYR A   7      19.956   6.743  16.507  1.00 50.00           C
ATOM     57  CZ  TYR A   7      19.099   5.823  15.922  1.00 50.00           C
ATOM     58  OH  TYR A   7      17.818   5.683  16.382  1.00 50.00           O
ATOM     59  OXT TYR A   7      25.410   8.093  13.703  1.00 50.00           O
TER
END
"""

pml_str="""

cmd.bg_color("white")

load ./vrm_example.pdb, ph
cmd.hide("everything","ph")
cmd.show("sticks"    ,"ph")
cmd.show("spheres"   ,"ph")
#color blue, ph

#
# ﻿unbond (resname D8U and name D), (name O4')
#  unbond (resname D8U and name D), (resi 46 and name OE2)
#

# MAIN
load ./vrm_example.map, map9, 1 , ccp4
isomesh mesh9, map9, 2.5, (all), 0, 1, 1.8
cmd.color("blue","mesh9")

#                                                 #     map(s)
#                                                 #
#                                                 #     ball&stic settings
set stick_radius,0.1
util.cbaw("pmod")
set_color cyan,[ 0.38, 0.78, 0.80]
set_color blue,[ 0.11, 0.25, 0.88]
set_color red,[ 1.00, 0.13, 0.13]
#color cyan, element h & mod6
#color grey, element c & mod6
#color red, all # for the final picture only
#color black, name ca | name n | name c | name O
set sphere_scale,0.075
#                                                  #    map settings
set mesh_width,0.3
set mesh_radius,0.001 #0.012
#                                                  #    perfomance
set direct,1.
set orthoscopic,on
set ray_trace_fog_start,0.5
util.performance(0)
util.ray_shadows('light')
cmd.space('cmyk')
set antialias,on
#                                                  #

set_view (     0.960335672,   -0.271858096,   -0.062032577,    -0.006946662,    0.199068844,   -0.979960561,     0.278759062,    0.941522539,    0.189284056,    -0.000007488,    0.000003186,  -51.567455292,    16.824468613,    7.543096542,   11.807725906,    30.981864929,   72.153221130,   20.000000000 )

rebuild
#
"""

def run(d_min=2, table="wk1995"):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.setup_scattering_dictionaries(
    scattering_table = table, d_min=d_min)
  with open("vrm_example.pdb","w") as fo:
    fo.write(model.model_as_pdb())
  with open("vrm_example.pml","w") as fo:
    fo.write(pml_str)

  cs = pdb_inp.crystal_symmetry()
  hierarchy = pdb_inp.construct_hierarchy()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = cs.unit_cell(),
    space_group_info = cs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.5
    )
  n_real = crystal_gridding.n_real()

  resolutions = flex.double(model.size(), -1)

  sel = model.selection(string="resseq 1:3")
  resolutions = resolutions.set_selected(sel, 1)
  sel = model.selection(string="resseq 4:7")
  resolutions = resolutions.set_selected(sel, 4)

  print(flex.min(resolutions), flex.max(resolutions))

  OmegaMap = qmap.compute(
    xray_structure = model.get_xray_structure(),
    n_real         = n_real,
    resolutions    = resolutions,
    use_exp_table  = False,
    debug          = False).map_data()
  mm = iotbx.map_manager.map_manager(
    map_data                   = OmegaMap,
    unit_cell_grid             = OmegaMap.all(),
    unit_cell_crystal_symmetry = cs,
    wrapping                   = True)
  mm.write_map("vrm_example.map")

if (__name__ == "__main__"):
  run()
