
from __future__ import division
from __future__ import print_function
import libtbx.phil
from cStringIO import StringIO

master_params = """
set_partial_occupancy = None
  .type = float
rmsd_min = 1.0
  .type = float
rigid_body_refine = False
  .type = bool
translation = 0.5
  .type = float
  .help = Specifies the radius for the optional translation search, in \
    Angstroms
translation_sampling = 10
  .type = int
  .optional = False
  .help = Number of distances within the translation search to sample.
torsion_search {
  include scope mmtbx.building.alternate_conformations.conformer_generation.torsion_search_params
}"""

def master_phil():
  return libtbx.phil.parse(master_params, process_includes=True)

def exercise():
  from mmtbx.building.alternate_conformations import conformer_generation
  from mmtbx.monomer_library import server
  import iotbx.pdb.hierarchy
  generate_inputs()
  params = master_phil().extract()
  mon_lib_srv = server.server()
  pdb_in = iotbx.pdb.hierarchy.input(file_name="ser_frag_single.pdb")
  hierarchy = pdb_in.hierarchy
  pdb_atoms = hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  sites_cart = pdb_atoms.extract_xyz()
  xrs = pdb_in.input.xray_structure_simple()
  models = []
  prev_res = next_res = next_next_res = None
  for chain in hierarchy.only_model().chains():
    residue_groups = chain.residue_groups()
    n_rg = len(residue_groups)
    for i_res, residue_group in enumerate(residue_groups):
      sites_orig = sites_cart.deep_copy()
      next_res = next_next_res = None
      if (i_res < (n_rg - 1)):
        next_res = residue_groups[i_res+1].atom_groups()[0]
      if (i_res < (n_rg - 2)):
        next_next_res = residue_groups[i_res+2].atom_groups()[0]
      atom_groups = residue_group.atom_groups()
      primary_conf = atom_groups[0]
      out = StringIO()
      confs = []
      for conf in conformer_generation.generate_single_residue_confs(
        atom_group=primary_conf,
        sites_cart=sites_cart.deep_copy(),
        mon_lib_srv=mon_lib_srv,
        params=params.torsion_search,
        prev_residue=prev_res,
        next_residue=next_res,
        next_next_residue=next_next_res,
        backrub=True,
        shear=False):
          conf.show_summary(out=out)
          confs.append(conf)
      prev_res = primary_conf
      if (confs is None):
        continue
      if (i_res == 1):
        assert ("""  SER A  99     20.0    None       t""" in out.getvalue())
      for conf in confs :
        sites_new = sites_cart.set_selected(conf.sites_selection,
          conf.sites_selected())
        pdb_atoms.set_xyz(sites_new)
        models.append(hierarchy.only_model().detached_copy())
      confs = []
      for conf in conformer_generation.generate_single_residue_confs(
        atom_group=primary_conf,
        sites_cart=sites_cart.deep_copy(),
        mon_lib_srv=mon_lib_srv,
        params=params.torsion_search,
        prev_residue=prev_res,
        next_residue=next_res,
        next_next_residue=next_next_res,
        backrub=True,
        shear=False):
          conf.show_summary(out=out)
          confs.append(conf)
      if (i_res == 1):
        print(len(confs))
  new_hierarchy = iotbx.pdb.hierarchy.root()
  for i_model, conf in enumerate(models):
    conf.id = str(i_model + 1)
    new_hierarchy.append_model(conf)
  open("ser_frag_naive_ensemble.pdb", "w").write(new_hierarchy.as_pdb_string())

def generate_inputs():
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(source_info=None, pdb_string="""\
REMARK derived from RT CypA structure (3k0n)
ATOM   1831  N  AALA A  98       5.400   4.097   9.945  0.63 11.08           N
ANISOU 1831  N  AALA A  98     1121   1491   1600    644    -58    486       N
ATOM   1832  N  BALA A  98       5.228   4.108  10.122  0.37 18.12           N
ANISOU 1832  N  BALA A  98     4482    783   1620  -1479    834   -410       N
ATOM   1833  CA AALA A  98       5.592   5.439  10.505  0.63 15.52           C
ANISOU 1833  CA AALA A  98     1431   1924   2544    946   1075    233       C
ATOM   1834  CA BALA A  98       5.523   5.512  10.398  0.37 17.11           C
ANISOU 1834  CA BALA A  98     3494   1013   1994  -1249  -1691    -89       C
ATOM   1835  C  AALA A  98       6.444   6.264   9.502  0.63 13.18           C
ANISOU 1835  C  AALA A  98     1271   2103   1635    252    567   -346       C
ATOM   1836  C  BALA A  98       6.278   6.121   9.224  0.37 11.74           C
ANISOU 1836  C  BALA A  98     1429    435   2596   -547   -978    803       C
ATOM   1837  O  AALA A  98       7.597   5.924   9.185  0.63  9.61           O
ANISOU 1837  O  AALA A  98     1179   1252   1220   -263   -271   -407       O
ATOM   1838  O  BALA A  98       7.194   5.512   8.675  0.37 10.05           O
ANISOU 1838  O  BALA A  98      982   1353   1482    -20   -562   -723       O
ATOM   1839  CB AALA A  98       6.272   5.307  11.875  0.63 24.37           C
ANISOU 1839  CB AALA A  98     3289   3773   2200   -572   -515    627       C
ATOM   1840  CB BALA A  98       6.403   5.648  11.644  0.37 13.06           C
ANISOU 1840  CB BALA A  98      810    880   3272    -73    748  -1092       C
ATOM   1847  H  AALA A  98       6.130   3.647   9.883  0.63 13.18           H
ATOM   1848  H  BALA A  98       5.806   3.571  10.465  0.37 21.63           H
ATOM   1849  HA AALA A  98       4.733   5.874  10.621  0.63 18.51           H
ATOM   1850  HA BALA A  98       4.699   6.005  10.537  0.37 20.41           H
ATOM   1851  HB2AALA A  98       5.692   4.768  12.435  0.63 29.13           H
ATOM   1852  HB2BALA A  98       7.161   5.054  11.529  0.37 15.55           H
ATOM   1853  HB3AALA A  98       7.104   4.828  11.740  0.63 29.13           H
ATOM   1854  HB3BALA A  98       6.727   6.562  11.671  0.37 15.55           H
ATOM      0  HB1AALA A  98       6.405   6.189  12.257  0.63 13.06           H
ATOM      0  HB1BALA A  98       6.587   6.586  11.810  0.37 13.06           H
ATOM   1869  N  ASER A  99       5.848   7.341   9.003  0.63 12.49           N
ANISOU 1869  N  ASER A  99     1356   1631   1760    944    173   -221       N
ATOM   1870  N  BSER A  99       5.897   7.337   8.855  0.37 15.41           N
ANISOU 1870  N  BSER A  99     3723    682   1450  -1448   -215    298       N
ATOM   1871  CA ASER A  99       6.394   8.125   7.899  0.63 13.11           C
ANISOU 1871  CA ASER A  99     2371   1161   1448    -17    736    285       C
ATOM   1872  CA BSER A  99       6.614   8.098   7.832  0.37 16.40           C
ANISOU 1872  CA BSER A  99     3570    146   2514    263  -1594   -223       C
ATOM   1873  C  ASER A  99       6.336   9.629   8.205  0.63 12.52           C
ANISOU 1873  C  ASER A  99     1739   1336   1681   -614   1240   -147       C
ATOM   1874  C  BSER A  99       6.423   9.596   8.097  0.37 19.14           C
ANISOU 1874  C  BSER A  99     1866   1643   3764   1274  -2024   -476       C
ATOM   1875  O  ASER A  99       5.485  10.103   8.943  0.63 12.68           O
ANISOU 1875  O  ASER A  99     1407   1874   1538   -306    723    -73       O
ATOM   1876  O  BSER A  99       5.483  10.013   8.762  0.37 17.56           O
ANISOU 1876  O  BSER A  99     3755    731   2185    732   -108      8       O
ATOM   1877  CB ASER A  99       5.606   7.828   6.617  0.63 13.57           C
ANISOU 1877  CB ASER A  99     1824   1577   1755    814    409   -259       C
ATOM   1878  CB BSER A  99       6.150   7.721   6.417  0.37 16.85           C
ANISOU 1878  CB BSER A  99     1411   3011   1981     33    685  -1259       C
ATOM   1879  OG ASER A  99       6.268   8.360   5.479  0.63 13.95           O
ANISOU 1879  OG ASER A  99     2763   1384   1151    -98    729    -86       O
ATOM   1880  OG BSER A  99       4.883   8.287   6.124  0.37 27.21           O
ANISOU 1880  OG BSER A  99     3140   3376   3823   1918   -918  -2027       O
ATOM   1881  H  ASER A  99       5.100   7.648   9.297  0.63 14.88           H
ATOM   1882  H  BSER A  99       5.218   7.751   9.183  0.37 18.37           H
ATOM   1883  HA ASER A  99       7.321   7.878   7.753  0.63 15.61           H
ATOM   1884  HA BSER A  99       7.561   7.900   7.898  0.37 19.56           H
ATOM   1885  HB2ASER A  99       5.523   6.867   6.513  0.63 16.17           H
ATOM   1886  HB2BSER A  99       6.798   8.051   5.775  0.37 20.11           H
ATOM   1887  HB3ASER A  99       4.726   8.231   6.687  0.63 16.17           H
ATOM   1888  HB3BSER A  99       6.084   6.755   6.355  0.37 20.11           H
ATOM   1889  HG ASER A  99       5.837   8.194   4.801  0.63 16.62           H
ATOM   1890  HG BSER A  99       4.615   8.722   6.766  0.37 32.53           H
ATOM   1891  N  AALA A 100       7.250  10.392   7.626  0.63 15.44           N
ANISOU 1891  N  AALA A 100     1295   2790   1782    395    525   -223       N
ATOM   1892  N  BALA A 100       7.321  10.404   7.558  0.37 10.82           N
ANISOU 1892  N  BALA A 100     1931    291   1887   -576     48   -278       N
ATOM   1893  CA AALA A 100       7.303  11.846   7.832  0.63 14.86           C
ANISOU 1893  CA AALA A 100     1854   1399   2393   -745    515   -174       C
ATOM   1894  CA BALA A 100       7.355  11.848   7.846  0.37 12.08           C
ANISOU 1894  CA BALA A 100     2649   1349    590   1069   -937   -333       C
ATOM   1895  C  AALA A 100       6.303  12.576   6.947  0.63 12.24           C
ANISOU 1895  C  AALA A 100     2461   1276    915   -128   -213   -380       C
ATOM   1896  C  BALA A 100       6.380  12.629   6.953  0.37 17.09           C
ANISOU 1896  C  BALA A 100      986   3806   1703    281   1020    529       C
ATOM   1897  O  AALA A 100       6.229  12.341   5.735  0.63 15.74           O
ANISOU 1897  O  AALA A 100     2994   1293   1692    436    108    -47       O
ATOM   1898  O  BALA A 100       6.383  12.466   5.727  0.37 14.24           O
ANISOU 1898  O  BALA A 100     1035   3111   1265     66     28   -581       O
ATOM   1899  CB AALA A 100       8.709  12.365   7.527  0.63 18.70           C
ANISOU 1899  CB AALA A 100     1387   3317   2403     40    724    108       C
ATOM   1900  CB BALA A 100       8.769  12.376   7.622  0.37 17.08           C
ANISOU 1900  CB BALA A 100     2824    713   2953  -1235   -857    141       C
ATOM   1907  H  AALA A 100       7.863  10.095   7.102  0.63 18.41           H
ATOM   1908  H  BALA A 100       7.934  10.146   7.012  0.37 12.86           H
ATOM   1909  HA AALA A 100       7.105  12.042   8.761  0.63 17.72           H
ATOM   1910  HA BALA A 100       7.121  11.992   8.776  0.37 14.37           H
ATOM   1911  HB2AALA A 100       8.991  12.026   6.663  0.63 22.33           H
ATOM   1912  HB2BALA A 100       9.090  12.052   6.766  0.37 20.38           H
ATOM   1913  HB3AALA A 100       8.692  13.335   7.513  0.63 22.33           H
ATOM   1914  HB3BALA A 100       8.741  13.346   7.615  0.37 20.38           H
ATOM      0  HB1AALA A 100       8.737  13.325   7.665  0.63 17.08           H
ATOM      0  HB1BALA A 100       8.793  13.327   7.811  0.37 17.08           H
END
""")
  xrs = pdb_in.input.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  f = open("ser_frag.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()
  for chain in pdb_in.hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      residue_group.remove_atom_group(atom_groups[1])
      atom_groups[0].altloc = " "
      for atom in residue_group.atoms():
        atom.occ = 1.0
  f = open("ser_frag_single.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()

if (__name__ == "__main__"):
  exercise()
  print("OK")
