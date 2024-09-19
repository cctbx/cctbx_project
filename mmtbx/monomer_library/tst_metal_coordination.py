from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run
from six.moves import range

pdbs = {
  'linking_test_F3S.pdb' : '''
ATOM      1 HG11 VAL S 187     -40.693  -9.791   9.764  0.83  6.10           H
ATOM      2 HG21 VAL S 187     -40.266 -10.791  11.945  1.00  7.89           H
ATOM      3 HG21 THR S 227     -39.621 -13.489  13.393  0.98  4.82           H
ATOM      4 HG23 THR S 227     -38.478 -13.500  12.303  0.90  4.74           H
ATOM      5  ND2 ASN S 229     -37.467 -10.074  16.225  1.00  5.64           N
ATOM      6  HB2 ASN S 229     -39.518 -11.317  15.512  1.00  7.73           H
ATOM      7 HD22 ASN S 229     -37.341 -10.485  15.480  0.91  5.69           H
ATOM      8  CB  CYS S 231     -39.137  -6.330  13.077  1.00  5.84           C
ATOM      9  SG  CYS S 231     -38.302  -7.647  14.022  1.00  5.36           S
ATOM     10  HB2 CYS S 231     -39.693  -6.723  12.387  1.00  8.29           H
ATOM     11  HB3 CYS S 231     -38.460  -5.757  12.678  1.00  6.29           H
ATOM     12  CZ  PHE S 236     -35.277  -5.448  15.177  1.00  6.18           C
ATOM     13  HE2 PHE S 236     -36.307  -5.299  13.454  1.00  7.47           H
ATOM     14  HZ  PHE S 236     -35.528  -6.321  15.380  0.92  7.31           H
ATOM     15  CE3 TRP S 241     -39.272  -5.809   9.179  1.00  5.95           C
ATOM     16  HE3 TRP S 241     -38.510  -6.063   9.649  1.00  6.21           H
ATOM     17  HZ3 TRP S 241     -39.926  -7.665   8.658  1.00  5.38           H
ATOM     18  HB2 PRO S 242     -35.950  -7.216   7.609  0.98  7.17           H
ATOM     19  HG2 PRO S 242     -34.665  -7.362   9.537  1.00  5.32           H
ATOM     20  HD2 PRO S 242     -36.252  -5.820   9.951  1.00  7.08           H
ATOM     21  CB  CYS S 249     -37.959 -11.517   7.758  1.00  5.36           C
ATOM     22  SG  CYS S 249     -37.935 -12.287   9.420  1.00  5.03           S
ATOM     23  HA  CYS S 249     -36.025 -10.764   7.557  1.00  6.71           H
ATOM     24  HB2 CYS S 249     -38.359 -10.636   7.819  1.00  7.35           H
ATOM     25  N   ILE S 250     -34.870 -12.949   7.821  1.00  4.98           N
ATOM     26  H   ILE S 250     -34.658 -12.357   8.401  0.82  3.62           H
ATOM     27 HG12 ILE S 250     -32.286 -12.355   9.066  1.00  5.88           H
ATOM     28 HG13 ILE S 250     -32.306 -13.743   9.849  1.00  7.20           H
ATOM     29  H   GLY S 251     -35.959 -13.893   9.709  1.00  7.07           H
ATOM     30  CB  CYS S 252     -32.706 -13.130  13.306  1.00  5.00           C
ATOM     31  SG  CYS S 252     -33.261 -11.513  13.953  1.00  5.02           S
ATOM     32  H   CYS S 252     -34.873 -13.554  12.121  1.00  4.62           H
ATOM     33  HB3 CYS S 252     -32.490 -13.045  12.358  1.00  4.19           H
TER
ATOM     34  HE2 LYS L 232     -32.566  -7.382  11.557  0.97  4.74           H
ATOM     35  NE2 GLN L 237     -33.155  -9.602   7.635  1.00  5.48           N
ATOM     36 HE22 GLN L 237     -33.535 -10.015   8.294  0.91  3.83           H
TER
HETATM   37  S1  F3S S1003     -37.373  -8.610  10.393  1.00  5.80           S
HETATM   38  S2  F3S S1003     -34.483 -11.028  10.412  1.00  5.12           S
HETATM   39  S3  F3S S1003     -37.067 -11.059  12.944  1.00  5.03           S
HETATM   40  S4  F3S S1003     -34.752  -8.309  13.065  1.00  5.35           S
HETATM   41 FE1  F3S S1003     -36.744 -10.733  10.683  1.00  4.94          Fe2+
HETATM   42 FE3  F3S S1003     -36.932  -8.794  12.568  1.00  5.27          Fe2+
HETATM   43 FE4  F3S S1003     -34.843 -10.514  12.587  1.00  4.74          Fe2+
TER
HETATM   44  O   HOH S1186     -31.997  -7.261  14.807  1.00  5.59           O
HETATM   45  H1  HOH S1186     -32.757  -7.470  14.513  0.97 11.34           H
HETATM   46  O   HOH S1306     -34.747  -8.995  16.422  1.00  6.68           O
HETATM   47  H2  HOH S1306     -34.624  -8.759  15.624  0.97 10.19           H
TER
HETATM   48  O   HOH L 974     -31.298 -10.022  11.614  1.00  5.46           O
HETATM   49  H1  HOH L 974     -31.498 -10.034  10.798  1.00 12.35           H
HETATM   50  H2  HOH L 974     -32.009 -10.160  12.036  1.00 11.06           H
END
''',
  'linking_test_SF4.pdb' : '''
ATOM      1  CB  CYS A  43       3.194   5.758   4.327  1.00  3.39           C
ATOM      2  SG  CYS A  43       1.682   6.600   4.910  1.00  3.32           S
ATOM      3  CB  CYS A  46       6.976   9.042   5.345  1.00  3.45           C
ATOM      4  SG  CYS A  46       6.736  10.109   6.828  1.00  3.39           S
ATOM      5  CB  CYS A  61       2.826   5.722  11.129  1.00  4.19           C
ATOM      6  SG  CYS A  61       3.191   7.503  11.225  1.00  3.86           S
ATOM      7  CD1 LEU A  63       3.399  11.917  12.377  1.00  4.96           C
ATOM      8  N   CYS A  75       0.009  11.484   5.007  1.00  3.12           N
ATOM      9  CB  CYS A  75      -0.759  11.967   7.348  1.00  3.25           C
ATOM     10  SG  CYS A  75       0.712  12.086   8.405  1.00  3.33           S
TER
HETATM   11  S1  SF4 A  84       3.034  10.140   5.716  1.00  3.17           S
HETATM   12  S2  SF4 A  84       4.314   7.265   7.546  1.00  3.24           S
HETATM   13  S3  SF4 A  84       0.963   8.387   8.105  1.00  3.23           S
HETATM   14  S4  SF4 A  84       3.811  10.428   9.197  1.00  3.34           S
HETATM   15 FE1  SF4 A  84       2.436   8.034   6.473  1.00  3.03          Fe
HETATM   16 FE2  SF4 A  84       4.585   9.549   7.307  1.00  3.10          Fe
HETATM   17 FE3  SF4 A  84       3.040   8.257   9.103  1.00  3.11          Fe
HETATM   18 FE4  SF4 A  84       2.046  10.381   7.746  1.00  3.04          Fe
''',
  'linking_test_SF4_FE.pdb' : '''
CRYST1  366.325   97.352  133.276  90.00 108.37  90.00 C 1 2 1
SCALE1      0.002730  0.000000  0.000907        0.00000
SCALE2      0.000000  0.010272  0.000000        0.00000
SCALE3      0.000000  0.000000  0.007906        0.00000
ATOM      1  CB  CYS A 197      21.884  17.895  41.629  1.00 18.47           C
ATOM      2  SG  CYS A 197      23.219  18.455  40.509  1.00 22.44           S
TER       3      CYS A 197
ATOM      3  CA  CYS G  11      18.621  25.073  31.833  1.00 23.66           C
ATOM      4  CB  CYS G  11      18.380  24.228  33.075  1.00 23.77           C
ATOM      5  SG  CYS G  11      19.849  24.012  34.110  1.00 27.62           S
ATOM      6  N   CYS G  13      22.408  22.968  32.362  1.00 19.54           N
ATOM      7  CA  CYS G  13      22.966  21.958  33.256  1.00 20.10           C
ATOM      8  C   CYS G  13      24.416  21.774  32.862  1.00 24.21           C
ATOM      9  CB  CYS G  13      22.188  20.647  33.204  1.00 20.72           C
ATOM     10  SG  CYS G  13      22.629  19.484  34.522  1.00 24.93           S
ATOM     26  N   CYS G  46      18.764  25.237  39.142  1.00 28.00           N
ATOM     27  CA  CYS G  46      18.564  24.859  40.542  1.00 27.69           C
ATOM     28  CB  CYS G  46      18.593  23.338  40.716  1.00 27.39           C
ATOM     29  SG  CYS G  46      20.152  22.551  40.217  1.00 30.82           S
ATOM     34  N   CYS G  71      15.783  19.220  34.124  1.00 27.22           N
ATOM     35  CA  CYS G  71      15.862  17.804  34.451  1.00 27.48           C
ATOM     36  CB  CYS G  71      16.848  17.607  35.600  1.00 28.10           C
ATOM     37  SG  CYS G  71      16.505  18.649  37.042  1.00 32.48           S
HETATM   52 FE    FE G 705      19.937  20.352  40.205  0.50 22.13          Fe
HETATM   53  S1  SF4 G 706      21.592  22.331  36.609  1.00 24.36           S
HETATM   54  S2  SF4 G 706      17.978  22.365  37.279  1.00 22.14           S
HETATM   55  S3  SF4 G 706      19.266  20.138  34.768  1.00 23.53           S
HETATM   56  S4  SF4 G 706      20.047  19.377  38.173  1.00 21.03           S
HETATM   57 FE1  SF4 G 706      18.287  20.136  36.841  1.00 22.20          Fe
HETATM   58 FE2  SF4 G 706      21.110  20.098  36.275  1.00 26.08          Fe
HETATM   59 FE3  SF4 G 706      20.045  21.717  38.135  1.00 23.49          Fe
HETATM   60 FE4  SF4 G 706      19.481  22.300  35.554  1.00 23.79          Fe
''',
        }

links = {
  'linking_test_F3S.pdb' : [0,3], # 3 bonds from F3S to SGs
  'linking_test_SF4.pdb' : [0,4],
  'linking_test_SF4_FE.pdb' : [0,5],
  }

def run_and_test(cmd, pdb, i):
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert (result.return_code == 0)
  for line in result.stdout_lines :
    if ("Write PDB file" in line):
      break
  else :
    raise RuntimeError("Missing expected log output")
  print("OK")
  # test .geo
  f=open(pdb.replace(".pdb", "_minimized.geo"), "r")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line.find("Bond | Metal coordination | restraints")>-1:
      bonds = int(line.split()[-1])
      break
  else:
    bonds = 0
  assert bonds == links[pdb][i], "found %d bonds but expected %s!" % (
    bonds,
    links[pdb][i],
    )
  new_geo = pdb.replace(".pdb", "_minimized_%d.geo" % i)
  if (os.path.isfile(new_geo)):
    os.remove(new_geo)
  os.rename(pdb.replace(".pdb", "_minimized.geo"), new_geo)
  print("OK")
  # test links
  if pdb in ["linking_test_LEU-CSY-VAL.pdb",
             ]:
    return
  number_of_links=0
  f=open(pdb.replace(".pdb", "_minimized.pdb"), "r")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line.find("LINK")>-1:
      number_of_links+=1
  if i==0:
    expected = 0
  else:
    expected = links[pdb][i]-links[pdb][0]
    if pdb in ["linking_test_Mg_EDT.pdb",
               "linking_test_Mg_HOH.pdb",
              ]:
      expected += 6
    elif pdb in ["linking_test_CD_GHE_A_B.pdb",
                ]:
      expected += 4
  if pdb in ["linking_test_HEM_TYR.pdb"]:
    expected += 1
  assert number_of_links == expected, "found %d LINK but expected %s!" % (
    number_of_links,
    expected,
    )
  if 0:
    cmd = "phenix.start_coot --no-guano --pdb %s" % pdb.replace(".pdb",
                                                                "_minimized.pdb"
                                                                )
    os.system(cmd)

def ideal_and_test(cmd, pdb, i):
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert (result.return_code == 0)
  for line in result.stdout_lines :
    if ("Write PDB file" in line):
      break
  else :
    raise RuntimeError("Missing expected log output")
  print("OK")
  for line in result.stdout_lines :
    if line.find('SF4 Regularisation')>-1: break
    if line.find('F3S Regularisation')>-1: break
  else:
    assert 0
  print('OK')

def run(only_i=None):
  try: only_i=int(only_i)
  except ValueError: only_i=None
  except TypeError: only_i=None
  cifs = ""
  for pdb in pdbs:
    f=open(pdb, "w")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
  j=0
  for pdb in sorted(pdbs):
    #break
    if pdb.endswith(".cif"): continue
    if pdb.endswith(".params"): continue
    #if pdb.find('F3S')==-1: continue
    print('pdb',pdb)
    j+=1
    if only_i is not None and only_i!=j: continue
    for i in range(2):
      log_filename = "%s_%d.log" % (pdb, i)
      cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
      cmd += " link_all=%d link_carbohydrate=%d %s" % (i, i, cifs)
      cmd += " link_metal=%d" % i
      if pdb.replace(".pdb", ".params") in pdbs:
        cmd += " %s" % pdb.replace(".pdb", ".params")
      print("test number:", j)
      print(cmd)
      run_and_test(cmd, pdb, i)

      cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
      cmd += ' superpose_ideal_ligand=all'
      ideal_and_test(cmd, pdb, i)


if __name__=="__main__":
  import sys
  run(*tuple(sys.argv[1:]))
