from __future__ import division
import os
from libtbx import easy_run

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
        }

links = {
  'linking_test_F3S.pdb' : [0,3], # 3 bonds from F3S to SGs
  }

def run_and_test(cmd, pdb, i):
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert (result.return_code == 0)
  for line in result.stdout_lines :
    if ("Write PDB file" in line) :
      break
  else :
    raise RuntimeError("Missing expected log output")
  print "OK"
  # test .geo
  f=file(pdb.replace(".pdb", "_minimized.geo"), "rb")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line.find("Metal coordination restraints:")>-1:
      bonds = int(line.split()[-1])
      break
  assert bonds == links[pdb][i], "found %d bonds but expected %s!" % (
    bonds,
    links[pdb][i],
    )
  new_geo = pdb.replace(".pdb", "_minimized_%d.geo" % i)
  if (os.path.isfile(new_geo)) :
    os.remove(new_geo)
  os.rename(pdb.replace(".pdb", "_minimized.geo"), new_geo)
  print "OK"
  # test links
  if pdb in ["linking_test_LEU-CSY-VAL.pdb",
             ]:
    return
  number_of_links=0
  f=file(pdb.replace(".pdb", "_minimized.pdb"), "rb")
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

def run(only_i=None):
  try: only_i=int(only_i)
  except ValueError: only_i=None
  except TypeError: only_i=None
  cifs = ""
  for pdb in pdbs:
    f=file(pdb, "wb")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
  j=0
  for pdb in sorted(pdbs):
    #break
    if pdb.endswith(".cif"): continue
    if pdb.endswith(".params"): continue
    #if pdb.find('F3S')==-1: continue
    print 'pdb',pdb
    j+=1
    if only_i is not None and only_i!=j: continue
    for i in range(2):
      log_filename = "%s_%d.log" % (pdb, i)
      cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
      cmd += " link_all=%d link_carbohydrate=%d %s" % (i, i, cifs)
      cmd += " link_metal=%d" % i
      if pdb.replace(".pdb", ".params") in pdbs:
        cmd += " %s" % pdb.replace(".pdb", ".params")
      print "test number:",j
      print cmd
      run_and_test(cmd, pdb,i)


if __name__=="__main__":
  import sys
  run(*tuple(sys.argv[1:]))
