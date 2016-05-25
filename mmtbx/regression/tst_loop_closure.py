from __future__ import division
import iotbx.pdb
from mmtbx.building import loop_idealization

def exercise_ligand_after_chain():
  tst_pdb_1 = """\
ATOM      1  N   MET A   1    -159.943 -66.661-130.499  1.00 48.44           N
ATOM      2  CA  MET A   1    -160.787 -67.832-130.110  1.00 48.44           C
ATOM      3  C   MET A   1    -159.910 -69.027-129.841  1.00 48.44           C
ATOM      4  O   MET A   1    -158.711 -68.884-129.606  1.00 48.44           O
ATOM      5  CB  MET A   1    -161.759 -68.208-131.231  1.00 79.81           C
ATOM      6  CG  MET A   1    -163.209 -68.175-130.820  1.00 79.81           C
ATOM      7  SD  MET A   1    -163.833 -66.479-130.739  1.00 79.81           S
ATOM      8  CE  MET A   1    -163.171 -65.924-129.172  1.00 79.81           C
ATOM      9  N   ARG A   2    -160.518 -70.210-129.875  1.00 38.99           N
ATOM     10  CA  ARG A   2    -159.772 -71.438-129.666  1.00 38.99           C
ATOM     11  C   ARG A   2    -158.918 -71.422-128.395  1.00 38.99           C
ATOM     12  O   ARG A   2    -158.210 -72.391-128.104  1.00 38.99           O
ATOM     13  CB  ARG A   2    -158.865 -71.688-130.880  1.00 51.43           C
ATOM     14  CG  ARG A   2    -159.570 -72.239-132.112  1.00 51.43           C
ATOM     15  CD  ARG A   2    -159.999 -73.663-131.835  1.00 51.43           C
ATOM     16  NE  ARG A   2    -158.879 -74.455-131.325  1.00 51.43           N
ATOM     17  CZ  ARG A   2    -159.007 -75.488-130.497  1.00 51.43           C
ATOM     18  NH1 ARG A   2    -157.929 -76.147-130.090  1.00 51.43           N
ATOM     19  NH2 ARG A   2    -160.212 -75.855-130.066  1.00 51.43           N
ATOM     20  H   ARG A   2    -161.482 -70.262-130.050  1.00 51.43           H
ATOM     21  N   CYS A   3    -158.986 -70.340-127.630  1.00152.30           N
ATOM     22  CA  CYS A   3    -158.172 -70.252-126.431  1.00152.30           C
ATOM     23  C   CYS A   3    -158.912 -69.638-125.264  1.00152.30           C
ATOM     24  O   CYS A   3    -158.554 -69.870-124.109  1.00152.30           O
ATOM     25  CB  CYS A   3    -156.939 -69.403-126.697  1.00 59.00           C
ATOM     26  SG  CYS A   3    -157.266 -67.650-126.486  1.00 59.00           S
ATOM     27  H   CYS A   3    -159.604 -69.611-127.838  1.00 59.00           H
ATOM     28  N   ILE A   4    -159.926 -68.834-125.553  1.00 29.62           N
ATOM     29  CA  ILE A   4    -160.661 -68.204-124.470  1.00 29.62           C
ATOM     30  C   ILE A   4    -161.101 -69.218-123.409  1.00 29.62           C
ATOM     31  O   ILE A   4    -162.062 -69.976-123.578  1.00 29.62           O
ATOM     32  CB  ILE A   4    -161.841 -67.403-125.017  1.00 17.37           C
ATOM     33  CG1 ILE A   4    -161.294 -66.286-125.911  1.00 17.37           C
ATOM     34  CG2 ILE A   4    -162.683 -66.833-123.870  1.00 17.37           C
ATOM     35  CD1 ILE A   4    -160.294 -65.374-125.206  1.00 17.37           C
ATOM     36  H   ILE A   4    -160.176 -68.660-126.488  1.00 17.37           H
ATOM     37  N   GLY A   5    -160.346 -69.216-122.314  1.00 52.08           N
ATOM     38  CA  GLY A   5    -160.597 -70.120-121.222  1.00 52.08           C
ATOM     39  C   GLY A   5    -159.905 -71.424-121.536  1.00 52.08           C
ATOM     40  O   GLY A   5    -160.566 -72.368-121.974  1.00 52.08           O
ATOM     41  H   GLY A   5    -159.576 -68.614-122.236  1.00 37.04           H
ATOM   4278  N   SER B   1    -138.139-105.043-112.117  1.00 28.15           N
ATOM   4279  CA  SER B   1    -138.513-105.322-113.546  1.00 28.15           C
ATOM   4280  C   SER B   1    -139.470-104.237-113.982  1.00 28.15           C
ATOM   4281  O   SER B   1    -139.162-103.042-113.916  1.00 28.15           O
ATOM   4282  CB  SER B   1    -139.226-106.673-113.658  1.00 31.14           C
ATOM   4283  OG  SER B   1    -140.493-106.625-113.024  1.00 31.14           O
ATOM   4284  N   VAL B   2    -140.632-104.660-114.452  1.00 60.77           N
ATOM   4285  CA  VAL B   2    -141.655-103.705-114.795  1.00 60.77           C
ATOM   4286  C   VAL B   2    -142.205-103.400-113.397  1.00 60.77           C
ATOM   4287  O   VAL B   2    -143.408-103.192-113.192  1.00 60.77           O
ATOM   4288  CB  VAL B   2    -142.741-104.315-115.717  1.00108.23           C
ATOM   4289  CG1 VAL B   2    -142.440-103.962-117.172  1.00108.23           C
ATOM   4290  CG2 VAL B   2    -142.790-105.830-115.547  1.00108.23           C
ATOM   4291  H   VAL B   2    -140.797-105.604-114.625  1.00108.23           H
ATOM   4292  N   ALA B   3    -141.279-103.408-112.436  1.00 28.18           N
ATOM   4293  CA  ALA B   3    -141.548-103.125-111.037  1.00 28.18           C
ATOM   4294  C   ALA B   3    -140.547-102.051-110.611  1.00 28.18           C
ATOM   4295  O   ALA B   3    -140.863-101.159-109.825  1.00 28.18           O
ATOM   4296  CB  ALA B   3    -141.362-104.382-110.208  1.00 75.10           C
ATOM   4297  H   ALA B   3    -140.342-103.592-112.556  1.00 75.10           H
ATOM   4298  N   LEU B   4    -139.339-102.128-111.157  1.00 35.75           N
ATOM   4299  CA  LEU B   4    -138.314-101.150-110.823  1.00 35.75           C
ATOM   4300  C   LEU B   4    -137.667-100.486-112.046  1.00 35.75           C
ATOM   4301  O   LEU B   4    -138.011-100.795-113.198  1.00 35.75           O
ATOM   4302  CB  LEU B   4    -137.245-101.810-109.967  1.00 58.18           C
ATOM   4303  CG  LEU B   4    -137.744-102.589-108.750  1.00 58.18           C
ATOM   4304  CD1 LEU B   4    -138.674-101.712-107.940  1.00 58.18           C
ATOM   4305  CD2 LEU B   4    -138.442-103.866-109.197  1.00 58.18           C
ATOM   4306  H   LEU B   4    -139.132-102.867-111.766  1.00 58.18           H
ATOM   4307  N   VAL B   5    -136.718 -99.586-111.779  1.00 72.31           N
ATOM   4308  CA  VAL B   5    -136.022 -98.838-112.829  1.00 72.31           C
ATOM   4309  C   VAL B   5    -137.120 -98.370-113.765  1.00 72.31           C
ATOM   4310  O   VAL B   5    -137.069 -98.595-114.969  1.00 72.31           O
ATOM   4311  CB  VAL B   5    -135.032 -99.731-113.606  1.00188.85           C
ATOM   4312  CG1 VAL B   5    -134.206 -98.880-114.561  1.00188.85           C
ATOM   4313  CG2 VAL B   5    -134.125-100.474-112.636  1.00188.85           C
ATOM   4314  H   VAL B   5    -136.485 -99.417-110.841  1.00188.85           H
HETATM 9811  C1  NAG A 501    -145.629-131.141-126.187  1.00105.75           C
HETATM 9812  C2  NAG A 501    -145.718-131.799-127.560  1.00105.75           C
HETATM 9813  C3  NAG A 501    -147.101-132.414-127.738  1.00105.75           C
HETATM 9814  C4  NAG A 501    -147.345-133.422-126.615  1.00105.75           C
HETATM 9815  C5  NAG A 501    -147.160-132.734-125.244  1.00105.75           C
HETATM 9816  C6  NAG A 501    -147.264-133.704-124.079  1.00105.75           C
HETATM 9817  C7  NAG A 501    -144.518-131.062-129.492  1.00105.75           C
HETATM 9818  C8  NAG A 501    -144.964-131.298-130.926  1.00105.75           C
HETATM 9819  N2  NAG A 501    -145.464-130.823-128.595  1.00105.75           N
HETATM 9820  O3  NAG A 501    -147.184-133.061-129.000  1.00105.75           O
HETATM 9821  O4  NAG A 501    -148.659-133.960-126.729  1.00105.75           O
HETATM 9822  O5  NAG A 501    -145.851-132.118-125.159  1.00105.75           O
HETATM 9823  O6  NAG A 501    -148.610-133.860-123.650  1.00105.75           O
HETATM 9824  O7  NAG A 501    -143.325-131.114-129.200  1.00105.75           O
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=tst_pdb_1)
  pdb_h = pdb_inp.construct_hierarchy()
  loop_ideal_params = loop_idealization.master_phil.extract()
  loop_ideal_params.loop_idealization.enabled=True
  loop_ideal_params.loop_idealization.variant_search_level=1
  loop_ideal_params.loop_idealization.variant_number_cutoff=10
  loop_ideal_params.loop_idealization.number_of_ccd_trials=1
  loop_ideal_params.loop_idealization.minimize_whole=False
  loop_ideal = loop_idealization.loop_idealization(
      pdb_hierarchy=pdb_h,
      params = loop_ideal_params,
      secondary_structure_annotation=None,
      verbose=False)


def exercise():
  exercise_ligand_after_chain()
  print "OK"

if (__name__ == "__main__"):
  exercise()
