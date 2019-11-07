from __future__ import absolute_import, division, print_function
from mmtbx.tls import tools
import time
import iotbx.pdb
import mmtbx.tls.tools
from scitbx import matrix
from mmtbx.tls import analysis
import math
from scitbx.array_family import flex
from six.moves import range

pdb_str= """
REMARK   3
REMARK   3   TLS GROUP : 6
REMARK   3    NUMBER OF COMPONENTS GROUP : 1
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI
REMARK   3    RESIDUE RANGE :   A    65        A    77
REMARK   3    ORIGIN FOR THE GROUP (A):   5.7798  -4.9101  20.8090
REMARK   3    T TENSOR
REMARK   3      T11:   0.1766 T22:   0.1958
REMARK   3      T33:   0.1284 T12:  -0.0119
REMARK   3      T13:  -0.0231 T23:   0.0335
REMARK   3    L TENSOR
REMARK   3      L11:   0.2047 L22:   0.2291
REMARK   3      L33:   0.4446 L12:   0.0057
REMARK   3      L13:   0.2320 L23:   0.2028
REMARK   3    S TENSOR
REMARK   3      S11:   0.0460 S12:   0.0164 S13:  -0.0519
REMARK   3      S21:   0.0242 S22:   0.0228 S23:  -0.0546
REMARK   3      S31:   0.0280 S32:   0.0569 S33:  -0.0688
REMARK   3
CRYST1   70.190   80.460  110.890  90.00  90.00  90.00 P 21 21 21
SCALE1      0.014247  0.000000  0.000000        0.00000
SCALE2      0.000000  0.012429  0.000000        0.00000
SCALE3      0.000000  0.000000  0.009018        0.00000
ATOM      1  N   ASP A  65      14.147  -9.221  30.281  1.00 41.09           N
ANISOU    1  N   ASP A  65     5329   5873   4409   -633  -1875    965       N
ATOM      2  CA  ASP A  65      15.252  -9.249  29.320  1.00 42.10           C
ANISOU    2  CA  ASP A  65     5520   5848   4628   -454  -1645    509       C
ATOM      3  C   ASP A  65      15.235 -10.537  28.501  1.00 40.51           C
ANISOU    3  C   ASP A  65     5402   5699   4291   -472  -2125    624       C
ATOM      4  O   ASP A  65      14.791 -11.574  28.982  1.00 39.28           O
ANISOU    4  O   ASP A  65     5521   5937   3467   -466  -2793    943       O
ATOM      5  CB  ASP A  65      16.599  -9.117  30.052  1.00 44.12           C
ANISOU    5  CB  ASP A  65     5578   6377   4807   -323  -1759    465       C
ATOM      6  CG  ASP A  65      16.805  -7.734  30.675  1.00 47.85           C
ANISOU    6  CG  ASP A  65     6211   6644   5326   -533  -1632    224       C
ATOM      7  OD1 ASP A  65      16.037  -6.794  30.359  1.00 51.97           O
ANISOU    7  OD1 ASP A  65     6898   6544   6303    -15   -411    599       O
ATOM      8  OD2 ASP A  65      17.748  -7.584  31.482  1.00 51.60           O
ANISOU    8  OD2 ASP A  65     6250   7759   5593   -575  -1695   -226       O
ATOM      9  N   GLY A  66      15.695 -10.453  27.254  1.00 43.04           N
ANISOU    9  N   GLY A  66     5541   6021   4789   -564  -1568    173       N
ATOM     10  CA  GLY A  66      15.771 -11.614  26.358  1.00 42.81           C
ANISOU   10  CA  GLY A  66     5506   6041   4718   -265  -1800    137       C
ATOM     11  C   GLY A  66      14.469 -11.953  25.645  1.00 42.08           C
ANISOU   11  C   GLY A  66     5383   5987   4617   -209  -1639     57       C
ATOM     12  O   GLY A  66      14.459 -12.761  24.719  1.00 44.05           O
ANISOU   12  O   GLY A  66     5670   6379   4687   -247  -2191   -161       O
ATOM     13  N   ALA A  67      13.374 -11.318  26.052  1.00 39.31           N
ANISOU   13  N   ALA A  67     5011   5498   4427   -513  -1743    269       N
ATOM     14  CA  ALA A  67      12.035 -11.711  25.610  1.00 34.72           C
ANISOU   14  CA  ALA A  67     4873   4883   3435   -220  -1759    605       C
ATOM     15  C   ALA A  67      11.747 -11.272  24.182  1.00 34.65           C
ANISOU   15  C   ALA A  67     4719   5070   3375   -468  -1676    626       C
ATOM     16  O   ALA A  67      12.465 -10.435  23.611  1.00 33.02           O
ANISOU   16  O   ALA A  67     4398   5086   3060   -343  -2075    876       O
ATOM     17  CB  ALA A  67      10.980 -11.142  26.553  1.00 32.29           C
ANISOU   17  CB  ALA A  67     4476   4659   3132   -444  -1922    573       C
ATOM     18  N   ILE A  68      10.695 -11.864  23.612  1.00 33.06           N
ANISOU   18  N   ILE A  68     4600   4702   3256   -324  -1651    731       N
ATOM     19  CA  ILE A  68      10.159 -11.446  22.326  1.00 31.69           C
ANISOU   19  CA  ILE A  68     4462   4593   2985   -258  -1327    634       C
ATOM     20  C   ILE A  68       8.993 -10.516  22.600  1.00 28.77           C
ANISOU   20  C   ILE A  68     4017   4426   2485   -420  -1554    771       C
ATOM     21  O   ILE A  68       8.102 -10.863  23.367  1.00 29.09           O
ANISOU   21  O   ILE A  68     3969   4498   2585   -247  -1612    861       O
ATOM     22  CB  ILE A  68       9.711 -12.658  21.487  1.00 33.33           C
ANISOU   22  CB  ILE A  68     4545   4700   3419   -415  -1435    543       C
ATOM     23  CG1 ILE A  68      10.950 -13.475  21.091  1.00 33.58           C
ANISOU   23  CG1 ILE A  68     4108   5209   3440   -686  -2083    -48       C
ATOM     24  CG2 ILE A  68       8.940 -12.220  20.244  1.00 33.93           C
ANISOU   24  CG2 ILE A  68     4626   4880   3386   -228  -1345    557       C
ATOM     25  CD1 ILE A  68      10.723 -14.965  21.116  1.00 35.38           C
ANISOU   25  CD1 ILE A  68     4273   5208   3961   -444  -2399   -116       C
ATOM     26  N   LEU A  69       9.020  -9.336  21.977  1.00 27.70           N
ANISOU   26  N   LEU A  69     3666   4241   2616   -537  -1183    594       N
ATOM     27  CA  LEU A  69       7.993  -8.311  22.153  1.00 25.86           C
ANISOU   27  CA  LEU A  69     3574   4169   2081   -532  -1071    868       C
ATOM     28  C   LEU A  69       7.413  -8.016  20.799  1.00 25.13           C
ANISOU   28  C   LEU A  69     3441   4069   2036   -439  -1021    734       C
ATOM     29  O   LEU A  69       8.156  -7.877  19.831  1.00 26.34           O
ANISOU   29  O   LEU A  69     3552   4318   2137   -465   -987   1038       O
ATOM     30  CB  LEU A  69       8.613  -7.021  22.657  1.00 27.46           C
ANISOU   30  CB  LEU A  69     3857   4296   2280   -543   -940    640       C
ATOM     31  CG  LEU A  69       7.850  -6.047  23.563  1.00 28.66           C
ANISOU   31  CG  LEU A  69     4059   4448   2381   -330   -907    634       C
ATOM     32  CD1 LEU A  69       8.040  -4.642  23.035  1.00 30.10           C
ANISOU   32  CD1 LEU A  69     4062   4597   2774   -143   -556    890       C
ATOM     33  CD2 LEU A  69       6.379  -6.345  23.802  1.00 27.33           C
ANISOU   33  CD2 LEU A  69     3995   4289   2097   -164  -1171    609       C
ATOM     34  N   ILE A  70       6.097  -7.893  20.738  1.00 24.58           N
ANISOU   34  N   ILE A  70     3388   4002   1947   -615   -905    886       N
ATOM     35  CA  ILE A  70       5.415  -7.521  19.506  1.00 22.86           C
ANISOU   35  CA  ILE A  70     3182   3660   1841   -534   -756    774       C
ATOM     36  C   ILE A  70       4.678  -6.222  19.803  1.00 22.24           C
ANISOU   36  C   ILE A  70     3080   3605   1762   -608   -710    693       C
ATOM     37  O   ILE A  70       3.979  -6.128  20.805  1.00 23.02           O
ANISOU   37  O   ILE A  70     3306   3725   1714   -573   -615    824       O
ATOM     38  CB  ILE A  70       4.422  -8.609  19.076  1.00 22.84           C
ANISOU   38  CB  ILE A  70     3207   3610   1862   -508   -784    792       C
ATOM     39  CG1 ILE A  70       5.161  -9.951  18.912  1.00 23.90           C
ANISOU   39  CG1 ILE A  70     3382   3677   2020   -472   -958    841       C
ATOM     40  CG2 ILE A  70       3.706  -8.190  17.801  1.00 22.46           C
ANISOU   40  CG2 ILE A  70     3286   3414   1832   -387   -736    764       C
ATOM     41  CD1 ILE A  70       4.279 -11.135  18.629  1.00 24.53           C
ANISOU   41  CD1 ILE A  70     3665   3599   2056   -465  -1065    836       C
ATOM     42  N   PHE A  71       4.860  -5.208  18.973  1.00 21.59           N
ANISOU   42  N   PHE A  71     3035   3501   1667   -533   -584    589       N
ATOM     43  CA  PHE A  71       4.050  -3.990  19.088  1.00 21.28           C
ANISOU   43  CA  PHE A  71     2952   3507   1624   -540   -545    564       C
ATOM     44  C   PHE A  71       2.699  -4.237  18.451  1.00 21.01           C
ANISOU   44  C   PHE A  71     2922   3406   1654   -432   -574    613       C
ATOM     45  O   PHE A  71       2.629  -4.859  17.388  1.00 20.73           O
ANISOU   45  O   PHE A  71     2910   3342   1624   -180   -454    644       O
ATOM     46  CB  PHE A  71       4.686  -2.798  18.386  1.00 21.03           C
ANISOU   46  CB  PHE A  71     2929   3405   1654   -406   -537    615       C
ATOM     47  CG  PHE A  71       6.072  -2.463  18.847  1.00 21.48           C
ANISOU   47  CG  PHE A  71     2935   3594   1629   -484   -470    627       C
ATOM     48  CD1 PHE A  71       6.392  -2.412  20.202  1.00 21.38           C
ANISOU   48  CD1 PHE A  71     2856   3630   1635   -537   -462    666       C
ATOM     49  CD2 PHE A  71       7.055  -2.167  17.930  1.00 21.64           C
ANISOU   49  CD2 PHE A  71     2901   3553   1766   -462   -479    813       C
ATOM     50  CE1 PHE A  71       7.668  -2.084  20.609  1.00 21.79           C
ANISOU   50  CE1 PHE A  71     2851   3665   1763   -604   -328    661       C
ATOM     51  CE2 PHE A  71       8.334  -1.823  18.329  1.00 22.27           C
ANISOU   51  CE2 PHE A  71     2976   3735   1750   -509   -523    728       C
ATOM     52  CZ  PHE A  71       8.647  -1.798  19.675  1.00 21.40           C
ANISOU   52  CZ  PHE A  71     2754   3687   1689   -412   -381    695       C
ATOM     53  N   SER A  72       1.633  -3.709  19.049  1.00 20.51           N
ANISOU   53  N   SER A  72     2877   3406   1506   -663   -396    573       N
ATOM     54  CA  SER A  72       0.286  -3.952  18.532  1.00 20.43           C
ANISOU   54  CA  SER A  72     2885   3288   1590   -494   -528    655       C
ATOM     55  C   SER A  72       0.022  -3.221  17.213  1.00 20.34           C
ANISOU   55  C   SER A  72     2949   3214   1563   -483   -511    618       C
ATOM     56  O   SER A  72       0.781  -2.352  16.810  1.00 19.70           O
ANISOU   56  O   SER A  72     2751   3296   1435   -422   -264    418       O
ATOM     57  CB  SER A  72      -0.781  -3.534  19.540  1.00 20.57           C
ANISOU   57  CB  SER A  72     2835   3406   1571   -642   -526    686       C
ATOM     58  OG  SER A  72      -0.997  -2.142  19.497  1.00 20.81           O
ANISOU   58  OG  SER A  72     2841   3421   1642   -600   -568    401       O
ATOM     59  N   ALA A  73      -1.106  -3.558  16.591  1.00 20.66           N
ANISOU   59  N   ALA A  73     3019   3240   1590   -496   -570    636       N
ATOM     60  CA  ALA A  73      -1.494  -3.007  15.297  1.00 20.78           C
ANISOU   60  CA  ALA A  73     3052   3312   1531   -448   -482    632       C
ATOM     61  C   ALA A  73      -1.626  -1.490  15.347  1.00 20.78           C
ANISOU   61  C   ALA A  73     3069   3352   1471   -345   -372    612       C
ATOM     62  O   ALA A  73      -1.458  -0.821  14.356  1.00 19.73           O
ANISOU   62  O   ALA A  73     2880   3148   1466   -287   -376    550       O
ATOM     63  CB  ALA A  73      -2.806  -3.638  14.814  1.00 20.21           C
ANISOU   63  CB  ALA A  73     3001   3095   1581   -346   -502    688       C
ATOM     64  N   HIS A  74      -1.886  -0.956  16.539  1.00 20.28           N
ANISOU   64  N   HIS A  74     2845   3313   1546   -227   -507    527       N
ATOM     65  CA  HIS A  74      -2.151   0.462  16.721  1.00 19.54           C
ANISOU   65  CA  HIS A  74     2627   3313   1483   -230   -398    516       C
ATOM     66  C   HIS A  74      -0.912   1.311  16.604  1.00 20.14           C
ANISOU   66  C   HIS A  74     2915   3371   1366   -422   -203    470       C
ATOM     67  O   HIS A  74      -1.040   2.504  16.391  1.00 20.61           O
ANISOU   67  O   HIS A  74     3071   3383   1374   -333   -275    486       O
ATOM     68  CB  HIS A  74      -2.819   0.697  18.100  1.00 19.53           C
ANISOU   68  CB  HIS A  74     2576   3370   1474   -350   -395    536       C
ATOM     69  CG  HIS A  74      -3.892  -0.291  18.378  1.00 20.17           C
ANISOU   69  CG  HIS A  74     2778   3481   1402   -502   -335    465       C
ATOM     70  ND1 HIS A  74      -3.648  -1.473  19.035  1.00 20.32           N
ANISOU   70  ND1 HIS A  74     2804   3504   1412   -589   -338    501       N
ATOM     71  CD2 HIS A  74      -5.167  -0.364  17.933  1.00 20.96           C
ANISOU   71  CD2 HIS A  74     2842   3608   1513   -472   -453    616       C
ATOM     72  CE1 HIS A  74      -4.753  -2.193  19.063  1.00 20.52           C
ANISOU   72  CE1 HIS A  74     2829   3533   1432   -616   -376    493       C
ATOM     73  NE2 HIS A  74      -5.686  -1.541  18.398  1.00 21.01           N
ANISOU   73  NE2 HIS A  74     2823   3553   1607   -639   -538    401       N
ATOM     74  N   GLY A  75       0.272   0.715  16.755  1.00 19.61           N
ANISOU   74  N   GLY A  75     2921   3134   1395   -466   -239    494       N
ATOM     75  CA  GLY A  75       1.529   1.447  16.690  1.00 19.84           C
ANISOU   75  CA  GLY A  75     2854   3233   1449   -434   -313    568       C
ATOM     76  C   GLY A  75       2.053   1.917  18.043  1.00 20.08           C
ANISOU   76  C   GLY A  75     2817   3329   1483   -566   -371    623       C
ATOM     77  O   GLY A  75       1.343   1.866  19.044  1.00 20.58           O
ANISOU   77  O   GLY A  75     2776   3576   1468   -608   -380    598       O
ATOM     78  N   VAL A  76       3.283   2.420  18.046  1.00 20.15           N
ANISOU   78  N   VAL A  76     2801   3358   1497   -543   -361    612       N
ATOM     79  CA  VAL A  76       3.920   2.883  19.288  1.00 20.26           C
ANISOU   79  CA  VAL A  76     2836   3378   1483   -616   -338    616       C
ATOM     80  C   VAL A  76       4.742   4.126  19.020  1.00 20.38           C
ANISOU   80  C   VAL A  76     2724   3525   1491   -675   -311    640       C
ATOM     81  O   VAL A  76       5.161   4.379  17.886  1.00 19.93           O
ANISOU   81  O   VAL A  76     2697   3305   1567   -666   -306    759       O
ATOM     82  CB  VAL A  76       4.819   1.795  19.881  1.00 20.01           C
ANISOU   82  CB  VAL A  76     2771   3329   1500   -616   -357    499       C
ATOM     83  CG1 VAL A  76       3.984   0.589  20.254  1.00 19.79           C
ANISOU   83  CG1 VAL A  76     2688   3366   1464   -624   -298    445       C
ATOM     84  CG2 VAL A  76       5.928   1.389  18.912  1.00 20.07           C
ANISOU   84  CG2 VAL A  76     2741   3346   1537   -552   -350    559       C
ATOM     85  N   SER A  77       4.954   4.909  20.069  1.00 20.67           N
ANISOU   85  N   SER A  77     2833   3545   1475   -699   -282    650       N
ATOM     86  CA  SER A  77       5.808   6.071  19.998  1.00 21.19           C
ANISOU   86  CA  SER A  77     3018   3514   1519   -731   -326    708       C
ATOM     87  C   SER A  77       7.265   5.677  19.776  1.00 21.98           C
ANISOU   87  C   SER A  77     2962   3688   1701   -726   -450    738       C
ATOM     88  O   SER A  77       7.671   4.545  20.040  1.00 21.51           O
ANISOU   88  O   SER A  77     2983   3649   1538   -724   -225    670       O
ATOM     89  CB  SER A  77       5.716   6.869  21.283  1.00 20.97           C
ANISOU   89  CB  SER A  77     2887   3496   1582   -716   -450    654       C
ATOM     90  OG  SER A  77       6.230   6.070  22.318  1.00 20.53           O
ANISOU   90  OG  SER A  77     2828   3335   1636   -702   -509    561       O
TER
END
"""

def exercise_00(pdb_str, formula):
  """
  TLS group 6 of 4muy.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  asc = pdb_hierarchy.atom_selection_cache()
  cs = pdb_inp.crystal_symmetry_from_cryst1()
  tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy    = pdb_hierarchy)
  #
  deg_to_rad_scale = math.pi/180
  tls_params_one_group = tls_extract.tls_params[0]
  T = matrix.sym(sym_mat3=tls_params_one_group.t)
  L = matrix.sym(sym_mat3=tls_params_one_group.l)
  S = matrix.sqr(tls_params_one_group.s)
  origin = tls_params_one_group.origin
  tlso = tools.tlso(
    t      = T.as_sym_mat3(),
    l      = L.as_sym_mat3(),
    s      = S,
    origin = origin)
  log = open("analysis.log","w")
  r = analysis.run(T=T, L=L*(deg_to_rad_scale**2), S=S*deg_to_rad_scale,
    log=log, find_t_S_using_formula=formula).self_check(show=False)
  log.close()
  #
  rs = flex.double()
  for trial in range(10):
    o = tools.u_tls_vs_u_ens(pdb_str=pdb_str,
      dx       = r.dx,
      dy       = r.dy,
      dz       = r.dz,
      sx       = r.sx,
      sy       = r.sy,
      sz       = r.sz,
      lx       = r.l_x,
      ly       = r.l_y,
      lz       = r.l_z,
      tx       = r.tx,
      ty       = r.ty,
      tz       = r.tz,
      vx       = r.v_x,
      vy       = r.v_y,
      vz       = r.v_z,
      w_M_lx   = r.w_M_lx,
      w_M_ly   = r.w_M_ly,
      w_M_lz   = r.w_M_lz,
      origin   = origin,
      n_models = 10000,
      assert_similarity=False)
    rs.append(o.r)
  return flex.mean(rs)

if (__name__ == "__main__"):
  t0 = time.time()
  for formula in ["10","11"]:
    print("formula:", formula)
    r = exercise_00(pdb_str=pdb_str, formula=formula)
    print("  ", r)
    if(formula=="10"): assert r>0.55
    if(formula=="11"): assert r<0.06
