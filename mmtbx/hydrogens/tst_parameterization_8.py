from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
import iotbx.phil
from mmtbx.hydrogens import connectivity
from libtbx.utils import null_out
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from six.moves import zip


def exercise1():

  pdb_str = """
CRYST1   21.850   24.325   24.603  90.00  90.00  90.00 P 1
SCALE1      0.045767  0.000000  0.000000        0.00000
SCALE2      0.000000  0.041110  0.000000        0.00000
SCALE3      0.000000  0.000000  0.040645        0.00000
ATOM      1  N   CYS D  37       8.132  18.794  15.343  1.00143.27      D    N
ANISOU    1  N   CYS D  37    26762  17628  10046  -3652   3033    538  D    N
ATOM      2  CA  CYS D  37       9.229  18.179  14.592  1.00143.16      D    C
ANISOU    2  CA  CYS D  37    26739  17512  10144  -3661   2916    574  D    C
ATOM      3  C   CYS D  37      10.551  18.924  14.744  1.00143.18      D    C
ANISOU    3  C   CYS D  37    26760  17474  10166  -3627   2815    455  D    C
ATOM      4  O   CYS D  37      11.615  18.295  14.696  1.00141.47      D    O
ANISOU    4  O   CYS D  37    26571  17211   9972  -3636   2700    489  D    O
ATOM      5  CB  CYS D  37       8.857  18.048  13.112  1.00144.61      D    C
ANISOU    5  CB  CYS D  37    26841  17595  10509  -3672   2948    608  D    C
ATOM      6  SG  CYS D  37       7.159  17.517  12.843  1.00142.64      D    S
ANISOU    6  SG  CYS D  37    26553  17390  10253  -3704   3086    716  D    S
ATOM      8  HA  CYS D  37       9.365  17.282  14.934  1.00171.95      D    H
ATOM      9  HB2 CYS D  37       8.970  18.910  12.683  1.00173.69      D    H
ATOM     10  HB3 CYS D  37       9.442  17.394  12.699  1.00173.69      D    H
ATOM     11  N   CYS D  40      13.775  17.698  16.644  1.00154.41      D    N
ANISOU   11  N   CYS D  40    28359  18912  11399  -3629   2475    474  D    N
ATOM     12  CA  CYS D  40      14.636  16.639  16.128  1.00149.28      D    C
ANISOU   12  CA  CYS D  40    27712  18183  10825  -3653   2369    564  D    C
ATOM     13  C   CYS D  40      15.084  16.830  14.685  1.00148.00      D    C
ANISOU   13  C   CYS D  40    27476  17888  10868  -3647   2342    540  D    C
ATOM     14  O   CYS D  40      16.158  16.339  14.335  1.00149.38      D    O
ANISOU   14  O   CYS D  40    27657  17989  11113  -3652   2233    563  D    O
ATOM     15  CB  CYS D  40      13.921  15.279  16.208  1.00140.80      D    C
ANISOU   15  CB  CYS D  40    26653  17143   9703  -3698   2398    729  D    C
ATOM     16  SG  CYS D  40      13.262  14.779  17.830  1.00136.64      D    S
ANISOU   16  SG  CYS D  40    26211  16771   8936  -3714   2439    794  D    S
ATOM     18  HA  CYS D  40      15.432  16.588  16.680  1.00179.29      D    H
ATOM     19  HB2 CYS D  40      13.174  15.295  15.590  1.00169.12      D    H
ATOM     20  HB3 CYS D  40      14.548  14.592  15.932  1.00169.12      D    H
HETATM   21  NB  HEC D1001       9.421  13.220  11.858  1.00140.59           N
ANISOU   21  NB  HEC D1001    26354  16953  10112  -3810   2765   1093       N
HETATM   22  ND  HEC D1001      12.886  11.034  12.753  1.00144.24           N
ANISOU   22  ND  HEC D1001    26954  17334  10517  -3839   2396   1214       N
HETATM   23  C1A HEC D1001      12.227  10.592   9.833  1.00139.00           C
ANISOU   23  C1A HEC D1001    26141  16470  10203  -3864   2452   1283       C
HETATM   24  C1B HEC D1001       8.860  13.248  10.595  1.00143.63           C
ANISOU   24  C1B HEC D1001    26666  17264  10645  -3817   2816   1110       C
HETATM   25  C1C HEC D1001      10.125  13.804  14.697  1.00139.27           C
ANISOU   25  C1C HEC D1001    26337  16981   9599  -3780   2714   1003       C
HETATM   26  C1D HEC D1001      13.435  10.960  14.028  1.00134.69           C
ANISOU   26  C1D HEC D1001    25819  16201   9156  -3833   2339   1204       C
HETATM   27  C2A HEC D1001      11.928  10.386   8.429  1.00142.28           C
ANISOU   27  C2A HEC D1001    26487  16789  10784  -3875   2476   1314       C
HETATM   28  C2B HEC D1001       7.789  14.214  10.617  1.00141.94           C
ANISOU   28  C2B HEC D1001    26417  17097  10418  -3798   2937   1044       C
HETATM   29  C2C HEC D1001      10.451  14.023  16.100  1.00140.94           C
ANISOU   29  C2C HEC D1001    26625  17289   9635  -3766   2687    965       C
HETATM   30  C2D HEC D1001      14.546  10.022  13.991  1.00142.48           C
ANISOU   30  C2D HEC D1001    26834  17133  10168  -3850   2214   1270       C
HETATM   31  C3A HEC D1001      10.804  11.059   8.149  1.00136.41           C
ANISOU   31  C3A HEC D1001    25701  16075  10056  -3866   2591   1279       C
HETATM   32  C3B HEC D1001       7.744  14.781  11.826  1.00143.82           C
ANISOU   32  C3B HEC D1001    26706  17431  10509  -3780   2957    983       C
HETATM   33  C3C HEC D1001      11.528  13.271  16.400  1.00137.87           C
ANISOU   33  C3C HEC D1001    26282  16879   9224  -3777   2567   1011       C
HETATM   34  C3D HEC D1001      14.632   9.534  12.538  1.00144.35           C
ANISOU   34  C3D HEC D1001    27005  17250  10590  -3867   2201   1319       C
HETATM   35  C4A HEC D1001      10.358  11.714   9.366  1.00144.69           C
ANISOU   35  C4A HEC D1001    26789  17236  10949  -3849   2645   1224       C
HETATM   36  C4B HEC D1001       8.768  14.145  12.649  1.00137.78           C
ANISOU   36  C4B HEC D1001    26012  16687   9651  -3787   2848   1015       C
HETATM   37  C4C HEC D1001      11.920  12.537  15.214  1.00134.63           C
ANISOU   37  C4C HEC D1001    25826  16356   8971  -3798   2515   1081       C
HETATM   38  C4D HEC D1001      13.570  10.218  11.849  1.00138.72           C
ANISOU   38  C4D HEC D1001    26232  16533   9944  -3858   2317   1280       C
HETATM   39  CAA HEC D1001      12.769   9.532   7.451  1.00132.04           C
ANISOU   39  CAA HEC D1001    25171  15378   9621  -3893   2385   1376       C
HETATM   40  CAB HEC D1001       6.719  15.866  12.240  1.00138.71           C
ANISOU   40  CAB HEC D1001    26042  16858   9803  -3757   3079    901       C
HETATM   41  CAC HEC D1001      12.297  13.190  17.738  1.00135.76           C
ANISOU   41  CAC HEC D1001    26101  16686   8794  -3768   2488    993       C
HETATM   42  CAD HEC D1001      15.637   8.535  11.925  1.00138.75           C
ANISOU   42  CAD HEC D1001    26295  16445   9980  -3888   2089   1394       C
HETATM   43  CBA HEC D1001      11.998   8.251   7.139  1.00132.66           C
ANISOU   43  CBA HEC D1001    25237  15456   9711  -3938   2416   1524       C
HETATM   44  CBB HEC D1001       5.650  15.212  13.142  1.00135.28           C
ANISOU   44  CBB HEC D1001    25642  16531   9227  -3784   3158    997       C
HETATM   45  CBC HEC D1001      11.427  12.956  18.987  1.00146.99           C
ANISOU   45  CBC HEC D1001    27580  18242  10026  -3778   2561   1040       C
HETATM   46  CBD HEC D1001      14.919   7.222  11.642  1.00139.54           C
ANISOU   46  CBD HEC D1001    26386  16545  10089  -3933   2120   1549       C
HETATM   47  CGA HEC D1001      12.527   7.624   5.875  1.00128.84           C
ANISOU   47  CGA HEC D1001    24714  14851   9389  -3953   2359   1570       C
HETATM   48  CGD HEC D1001      15.891   6.207  11.098  1.00141.03           C
ANISOU   48  CGD HEC D1001    26575  16642  10368  -3955   2011   1625       C
HETATM   49  CHA HEC D1001      13.299  10.063  10.513  1.00146.02           C
ANISOU   49  CHA HEC D1001    27088  17363  11029  -3868   2342   1305       C
HETATM   50  CHB HEC D1001       9.227  12.490   9.500  1.00144.37           C
ANISOU   50  CHB HEC D1001    26723  17255  10876  -3837   2764   1179       C
HETATM   51  CHC HEC D1001       9.083  14.384  13.975  1.00144.74           C
ANISOU   51  CHC HEC D1001    26964  17660  10370  -3774   2826    976       C
HETATM   52  CHD HEC D1001      12.978  11.644  15.133  1.00137.25           C
ANISOU   52  CHD HEC D1001    26183  16634   9330  -3815   2394   1143       C
HETATM   53  CMA HEC D1001      10.090  11.136   6.782  1.00134.85           C
ANISOU   53  CMA HEC D1001    25424  15802  10009  -3872   2656   1294       C
HETATM   54  CMB HEC D1001       6.909  14.580   9.402  1.00138.00           C
ANISOU   54  CMB HEC D1001    25835  16541  10059  -3798   3022   1038       C
HETATM   55  CMC HEC D1001       9.714  15.017  17.030  1.00140.17           C
ANISOU   55  CMC HEC D1001    26549  17298   9410  -3741   2781    878       C
HETATM   56  CMD HEC D1001      15.472   9.587  15.143  1.00143.64           C
ANISOU   56  CMD HEC D1001    27061  17330  10186  -3851   2112   1287       C
HETATM   57  NA  HEC D1001      11.246  11.403  10.381  1.00142.02           N
ANISOU   57  NA  HEC D1001    26521  16938  10501  -3848   2558   1228       N
HETATM   58  NC  HEC D1001      11.031  12.874  14.214  1.00135.40           N
ANISOU   58  NC  HEC D1001    25852  16413   9183  -3799   2607   1074       N
HETATM   59  O1A HEC D1001      13.768   7.499   5.723  1.00127.91           O
ANISOU   59  O1A HEC D1001    24612  14669   9319  -3944   2254   1548       O
HETATM   60  O1D HEC D1001      15.705   5.000  11.405  1.00143.73           O
ANISOU   60  O1D HEC D1001    26943  17007  10662  -3991   1995   1754       O
HETATM   61  O2A HEC D1001      11.708   7.243   5.000  1.00126.20           O
ANISOU   61  O2A HEC D1001    24330  14482   9138  -3973   2418   1629       O
HETATM   62  O2D HEC D1001      16.850   6.595  10.375  1.00138.38           O
ANISOU   62  O2D HEC D1001    26214  16212  10153  -3936   1943   1558       O
HETATM   63 FE   HEC D1001      11.115  12.157  12.415  1.00158.39          Fe
ANISOU   63 FE   HEC D1001    28683  19178  12321  -3824   2585   1151      Fe
HETATM   64  HAB HEC D1001       6.266  16.109  11.430  1.00166.61           H
HETATM   65  HAC HEC D1001      13.016  12.556  17.688  1.00163.06           H
HETATM   66  HHA HEC D1001      13.931   9.523   9.993  1.00175.38           H
HETATM   67  HHB HEC D1001       8.617  12.507   8.733  1.00173.40           H
HETATM   68  HHC HEC D1001       8.518  15.025  14.454  1.00173.84           H
HETATM   69  HHD HEC D1001      13.462  11.479  15.969  1.00164.85           H
HETATM   70 HAA1 HEC D1001      13.619   9.309   7.860  1.00158.61           H
HETATM   71 HAA2 HEC D1001      12.922  10.028   6.632  1.00158.61           H
HETATM   72 HAD1 HEC D1001      16.364   8.380  12.549  1.00166.66           H
HETATM   73 HAD2 HEC D1001      15.991   8.897  11.097  1.00166.66           H
HETATM   74 HBA1 HEC D1001      11.058   8.462   7.024  1.00159.35           H
HETATM   75 HBA2 HEC D1001      12.100   7.627   7.874  1.00159.35           H
HETATM   76 HBB1 HEC D1001       5.203  14.504  12.651  1.00162.49           H
HETATM   77 HBB2 HEC D1001       5.000  15.880  13.409  1.00162.49           H
HETATM   78 HBB3 HEC D1001       6.076  14.841  13.930  1.00162.49           H
HETATM   79 HBC1 HEC D1001      11.697  12.130  19.419  1.00176.54           H
HETATM   80 HBC2 HEC D1001      10.495  12.894  18.725  1.00176.54           H
HETATM   81 HBC3 HEC D1001      11.541  13.696  19.603  1.00176.54           H
HETATM   82 HBD1 HEC D1001      14.216   7.374  10.991  1.00167.61           H
HETATM   83 HBD2 HEC D1001      14.529   6.885  12.463  1.00167.61           H
HETATM   84 HMA1 HEC D1001       9.181  10.779   6.867  1.00161.97           H
HETATM   85 HMA2 HEC D1001      10.588  10.607   6.123  1.00161.97           H
HETATM   86 HMA3 HEC D1001      10.048  12.069   6.487  1.00161.97           H
HETATM   87 HMB1 HEC D1001       7.003  15.525   9.208  1.00165.76           H
HETATM   88 HMB2 HEC D1001       5.981  14.382   9.604  1.00165.76           H
HETATM   89 HMB3 HEC D1001       7.190  14.061   8.632  1.00165.76           H
HETATM   90 HMC1 HEC D1001       9.327  14.535  17.777  1.00168.36           H
HETATM   91 HMC2 HEC D1001       9.010  15.463  16.534  1.00168.36           H
HETATM   92 HMC3 HEC D1001      10.343  15.677  17.362  1.00168.36           H
HETATM   93 HMD1 HEC D1001      16.399   9.817  14.922  1.00172.53           H
HETATM   94 HMD2 HEC D1001      15.399   8.618  15.274  1.00172.53           H
HETATM   95 HMD3 HEC D1001      15.209  10.048  15.966  1.00172.53           H
TER
END
  """

  edits = """
geometry_restraints.edits {
  bond {
    atom_selection_1 = chain D and resseq 40 and name SG
    atom_selection_2 = chain D and resseq 1001 and name CAC
    distance_ideal = 1.81
    sigma = 0.05
  }
  bond {
    atom_selection_1 = chain D and resseq 37 and name SG
    atom_selection_2 = chain D and resseq 1001 and name CAB
    distance_ideal = 1.81
    sigma = 0.05
  }
  angle {
    atom_selection_1 = chain D and resseq 40 and name SG
    atom_selection_2 = chain D and resseq 1001 and name CAC
    atom_selection_3 = chain D and resseq 1001 and name CBC
    angle_ideal = 109
    sigma = 2
  }
  angle {
    atom_selection_1 = chain D and resseq 40 and name SG
    atom_selection_2 = chain D and resseq 1001 and name CAC
    atom_selection_3 = chain D and resseq 1001 and name C3C
    angle_ideal = 109
    sigma = 2
  }
  angle {
    atom_selection_1 = chain D and resseq 37 and name SG
    atom_selection_2 = chain D and resseq 1001 and name CAB
    atom_selection_3 = chain D and resseq 1001 and name CBB
    angle_ideal = 109
    sigma = 2
  }
  angle {
    atom_selection_1 = chain D and resseq 1001 and name HAB
    atom_selection_2 = chain D and resseq 1001 and name CAB
    atom_selection_3 = chain D and resseq 37 and name SG
    angle_ideal = 109
    sigma = 1.1
  }
  angle {
    atom_selection_1 = chain D and resseq 1001 and name HAC
    atom_selection_2 = chain D and resseq 1001 and name CAC
    atom_selection_3 = chain D and resseq 40 and name SG
    angle_ideal = 109
    sigma = 1.1
  }
}
  """

  gm_phil = iotbx.phil.parse(
      input_string     = grand_master_phil_str,
      process_includes = True)
  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()

  # Make sure the angle edit is present
  assert (params.geometry_restraints.edits.angle[3].atom_selection_1 == \
    "chain D and resseq 1001 and name HAB")

  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)

  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params=params,
    make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = sites_cart,
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  type_list     = diagnostics.type_list

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_para) - h_para.count(None)

  for rc in h_para:
    if rc:
      assert(rc.ih != 61), 'Wrong atom not recognized.'

# Test if number of paramterized H atoms is correct
  assert (number_h == number_h_para + 1), 'Not all H atoms are parameterized'

  type_list_known = ['3neigbs', '2tetra', '2tetra', '3neigbs',
    '2tetra', '2tetra', '3neigbs', 'flat_2neigbs', 'flat_2neigbs',
    'flat_2neigbs', 'flat_2neigbs', '2tetra', '2tetra', '2tetra', '2tetra',
    '2tetra', '2tetra', 'prop', 'prop', 'prop', 'prop', 'prop', 'prop',
    '2tetra', '2tetra', 'prop', 'prop', 'prop', 'prop', 'prop', 'prop',
    'prop', 'prop', 'prop', 'prop', 'prop', 'prop']

  for ih in h_distances:
    # One H atom is expected to be far (HAC)
    if (ih == 62):
      continue
    labels = atoms[ih].fetch_labels()
    if (h_distances[ih] > 0.1):
      assert (h_distances[ih] < 0.1), \
        'distance too large: %s  atom: %s (%s) residue: %s ' \
        % (h_para[ih].htype, atoms[ih].name, ih, labels.resseq.strip())
#
  for type1, type2 in zip(type_list, type_list_known):
    assert (type1 == type2)
    #print "'%s'," % type1,

def exercise2():
  pdb_str = """
CRYST1   16.660   12.742   18.240  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.078481  0.000000        0.00000
SCALE3      0.000000  0.000000  0.054825        0.00000
ATOM      1  N   TYR A   7       9.837   5.000   6.625  1.00 15.00           N
ATOM      2  CA  TYR A   7      10.084   6.426   6.798  1.00 15.00           C
ATOM      3  C   TYR A   7      11.431   6.813   6.197  1.00 15.00           C
ATOM      4  O   TYR A   7      11.660   6.642   5.000  1.00 15.00           O
ATOM      5  CB  TYR A   7      10.042   6.803   8.281  1.00 15.00           C
ATOM      6  CG  TYR A   7       8.697   6.593   8.948  1.00 15.00           C
ATOM      7  CD1 TYR A   7       7.540   6.413   8.198  1.00 15.00           C
ATOM      8  CD2 TYR A   7       8.586   6.575  10.332  1.00 15.00           C
ATOM      9  CE1 TYR A   7       6.315   6.222   8.807  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.364   6.384  10.950  1.00 15.00           C
ATOM     11  CZ  TYR A   7       6.233   6.208  10.183  1.00 15.00           C
ATOM     12  OH  TYR A   7       5.015   6.018  10.794  1.00 15.00           O
ATOM     13  HA  TYR A   7       9.398   6.930   6.331  1.00 15.00           H
ATOM     14  HB2 TYR A   7      10.693   6.264   8.757  1.00 15.00           H
ATOM     15  HB3 TYR A   7      10.270   7.742   8.369  1.00 15.00           H
ATOM     16  HD1 TYR A   7       7.589   6.422   7.269  1.00 15.00           H
ATOM     17  HD2 TYR A   7       9.347   6.694  10.853  1.00 15.00           H
ATOM     18  HE1 TYR A   7       5.550   6.103   8.292  1.00 15.00           H
ATOM     19  HE2 TYR A   7       7.306   6.375  11.878  1.00 15.00           H
ATOM     20  HH  TYR A   7       5.000   6.415  11.534  1.00 15.00           H
TER
HETATM   21  O   HOH B   1       5.307   7.545  13.240  1.00 30.00           O
TER
END
  """

  edits = """
geometry_restraints.edits {
  bond {
    atom_selection_1 = chain A and resseq 7 and name HH
    atom_selection_2 = chain B and resseq 1 and name O
    distance_ideal = 1.81
    sigma = 0.05
  }
}
  """

  type_list_known = ['3neigbs', '2tetra', '2tetra', 'flat_2neigbs',
    'flat_2neigbs', 'flat_2neigbs', 'flat_2neigbs', 'alg1b']


  gm_phil = iotbx.phil.parse(
      input_string     = grand_master_phil_str,
      process_includes = True)
  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()

  # Make sure the angle edit is present
  assert (params.geometry_restraints.edits.bond[0].atom_selection_1 == \
    "chain A and resseq 7 and name HH")

  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)

  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params=params,
    make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = sites_cart,
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  type_list     = diagnostics.type_list

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_para) - h_para.count(None)

  connectivity_manager = connectivity.determine_connectivity(
    pdb_hierarchy       = pdb_hierarchy,
    geometry_restraints = model.get_restraints_manager().geometry)
  double_H = connectivity_manager.double_H

# Test if number of paramterized H atoms is correct
  assert (number_h == number_h_para), 'Not all H atoms are parameterized'
  assert (double_H[19] == [11, 20]), 'H bound to two atoms wrongly recognized'
  assert (number_h_para == 8), 'Not all H atoms are parameterized'

  for ih in h_distances:
   labels = atoms[ih].fetch_labels()
   if (h_distances[ih] > 0.1):
     assert (h_distances[ih] < 0.1), \
       'distance too large: %s  atom: %s (%s) residue: %s ' \
       % (h_para[ih].htype, atoms[ih].name, ih, labels.resseq.strip())

  for type1, type2 in zip(type_list, type_list_known):
    assert (type1 == type2)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise1()
  exercise2()
  print("OK. Time: %8.3f"%(time.time()-t0))
