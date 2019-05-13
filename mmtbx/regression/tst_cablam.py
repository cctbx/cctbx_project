from __future__ import division
from __future__ import print_function
from mmtbx.cablam import cablam_validate
from libtbx.test_utils import show_diff
from iotbx import pdb
import libtbx.load_env
import os

ref_cablam_give_text = """residue : outlier_type : contour_level : ca_contour_level : sec struc recommendation : alpha score : beta score : three-ten score
     3  ILE:                    :0.52764:0.49162:                 :0.06929:0.00000:0.00000
     4  PHE:                    :0.83948:0.80426: try alpha helix :0.34122:0.00000:0.00936
     5  GLU:                    :0.81271:0.79155: try alpha helix :0.31550:0.00000:0.00445
     6  MET:                    :0.77295:0.73578: try alpha helix :0.20871:0.00000:0.01115
     7  LEU:                    :0.70530:0.69504: try alpha helix :0.12026:0.00000:0.00484
     8  ARG:                    :0.71859:0.72353: try alpha helix :0.18161:0.00000:0.08888
     9  ILE:                    :0.84025:0.82510: try alpha helix :0.43515:0.00000:0.00548
    10  ASP:                    :0.65772:0.66407: try alpha helix :0.11727:0.00000:0.00000
    11  GLU:                    :0.17045:0.15375:                 :0.01055:0.00000:0.00000
    12  GLY:                    :0.35640:0.49471:                 :0.00000:0.00000:0.00000
    13  LEU:                    :0.07149:0.09040:                 :0.00000:0.00102:0.00000
    14  ARG:                    :0.41594:0.41432: try beta sheet  :0.00000:0.15626:0.00000
    15  LEU:                    :0.37150:0.44845: try beta sheet  :0.01142:0.00202:0.00000
    16  LYS:                    :0.14341:0.12064: try beta sheet  :0.00034:0.00045:0.00000
    17  ILE: CaBLAM Outlier     :0.00777:0.25439: try beta sheet  :0.00000:0.00123:0.00000
    18  TYR:                    :0.38403:0.26228: try beta sheet  :0.00000:0.00048:0.00000
    19  LYS:                    :0.23229:0.31769:                 :0.00000:0.00601:0.00000
    20  ASP:                    :0.24246:0.15915:                 :0.00000:0.00000:0.00000
    21  THR:                    :0.07474:0.10616:                 :0.00639:0.00000:0.00000
    22  GLU:                    :0.11848:0.17442:                 :0.00163:0.00000:0.00000
    23  GLY:                    :0.87986:0.86098:                 :0.00000:0.00000:0.00000
    24  TYR:                    :0.26671:0.25287:                 :0.00000:0.00975:0.00000
    25  TYR:                    :0.17809:0.17537: try beta sheet  :0.00000:0.00060:0.00000
    26  THR:                    :0.26444:0.16084:                 :0.00000:0.04956:0.00000
    27  ILE:                    :0.13051:0.05851:                 :0.00003:0.00000:0.00000
    28  GLY:                    :0.57593:0.34210:                 :0.00218:0.00000:0.00000
    29  ILE: CaBLAM Outlier     :0.00287:0.02229:                 :0.00006:0.00000:0.00000
    30  GLY:                    :0.27781:0.21754:                 :0.00000:0.00000:0.00000
    31  HIS:                    :0.38269:0.14269:                 :0.00000:0.05822:0.00000
    32  LEU:                    :0.34286:0.57402:                 :0.00000:0.07717:0.00000
    40C ASP:                    :0.94906:0.96373:                 :0.76901:0.00000:0.01862
    41  ALA:                    :0.92911:0.76129: try alpha helix :0.63484:0.00000:0.00712
    42  ALA:                    :0.96131:0.96875: try alpha helix :0.77091:0.00000:0.01214
    43  LYS:                    :0.97552:0.87406: try alpha helix :0.83796:0.00000:0.03254
    44  SER:                    :0.78399:0.86899: try alpha helix :0.79994:0.00000:0.01036
    45  GLU:                    :0.87884:0.86969: try alpha helix :0.73818:0.00000:0.01058
    46  LEU:                    :0.91916:0.88671: try alpha helix :0.74961:0.00000:0.03756
    47  ASP:                    :0.97483:0.96910: try alpha helix :0.90420:0.00000:0.01721
    48  LYS:                    :0.91878:0.96392: try alpha helix :0.80659:0.00000:0.01638
    49  ALA:                    :0.72801:0.71401: try alpha helix :0.14364:0.00000:0.00469
    50  ILE:                    :0.37448:0.38935:                 :0.04203:0.00000:0.00000
    51  GLY:                    :0.53355:0.56969:                 :0.00000:0.00000:0.00000
    52  ARG:                    :0.11630:0.07344:                 :0.00000:0.00736:0.00000
    53  ASN:                    :0.11365:0.08824: try beta sheet  :0.00000:0.00088:0.00000
    54  THR:                    :0.05756:0.15570:                 :0.00000:0.07997:0.00000
    55  ASN: CaBLAM Disfavored  :0.02086:0.00975:                 :0.00010:0.00000:0.00000
    56  GLY:                    :0.64275:0.61690:                 :0.00000:0.00000:0.00000
    57  VAL:                    :0.26151:0.15786:                 :0.00000:0.00000:0.00000
    58  ILE:                    :0.44409:0.45268:                 :0.00000:0.00268:0.00000
    59  THR:                    :0.49656:0.40961:                 :0.00000:0.00063:0.00000
    60  LYS:                    :0.57162:0.65356:                 :0.12322:0.00000:0.00000
    61  ASP:                    :0.87174:0.84785: try alpha helix :0.45478:0.00000:0.04867
    62  GLU:                    :0.94925:0.94200: try alpha helix :0.91876:0.00000:0.01601
    63  ALA:                    :0.78514:0.94158: try alpha helix :0.89034:0.00000:0.02358
    64  GLU:                    :0.84156:0.78750: try alpha helix :0.35344:0.00000:0.04013
    65  LYS:                    :0.80299:0.80830: try alpha helix :0.36710:0.00000:0.03984
    66  LEU:                    :0.86889:0.86756: try alpha helix :0.48381:0.00000:0.00890
    67  PHE:                    :0.90189:0.88748: try alpha helix :0.54551:0.00000:0.04663
    68  ASN:                    :0.98865:0.95803: try alpha helix :0.87462:0.00000:0.01623
    69  GLN:                    :0.97655:0.97438: try alpha helix :0.88753:0.00000:0.01169
    70  ASP:                    :0.78162:0.77513: try alpha helix :0.96064:0.00000:0.01640
    71  VAL:                    :0.75300:0.94269: try alpha helix :0.96236:0.00000:0.01399
    72  ASP:                    :0.96009:0.96323: try alpha helix :0.89004:0.00000:0.01905
    73  ALA:                    :0.89729:0.86971: try alpha helix :0.74118:0.00000:0.02781
    74  ALA:                    :0.96049:0.91354: try alpha helix :0.79346:0.00000:0.02381
    75  VAL:                    :0.90762:0.90773: try alpha helix :0.63817:0.00000:0.00743
    76  ARG:                    :0.85886:0.91395: try alpha helix :0.72908:0.00000:0.01009
    77  GLY:                    :0.96700:0.98047: try alpha helix :0.73374:0.00000:0.03867
    78  ILE:                    :0.98019:0.98238: try alpha helix :0.87047:0.00000:0.01631
    79  LEU:                    :0.58012:0.90368: try alpha helix :0.72839:0.00000:0.00859
    80  ARG:                    :0.44662:0.55775:                 :0.01191:0.00000:0.00000
    81  ASN:                    :0.34785:0.31423:                 :0.00000:0.00031:0.00000
    82  ALA:                    :0.29269:0.38361:                 :0.02096:0.00000:0.00000
    83  LYS:                    :0.58211:0.71041: try alpha helix :0.14248:0.00000:0.00001
    84  LEU:                    :0.41959:0.32368: try alpha helix :0.05148:0.00000:0.00626
    85  LYS:                    :0.66986:0.51927: try three-ten   :0.09610:0.00000:0.21238
    86  PRO:                    :0.98883:0.99594: try alpha helix :0.86833:0.00000:0.01204
    87  VAL:                    :0.75880:0.95367: try alpha helix :0.78186:0.00000:0.00912
    88  TYR:                    :0.91356:0.90539: try alpha helix :0.58192:0.00000:0.02264
    89  ASP:                    :0.73904:0.80679: try alpha helix :0.35668:0.00000:0.05555
    90  SER:                    :0.55972:0.60595: try three-ten   :0.00733:0.00000:0.01601
    91  LEU:                    :0.27769:0.22558:                 :0.00000:0.00107:0.00000
    92  ASP:                    :0.42270:0.36276:                 :0.00000:0.00000:0.00000
    93  ALA:                    :0.45818:0.40113:                 :0.07442:0.00000:0.00000
    94  VAL:                    :0.77606:0.76480: try alpha helix :0.51294:0.00000:0.00679
    95  ARG:                    :0.73578:0.72447: try alpha helix :0.16465:0.00000:0.02319
    96  ARG:                    :0.73291:0.70074: try three-ten   :0.14308:0.00000:0.15348
    97  ALA:                    :0.83719:0.80947: try alpha helix :0.32626:0.00000:0.05371
    98  ALA:                    :0.86385:0.88142: try alpha helix :0.53840:0.00000:0.04449
    99  LEU:                    :0.97900:0.93347: try alpha helix :0.84827:0.00000:0.01657
   100  ILE:                    :0.83575:0.81549: try alpha helix :0.35748:0.00000:0.04442
   101  ASN:                    :0.87477:0.84323: try alpha helix :0.45214:0.00000:0.04502
   102  MET:                    :0.96302:0.94884: try alpha helix :0.86991:0.00000:0.01048
   103  VAL:                    :0.96613:0.92822: try alpha helix :0.78662:0.00000:0.01923
   104  PHE:                    :0.77634:0.74271: try alpha helix :0.21660:0.00000:0.00765
   105  GLN:                    :0.56557:0.55357: try alpha helix :0.09128:0.00000:0.00000
   106  MET:                    :0.17345:0.14246:                 :0.01504:0.00000:0.00000
   107  GLY:                    :0.45772:0.54732:                 :0.00000:0.00000:0.00000
   108  GLU:                    :0.46203:0.46061:                 :0.05310:0.00000:0.00000
   109  THR:                    :0.92592:0.77315: try alpha helix :0.79332:0.00000:0.03679
   110  GLY:                    :0.94492:0.92742: try alpha helix :0.34326:0.00000:0.00675
   111  VAL:                    :0.65875:0.68327: try alpha helix :0.11230:0.00000:0.01439
   112  ALA:                    :0.68953:0.69294: try three-ten   :0.05884:0.00000:0.37750
   113  GLY:                    :0.35710:0.61555: try three-ten   :0.00122:0.00000:0.06679
   114  PHE: CaBLAM Disfavored  :0.02716:0.02012:                 :0.00000:0.00000:0.00000
   115  THR:                    :0.21716:0.25651:                 :0.02337:0.00038:0.00000
   116  ASN:                    :0.93064:0.87429: try alpha helix :0.62433:0.00000:0.00578
   117  SER:                    :0.84801:0.85414: try alpha helix :0.42589:0.00000:0.02916
   118  LEU:                    :0.82109:0.80540: try alpha helix :0.39810:0.00000:0.05275
   119  ARG:                    :0.97488:0.91843: try alpha helix :0.84079:0.00000:0.02236
   120  MET:                    :0.97381:0.94393: try alpha helix :0.85304:0.00000:0.02222
   121  LEU:                    :0.96664:0.96247: try alpha helix :0.76570:0.00000:0.01682
   122  GLN:                    :0.81647:0.96855: try alpha helix :0.78504:0.00000:0.03365
   123  GLN:                    :0.38558:0.54762:                 :0.01131:0.00000:0.00000
   124  LYS:                    :0.27195:0.31297:                 :0.00000:0.00063:0.00000
   125  ARG:                    :0.15025:0.03510: try beta sheet  :0.00000:0.01617:0.00000
   126  TRP:                    :0.36307:0.24379:                 :0.01956:0.00059:0.00000
   127  ASP:                    :0.90607:0.87620: try alpha helix :0.66236:0.00000:0.00567
   128  GLU:                    :0.94643:0.72509: try alpha helix :0.87977:0.00000:0.01864
   129  ALA:                    :0.78394:0.78865: try alpha helix :0.29204:0.00000:0.04696
   130  ALA:                    :0.83281:0.80403: try alpha helix :0.35091:0.00000:0.04734
   131  VAL:                    :0.89150:0.84481: try alpha helix :0.42877:0.00000:0.00505
   132  ASN:                    :0.84250:0.79163: try alpha helix :0.36966:0.00000:0.01480
   133  LEU:                    :0.72732:0.71627: try three-ten   :0.12666:0.00000:0.12935
   134  ALA:                    :0.58761:0.61900: try three-ten   :0.00324:0.00000:0.37544
   135  LYS: CaBLAM Disfavored  :0.04655:0.16709: try three-ten   :0.00000:0.00000:0.02864
   136  SER:                    :0.21418:0.10964:                 :0.00000:0.00275:0.00000
   137  ARG:                    :0.44964:0.57604:                 :0.05596:0.00000:0.00000
   138  TRP:                    :0.72494:0.72727: try alpha helix :0.15493:0.00000:0.07705
   139  TYR:                    :0.84309:0.82430: try alpha helix :0.39345:0.00000:0.00260
   140  ASN:                    :0.72053:0.79855: try alpha helix :0.36505:0.00000:0.00554
   141  GLN:                    :0.53741:0.51793:                 :0.07728:0.00000:0.00000
   142  THR:                    :0.19896:0.17954:                 :0.00051:0.00000:0.00000
   143  PRO:                    :0.66687:0.52699:                 :0.04263:0.00083:0.00000
   144  ASN:                    :0.77072:0.77580: try alpha helix :0.27392:0.00000:0.00914
   145  ARG:                    :0.74389:0.76664: try alpha helix :0.24488:0.00000:0.00073
   146  ALA:                    :0.91235:0.91083: try alpha helix :0.60613:0.00000:0.03068
   147  LYS:                    :0.75173:0.81694: try alpha helix :0.35878:0.00000:0.05334
   148  ARG:                    :0.79680:0.75473: try alpha helix :0.31349:0.00000:0.01607
   149  VAL:                    :0.82542:0.80165: try alpha helix :0.33416:0.00000:0.00323
   150  ILE:                    :0.90417:0.91613: try alpha helix :0.63958:0.00000:0.03949
   151  THR:                    :0.81490:0.82636: try alpha helix :0.48350:0.00000:0.03861
   152  THR:                    :0.86653:0.84523: try alpha helix :0.45410:0.00000:0.02507
   153  PHE:                    :0.88416:0.85966: try alpha helix :0.48595:0.00000:0.00357
   154  ARG:                    :0.58641:0.58326: try alpha helix :0.07799:0.00000:0.00000
   155  THR:                    :0.24379:0.24630:                 :0.02083:0.00000:0.00000
   156  GLY:                    :0.59583:0.57219:                 :0.00000:0.00063:0.00000
   157  THR:                    :0.32804:0.23017:                 :0.00000:0.00482:0.00000
   158  TRP:                    :0.08355:0.28156:                 :0.00000:0.00000:0.00000
   159  ASP:                    :0.40617:0.36789:                 :0.01844:0.00000:0.00000
   160  ALA:                    :0.66012:0.75135:                 :0.20693:0.00000:0.06554
"""

ref_cablam_give_oneline = """pdb103l.ent:151:3.3:1.3:0.00
"""

class cablam_test_string():
  #I wrote the regression test to use a class with a custom .write() method as a
  #  proof of principle for learning OOP and to see if I could. Possible because
  #  all my print functions accept an optional writeto= variable.
  def write(self,string):
    self.output += str(string)
  def __init__(self):
    self.output = ""

def exercise_cablam():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb103l.ent",
    test=os.path.isfile) #This is the same file used for tst_kinemage.py
  if (regression_pdb is None):
    print("Skipping exercise_cablam(): input pdb (pdb103l.ent) not available")
    return
  #-----
  pdb_io = pdb.input(regression_pdb)
  pdbid = os.path.basename(regression_pdb)
  hierarchy = pdb_io.construct_hierarchy()

  oneline_test = cablam_test_string()
  cablam_validate.oneline(hierarchy, peptide_cutoff=0.05,
    peptide_bad_cutoff=0.01, ca_cutoff=0.005, pdbid=pdbid,
    writeto=oneline_test)

  text_test = cablam_test_string()
  outliers = cablam_validate.analyze_pdb(hierarchy, outlier_cutoff=0.05, pdbid=pdbid)
  cablam_validate.give_text(outliers, writeto=text_test)

  #print '|',oneline_test.output,'|'
  #print '|',ref_cablam_give_oneline,'|'
  #print '|',text_test.output,'|'
  #print '|',ref_cablam_give_text,'|'

  assert not show_diff(oneline_test.output , ref_cablam_give_oneline)
  assert not show_diff(text_test.output , ref_cablam_give_text)

def run():
  ##if (not libtbx.env.has_module(name="probe")):
  ##  print \
  ##    "Skipping kinemage test:" \
  ##    " probe not available"
  ##elif (not libtbx.env.has_module(name="reduce")):
  ##  print \
  ##    "Skipping kinemage test:" \
  ##    " reduce not available"
  ##else:
  ##  exercise_kinemage()
  ##print "OK"
  exercise_cablam()
  print("OK")

if (__name__ == "__main__"):
  run()
