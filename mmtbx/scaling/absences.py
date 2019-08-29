
from __future__ import absolute_import, division, print_function
import mmtbx.scaling
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from cctbx import sgtbx
from libtbx import table_utils
import math
from six.moves import zip

# name     hkl selector  condition
absence_and_conditions = {
  "2_0 (a)" : [(None,0,0),  (0 ,0 ,0 ) ],
  "2_0 (b)" : [(0,None,0),  (0,0,0)],
  "2_0 (c)" : [(0,0,None),  (0,0,0)],


  "2_1 (a)" : [(None,0,0),  (1.0/2.0,0,0)],
  "2_1 (b)" : [(0,None,0),  (0,1.0/2.0,0)],
  "2_1 (c)" : [(0,0,None),  (0,0,1.0/2.0)],


  "3_0 (c)" : [(0,0,None),  (0,0,0)],
  "3_1 (c)" : [(0,0,None),  (0,0,1.0/3.0)],
  "3_2 (c)" : [(0,0,None),  (0,0,1.0/3.0)],

  "4_0 (c)" : [(0,0,None),  (0,0,0)],
  "4_1 (c)" : [(0,0,None),  (0,0,1.0/4.0)],
  "4_2 (c)" : [(0,0,None),  (0,0,2.0/4.0)],
  "4_3 (c)" : [(0,0,None),  (0,0,3.0/4.0)],


  "4_0 (a)" : [(None,0,0),  (0,0,0)],
  "4_1 (a)" : [(None,0,0),  (1.0/4.0,0,0)],
  "4_2 (a)" : [(None,0,0),  (2.0/4.0,0,0)],
  "4_3 (a)" : [(None,0,0),  (3.0/4.0,0,0)],


  "4_0 (b)" : [(0,None,0),  (0,0,0)],
  "4_1 (b)" : [(0,None,0),  (0,1.0/4.0,0)],
  "4_2 (b)" : [(0,None,0),  (0,2.0/4.0,0)],
  "4_3 (b)" : [(0,None,0),  (0,3.0/4.0,0)],


  "6_0 (c)" : [(0,0,None),  (0,0,0)],
  "6_1 (c)" : [(0,0,None),  (0,0,1.0/6.0)],
  "6_2 (c)" : [(0,0,None),  (0,0,2.0/6.0)],
  "6_3 (c)" : [(0,0,None),  (0,0,3.0/6.0)],
  "6_4 (c)" : [(0,0,None),  (0,0,4.0/6.0)],
  "6_5 (c)" : [(0,0,None),  (0,0,5.0/6.0)],

  "b (a)"   : [(0,None,None),  (0,1.0/2.0,0)],
  "c (a)"   : [(0,None,None),  (0,0,1.0/2.0)],
  "n (a)"   : [(0,None,None),  (0,1.0/2.0,1.0/2.0)],
  "d (a)"   : [(0,None,None),  (0,1.0/4.0,1.0/4.0)],

  "a (b)"   : [(None,0,None),  (1.0/2.0,0,0)],
  "c (b)"   : [(None,0,None),  (0,0,1.0/2.0)],
  "n (b)"   : [(None,0,None),  (1.0/2.0,0,1.0/2.0)],
  "d (b)"   : [(None,0,None),  (1.0/4.0,0,1.0/4.0)],

  "a (c)"   : [(None,None,0),  (1.0/2.0,0,0)],
  "b (c)"   : [(None,None,0),  (0,1.0/2.0,0)],
  "n (c)"   : [(None,None,0),  (1.0/2.0,1.0/2.0,0)],
  "d (c)"   : [(None,None,0),  (1.0/4.0,1.0/4.0,0)]
}


absence_classes = {
  "along a 2" : ["2_0 (a)", "2_1 (a)"],
  "along a"   : ["b (a)", "c (a)", "n (a)", "d (a)"],
  "along b 4" : ["4_0 (b)", "4_1 (b)", "4_2 (b)", "4_3 (b)"],
  "along b 2" : ["2_0 (b)", "2_1 (b)"],
  "along b"   : ["a (b)", "c (b)", "n (b)", "d (b)"],
  "along c"   : ["a (c)", "b (c)", "n (c)", "d (c)"],
  "along c 2" : ["2_0 (c)", "2_1 (c)"],
  "along c 3" : ["3_0 (c)", "3_1 (c)", "3_2 (c)"],
  "along c 4" : ["4_0 (c)", "4_1 (c)", "4_2 (c)", "4_3 (c)"],
  "along c 6" : ["6_0 (c)", "6_1 (c)", "6_2 (c)", "6_3 (c)", "6_4 (c)", "6_5 (c)"] }


along_a   = absence_classes["along a"]
along_a_2 = absence_classes["along a 2"]

along_b_4 = absence_classes["along b 4"]
along_b   = absence_classes["along b"]
along_b_2 = absence_classes["along b 2"]

along_c   = absence_classes["along c"]
along_c_2 = absence_classes["along c 2"]
along_c_3 = absence_classes["along c 3"]
along_c_4 = absence_classes["along c 4"]
along_c_6 = absence_classes["along c 6"]

absences_via_intensity_symmetry = {
"P -1"       : []                      ,# l>0 or (l==0 and (h>0 or (h==0 and k>=0)))

"P 1 2/m 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"C 1 2/m 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"A 1 2/m 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"I 1 2/m 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))

"P 1 1 2/m"  : along_c+along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"A 1 1 2/m"  : along_c+along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"B 1 1 2/m"  : along_c+along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"I 1 1 2/m"  : along_c+along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))

"P 2/m 1 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"B 2/m 1 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"C 2/m 1 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"I 2/m 1 1"  : along_b+along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))

"P m m m"    : along_a+along_b+along_c+along_c_2 ,# h>=0 and k>=0 and l>=0
"C m m m"    : along_a+along_b+along_c+along_c_2 ,# h>=0 and k>=0 and l>=0
"A m m m"    : along_a+along_b+along_c+along_c_2 ,# h>=0 and k>=0 and l>=0
"B m m m"    : along_a+along_b+along_c+along_c_2 ,# h>=0 and k>=0 and l>=0
"F m m m"    : along_a+along_b+along_c+along_c_2 ,# h>=0 and k>=0 and l>=0
"I m m m"    : along_a+along_b+along_c+along_c_2 ,# h>=0 and k>=0 and l>=0

"P 4/m"      : along_c_4+along_c ,# l>=0 and ((h>=0 and k>0) or (h==0 and k==0))
"I 4/m"      : along_c_4+along_c         ,# l>=0 and ((h>=0 and k>0) or (h==0 and k==0))

"P 4/m m m"  : along_a+along_c_4+along_c         ,# h>=k and k>=0 and l>=0
"I 4/m m m"  : along_a+along_c_4+along_c         ,# h>=k and k>=0 and l>=0

"P -3"       : along_b+along_c_3+along_c         ,# (h>=0 and k>0) or (h==0 and k==0 and l>=0)
"R -3 :H"    : along_b+along_c_3+along_c        ,# (h>=0 and k>0) or (h==0 and k==0 and l>=0)
"R -3 :R"    : along_b+along_c_3+along_c         ,# (h>=0 and k>0) or (h==0 and k==0 and l>=0)

"P -3 1 m"   : along_b+along_c_3+along_c+along_b_2        ,# h>=k and k>=0 and (k>0 or l>=0)
"P -3 m 1"   : along_b+along_c_3+along_c+along_b_2        ,# h>=k and k>=0 and (h>k or l>=0)
"R -3 m :H"  : along_b+along_c_3+along_c+along_b_2        ,# h>=k and k>=0 and (h>k or l>=0)
"R -3 m :R"  : along_b+along_c_3+along_c+along_b_2        ,# h>=k and k>=0 and (h>k or l>=0)

"P 6/m"      : along_c+along_b+along_c_6         ,# l>=0 and ((h>=0 and k>0) or (h==0 and k==0))
"P 6/m m m"  : along_c+along_b+along_c_6         ,# h>=k and k>=0 and l>=0

"P m -3"     : along_c+along_c_2+along_c_4+along_b_2                ,# h>=0 and ((l>=h and k>h) or (l==h and k==h))
"F m -3"     : along_c+along_c_2+along_c_4+along_b_2                 ,# h>=0 and ((l>=h and k>h) or (l==h and k==h))
"I m -3"     : along_c+along_c_2+along_c_4+along_b_2                 ,# h>=0 and ((l>=h and k>h) or (l==h and k==h))

"P m -3 m"   : along_b+along_b_4+along_b_2                 ,# k>=l and l>=h and h>=0
"F m -3 m"   : along_b+along_b_4+along_b_2                 ,# k>=l and l>=h and h>=0
"I m -3 m"   : along_b+along_b_4+along_b_2                 ,# k>=l and l>=h and h>=0
}


absences_via_intensity_symmetry_prot = {
"P -1"       : []                      ,# l>0 or (l==0 and (h>0 or (h==0 and k>=0)))

"P 1 2/m 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"C 1 2/m 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"A 1 2/m 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"I 1 2/m 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))

"P 1 1 2/m"  : along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"A 1 1 2/m"  : along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"B 1 1 2/m"  : along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"I 1 1 2/m"  : along_c_2       ,# k>=0 and (l>0 or (l==0 and h>=0))

"P 2/m 1 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"B 2/m 1 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"C 2/m 1 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))
"I 2/m 1 1"  : along_b_2       ,# k>=0 and (l>0 or (l==0 and h>=0))

"P m m m"    : along_a_2+along_b_2+along_c_2 ,# h>=0 and k>=0 and l>=0
"C m m m"    : along_c_2 ,# h>=0 and k>=0 and l>=0
"A m m m"    : along_a_2 ,# h>=0 and k>=0 and l>=0
"B m m m"    : along_b_2 ,# h>=0 and k>=0 and l>=0
"F m m m"    : along_a_2+along_b_2+along_c_2 ,# h>=0 and k>=0 and l>=0
"I m m m"    : along_a_2+along_b_2+along_c_2 ,# h>=0 and k>=0 and l>=0

"P 4/m"      : along_c_4 ,# l>=0 and ((h>=0 and k>0) or (h==0 and k==0))
"I 4/m"      : along_c_4         ,# l>=0 and ((h>=0 and k>0) or (h==0 and k==0))

"P 4/m m m"  : along_a_2+along_c_4         ,# h>=k and k>=0 and l>=0
"I 4/m m m"  : along_a_2+along_c_4         ,# h>=k and k>=0 and l>=0

"P -3"       : along_c_3        ,# (h>=0 and k>0) or (h==0 and k==0 and l>=0)
"R -3 :H"    : along_c_3        ,# (h>=0 and k>0) or (h==0 and k==0 and l>=0)
"R -3 :R"    : along_c_3         ,# (h>=0 and k>0) or (h==0 and k==0 and l>=0)

"P -3 1 m"   : along_c_3        ,# h>=k and k>=0 and (k>0 or l>=0)
"P -3 m 1"   : along_c_3        ,# h>=k and k>=0 and (h>k or l>=0)
"R -3 m :H"  : along_c_3        ,# h>=k and k>=0 and (h>k or l>=0)
"R -3 m :R"  : along_c_3        ,# h>=k and k>=0 and (h>k or l>=0)

"P 6/m"      : along_c_6         ,# l>=0 and ((h>=0 and k>0) or (h==0 and k==0))
"P 6/m m m"  : along_c_6         ,# h>=k and k>=0 and l>=0

"P m -3"     : along_c_2+along_c_4+along_b_2                ,# h>=0 and ((l>=h and k>h) or (l==h and k==h))
"F m -3"     : along_c_2+along_c_4+along_b_2                 ,# h>=0 and ((l>=h and k>h) or (l==h and k==h))
"I m -3"     : along_c_2+along_c_4+along_b_2                 ,# h>=0 and ((l>=h and k>h) or (l==h and k==h))

"P m -3 m"   : along_b_4+along_b_2                 ,# k>=l and l>=h and h>=0
"F m -3 m"   : along_b_4+along_b_2                 ,# k>=l and l>=h and h>=0
"I m -3 m"   : along_b_4+along_b_2                 ,# k>=l and l>=h and h>=0
}


lattice_from_string = { "1 (0, 0, 0) 0,1/2,1/2"   : "A",
                        "1 (0, 0, 0) 1/2,0,1/2"   : "B",
                        "1 (0, 0, 0) 1/2,1/2,0"   : "C",
                        "1 (0, 0, 0) 1/2,1/2,1/2" : "I",
                      }


symbol_from_string = { "-2 (1, 0, 0) 0,1/2,0" : "b (a)",
                       "-2 (1, 0, 0) 0,0,1/2" : "c (a)",
                       "-2 (1, 0, 0) 0,1/2,1/2" : "n (a)",
                       "2 (1, 0, 0) 0,0,0" : "2_0 (a)",     # 2a
                       "2 (1, 0, 0) 1/2,0,0" : "2_1 (a)",
                       "-2 (0, 0, 1) 1/2,1/2,0" : "n (c)",
                       "6 (0, 0, 1) 0,0,0" : "6_0 (c)",     # 6c
                       "6 (0, 0, 1) 0,0,2/3" : "6_4 (c)",
                       "-2 (1, 0, 0) 0,1/4,1/4" : "d (a)",
                       "6 (0, 0, 1) 0,0,5/6" : "6_5 (c)",
                       "6 (0, 0, 1) 0,0,1/2" : "6_3 (c)",
                       "-2 (0, 1, 0) 1/2,0,0" : "a (b)",
                       "2 (0, 0, 1) 0,0,0"  : "2_0 (c)",    #2c
                       "2 (0, 0, 1) 0,0,1/2" : "2_1 (c)",
                       "4 (0, 0, 1) 0,0,0" : "4_0 (c)",     #4c
                       "4 (0, 0, 1) 0,0,1/4" : "4_1 (c)",
                       "-2 (0, 0, 1) 1/4,1/4,0" : "d (c)",
                       "4 (1, 0, 0) 1/4,0,0" : "4_1 (a)",
                       "4 (1, 0, 0) 0,0,0" : "4_0 (a)",     #4a
                       "3 (0, 0, 1) 0,0,0" : "3_0 (c)",     #3c
                       "3 (0, 0, 1) 0,0,1/3" : "3_1 (c)",
                       "4 (0, 0, 1) 0,0,3/4" : "4_3 (c)",
                       "-2 (0, 0, 1) 0,1/2,0" : "b (c)",
                       "6 (0, 0, 1) 0,0,1/3" : "6_2 (c)",
                       "-2 (0, 1, 0) 1/2,0,1/2" : "n (b)",
                       "4 (1, 0, 0) 3/4,0,0" : "4_3 (a)",
                       "-2 (0, 1, 0) 0,0,1/2" : "c (b)",
                       "6 (0, 0, 1) 0,0,1/6" : "6_1 (c)",
                       "-2 (0, 0, 1) 1/2,0,0" : "a (c)",
                       "2 (0, 1, 0) 0,0,0" : "2_0 (b)",     #2b
                       "2 (0, 1, 0) 0,1/2,0" : "2_1 (b)",
                       "-2 (0, 1, 0) 1/4,0,1/4" : "d (b)",
                       "4 (0, 0, 1) 0,0,1/2" : "4_2 (c)",
                       "4 (1, 0, 0) 2/4,0,0" : "4_2 (a)",
                       "3 (0, 0, 1) 0,0,2/3" : "3_2 (c)" }

string_from_symbol ={"c (b)"  :  "-2 (0, 1, 0) 0,0,1/2",
                     "n (b)"  :  "-2 (0, 1, 0) 1/2,0,1/2",
                     "a (b)"  :  "-2 (0, 1, 0) 1/2,0,0",
                     "a (c)"  :  "-2 (0, 0, 1) 1/2,0,0",
                     "n (c)"  :  "-2 (0, 0, 1) 1/2,1/2,0",
                     "b (c)"  :  "-2 (0, 0, 1) 0,1/2,0",
                     "b (a)"  :  "-2 (1, 0, 0) 0,1/2,0",
                     "n (a)"  :  "-2 (1, 0, 0) 0,1/2,1/2",
                     "c (a)"  :  "-2 (1, 0, 0) 0,0,1/2",
                     "2_1 (b)":   "2 (0, 1, 0) 0,1/2,0",
                     "2_1 (c)":   "2 (0, 0, 1) 0,0,1/2",
                     "2_1 (a)":   "2 (1, 0, 0) 1/2,0,0",
                     "2_0 (b)":   "2 (0, 1, 0) 0,0,0",
                     "2_0 (c)":   "2 (0, 0, 1) 0,0,0",
                     "2_0 (a)":   "2 (1, 0, 0) 0,0,0",
                     "d (b)"  :  "-2 (0, 1, 0) 1/4,0,1/4",
                     "d (c)"  :  "-2 (0, 0, 1) 1/4,1/4,0",
                     "d (a)"  :  "-2 (1, 0, 0) 0,1/4,1/4",
                     "4_0 (c)":   "4 (0, 0, 1) 0,0,0",
                     "4_1 (c)":   "4 (0, 0, 1) 0,0,1/4",
                     "4_2 (c)":   "4 (0, 0, 1) 0,0,1/2",
                     "4_3 (c)":   "4 (0, 0, 1) 0,0,3/4",
                     "4_0 (a)":   "4 (1, 0, 0) 0,0,0",
                     "4_1 (a)":   "4 (1, 0, 0) 1/4,0,0",
                     "4_2 (a)":   "4 (1, 0, 0) 2/4,0,0",
                     "4_3 (a)":   "4 (1, 0, 0) 3/4,0,0",
                     "3_0 (c)":   "3 (0, 0, 1) 0,0,0",
                     "3_1 (c)":   "3 (0, 0, 1) 0,0,1/3",
                     "3_2 (c)":   "3 (0, 0, 1) 0,0,2/3",
                     "6_0 (c)":   "6 (0, 0, 1) 0,0,0",
                     "6_1 (c)":   "6 (0, 0, 1) 0,0,1/6",
                     "6_5 (c)":   "6 (0, 0, 1) 0,0,5/6",
                     "6_4 (c)":   "6 (0, 0, 1) 0,0,2/3",
                     "6_2 (c)":   "6 (0, 0, 1) 0,0,1/3",
                     "6_3 (c)":   "6 (0, 0, 1) 0,0,1/2" }

equivs = { "6_1 (c)":"6_5 (c)", "6_2 (c)":"6_4 (c)",
           "6_5 (c)":"6_1 (c)", "6_4 (c)":"6_2 (c)",
           "4_1 (c)":"4_3 (c)",
           "4_3 (c)":"4_1 (c)",
           "3_1 (c)":"3_2 (c)",
           "3_2 (c)":"3_1 (c)" }


class conditions_for_operator(object):
  def __init__(self, s ):
    r_info = sgtbx.rot_mx_info( s.r() )
    t_info = sgtbx.translation_part_info( s )

    self.type = r_info.type()
    self.ev   = r_info.ev()
    self.trans= t_info.intrinsic_part()

  def condition(self):
    cond = []
    # get conditions for hkl
    if self.type < 0:
      if ( list(self.trans.as_double()) ).count(0) <3 :
        for ii in ev:
          if ii == 0:
            cond.append( None )
          else:
            cond.append( 0 )
      else:
        cond = (None,None,None)
    if type > 0:
      if ( list(self.trans.as_double()) ).count(0) <3 :
        for ii in ev:
          if ii != 0:
            cond.append(None)
          else:
            cond.append( 0 )

    if len(cond) == 0:
      cond = [ None,None,None ]

    return [ cond, list(t) ]

  def absence_type(self):
    #get conditions for translations
    t = self.trans.as_double()
    id_string = str(self.type)+" "+str(self.ev)+" "+str(self.trans)
    op_name = "Non absent"
    if id_string in symbol_from_string:
      op_name = symbol_from_string[ id_string ]
    return op_name

class absences(object):
  def __init__(self, mult=2.0, threshold=0.95):
    self.lib = absence_and_conditions
    self.absence_classes = absences_via_intensity_symmetry_prot
    self.mult = mult
    self.threshold = threshold

  def check(self, abs_type, hkl, return_bool=False ):
    # check if it is there
    message = "The reflection with index %s is "%(str(hkl))
    if abs_type in self.lib:
      mask, condition = self.lib[ abs_type ]
      mc = self.check_mask( hkl, mask )
      cc = self.check_condition(hkl, condition )
      if mc:
        if cc:
          message += " PRESENT under the %s operator"%(abs_type)
        else:
          message += " ABSENT under the %s operartor"%(abs_type)
      else:
        message += " not of the type that is absent under %s"%(abs_type)
    else:
      message =  "No such type of absence"
    if return_bool:
      return mc,cc
    else:
      return message

  def check_mask(self, hkl, mask):
    mask_none_count =0
    for ii in mask:
      if ii is None:
        mask_none_count += 1
    count = 0
    result = False
    for m,h in zip(mask,hkl):
      if m == h:
        count += 1
    if count == 3-mask_none_count:
      result = True
    return result

  def check_condition(self, hkl, condition):
    """magic"""
    result = False
    x = hkl[0]*condition[0] + hkl[1]*condition[1] + hkl[2]*condition[2]
    y = math.cos( x*2.0*math.pi )
    z = 0.5*(math.tanh(self.mult*y)+1)
    if z > self.threshold:
      result = True
    return result

def likelihood(z,sigz,absent_or_centric_or_acentric,sigma_inflation=1.0):
  from mmtbx.scaling import absence_likelihood
  flag = absent_or_centric_or_acentric
  result = absence_likelihood.log_p( z,sigz*sigma_inflation,flag )
  return result


class analyze_absences(mmtbx.scaling.xtriage_analysis):
  def __init__(self, miller_array, isigi_cut=3, sigma_inflation=1.0):
    self.cut = isigi_cut
    self.sigma_inflation=sigma_inflation
    self.miller_array = miller_array.deep_copy()

    self.n_abs        = []
    self.n_n_abs      = []
    self.n_tot        = []
    self.n_abs_viol   = []
    self.n_n_abs_viol = []
    self.n_tot_viol   = []

    self.isi_abs      = []
    self.isi_n_abs    = []
    self.isi_tot      = []

    self.i_abs        = []
    self.i_n_abs      = []
    self.i_tot        = []

    self.op_name      = []
    self.score        = []
    self.present      = []

    assert self.miller_array.sigmas() is not None
    #we need to have this in the standard setting please
    self.sg = sgtbx.space_group_info( group = self.miller_array.space_group() )
    self.cb_op = self.sg.change_of_basis_op_to_reference_setting()
    self.miller_array = self.miller_array.change_basis( self.cb_op )
    self.abs_check = absences()
    self.check_conditions()

  def _show_impl(self, out):
    out.show_sub_header("Table of systematic absence rules")
    out.show("""\
 The following table gives information about systematic absences allowed for
 the specified intensity point group.

 For each operator, the reflections are split in three classes:
""")
    out.show_preformatted_text("""
  Systematic absence: Reflections that are absent for this operator.
  Non absence       : Reflections of the same type (i.e. (0,0,l)) as above, but they
                      should be present.
  Other reflections : All other reflections.
""")
    out.show("""\
For each class, the <I/sigI> is reported, as well as the number of
violations. A violation is a reflection that is absent when it is expected
to be present for a particular space group, or present when it is
expected to be absent. The criteria are:""")
    out.show_preformatted_text("""
  Systematic absence violation: I/sigI > %(cut)2.1f
  Non absence violation       : I/sigI < %(cut)2.1f
  Other relections violation  : I/sigI < %(cut)2.1f
""" % {"cut":self.cut})
    out.show_text("""\
Operators with low associated violations for *both* systematically absent and
non absent reflections, are likely to be true screw axis or glide planes. Both
the number of violations and their percentages are given.  The number of
violations within the 'other reflections' class, can be used as a comparison
for the number of violations in the non-absent class.
""")
    out.show_table(self.table)

  def score_isigi(self,isig, absent=False, a=30.0):
    tmp = 0.5*(1+math.tanh( a*(isig-self.cut) ) )
    if not absent:
      return abs(math.log(tmp+1e-8))
    else:
      return abs(math.log(1-tmp+1e-8) )

  def propose(self, ops, thres=1):
    in_sg_and_seemingly_correct       = []
    not_in_sg_and_seemingly_incorrect = []
    in_sg_and_seemingly_incorrect     = []
    observed_but_not_in_sg            = []
    in_sg_but_no_observations         = []

    absent_class_violations     = 0
    not_absent_class_violations = 0

    total_score = 0

    for op, sc, n, n_n, n_v, n_n_v  in zip(self.op_name, self.score,
        self.n_abs, self.n_n_abs, self.n_abs_viol, self.n_n_abs_viol):
      if op in ops:
        total_score += sc
        absent_class_violations += n_v
        not_absent_class_violations += n_n_v

      if n > 0: # observed
        if sc <= thres: # correct
          if op in ops: # in proposed sg
            in_sg_and_seemingly_correct.append( op )
          else: # correct but not in proposed sg
            observed_but_not_in_sg.append( op )
        else: # not correct
          if op not in ops: # not in sg, and not correct
            not_in_sg_and_seemingly_incorrect.append( op )
      else: # not observed
        if op in ops: #in proposed sg
          in_sg_but_no_observations.append( op )

    pos = len(in_sg_and_seemingly_correct) + \
          len(not_in_sg_and_seemingly_incorrect)
    neg = len(in_sg_and_seemingly_incorrect) + len(observed_but_not_in_sg)
    abstain = len(in_sg_but_no_observations)

    return total_score,absent_class_violations, not_absent_class_violations

  def check_conditions(self,abs_lower_i_threshold=1e-6):
    table_labels = ('Operator',
      "# expected systematic absences",
      "<I/sigI> (violations)",
      "# expected non absences",
      "<I/sigI> (violations)",
      "# other reflections",
      "<I/sigI> (violations)",
      "Score")
    for  item in [0]: # absence_class in self.abs_check.absence_classes[ self.sg.group().crystal_system() ]:
      table_rows = []
      for condition in self.abs_check.absence_classes[
        str(sgtbx.space_group_info(
          group=self.sg.group().build_derived_reflection_intensity_group(False))\
            .as_reference_setting()
            ) ] : # crystal_system() ]:
        n_abs   = 0
        n_n_abs = 0
        n_tot   = 0
        n_abs_viol   = 0
        n_n_abs_viol = 0
        n_tot_viol   = 0

        isi_abs     = 0
        isi_n_abs   = 0
        isi_tot = 0

        i_abs     = 0
        i_n_abs   = 0
        i_tot  = 0

        score = 0

        for hkl, centric_flag, i, sigi in zip(self.miller_array.indices(), self.miller_array.centric_flags(), self.miller_array.data(), self.miller_array.sigmas() ):
          mc, cc = self.abs_check.check(condition,hkl, return_bool=True)
          if abs(i) < abs_lower_i_threshold:
            sigi=max(sigi,abs_lower_i_threshold)
          if mc: # mask checks out
            if cc: # not absent
              n_n_abs += 1
              isi_n_abs += i/sigi
              i_n_abs   += i
              # should be present. flag if not significant
              if i/sigi < self.cut:
                n_n_abs_viol += 1
              score += likelihood(i,sigi,centric_flag[1],self.sigma_inflation)
            else: #absent
              n_abs += 1
              isi_abs += i/sigi
              i_abs   += i
              # should be absent: flag if significant
              if i/sigi > self.cut:
                n_abs_viol += 1
              score += likelihood( i,sigi,None)
          else:
            n_tot +=1
            isi_tot += i/sigi
            i_tot += i
            if i/sigi <  self.cut:
              n_tot_viol += 1
        if n_abs > 0:
          isi_abs   = isi_abs/n_abs
          i_abs   = i_abs/n_abs
        if n_n_abs > 0:
          isi_n_abs = isi_n_abs/n_n_abs
          i_n_abs = i_n_abs/n_n_abs
        if n_tot > 0:
          isi_tot   = isi_tot/n_tot
          i_tot   = i_tot/n_tot

        self.n_abs.append(n_abs)
        self.n_n_abs.append(n_n_abs)
        self.n_tot.append(n_tot)
        self.n_abs_viol.append(n_abs_viol)
        self.n_n_abs_viol.append(n_n_abs_viol)
        self.n_tot_viol.append(n_tot_viol)

        self.isi_abs.append(isi_abs)
        self.isi_n_abs.append(isi_n_abs)
        self.isi_tot.append(isi_tot)

        self.i_abs.append(i_abs)
        self.i_n_abs.append(i_n_abs)
        self.i_tot.append(i_tot)

        self.op_name.append( condition )
        score = float(score)/max(1,n_abs+n_n_abs)
        self.score.append( score )

        table_rows.append( [condition,
          str("%8.0f"%(n_abs)),
          str("%8.2f  (%i, %4.1f%%)" % (isi_abs, n_abs_viol,
            100.0*float(n_abs_viol)/max(1,n_abs))),
          str("%8.0f"%(n_n_abs)),
          str("%8.2f  (%i, %4.1f%%)" % (isi_n_abs, n_n_abs_viol,
            100.0*float(n_n_abs_viol)/max(1,n_n_abs))),
          str("%8.0f"%(n_tot)),
          str("%8.2f  (%i, %4.1f%%)" % (isi_tot, n_tot_viol,
            100.0*float(n_tot_viol)/max(1,n_tot))),
          str("%8.2e"%(abs(score)))
        ])
      self.table = table_utils.simple_table(
        column_headers=table_labels,
        table_rows=table_rows)

class sgi_iterator(object):
  def __init__(self,chiral=True, crystal_system=None, intensity_symmetry=None):
    self.chiral = chiral
    self.crystal_system = crystal_system
    self.intensity_symmetry = intensity_symmetry
    assert ( [self.crystal_system , self.intensity_symmetry] ).count(None) != 0

  def comparator(self, sgi):
    if self.crystal_system is not None:
      return ((self.crystal_system is None) or
              (self.crystal_system == sgi.group().crystal_system()))
    else:
      return ((self.intensity_symmetry is None) or (self.intensity_symmetry ==
        sgi.group().build_derived_reflection_intensity_group(False)))

  def list(self):
    for symbols in sgtbx.space_group_symbol_iterator():
      sgi = sgtbx.space_group_info(group=sgtbx.space_group(
        space_group_symbols=symbols))
      if self.comparator(sgi):
        if (self.chiral is None) or (self.chiral == sgi.group().is_chiral()):
          yield sgi


class protein_space_group_choices(mmtbx.scaling.xtriage_analysis):
  def __init__(self,
      miller_array,
      threshold = 3,
      protein=True,
      print_all=True,
      sigma_inflation=1.0,
      original_data=None):
    self.threshold = 3.0
    assert miller_array.is_xray_intensity_array()
    self.miller_array = miller_array.deep_copy().f_sq_as_f(
      ).average_bijvoet_mates().f_as_f_sq().map_to_asu()
    space_group = self.miller_array.space_group()

    self.absences_table = analyze_absences(
      miller_array=self.miller_array,
      isigi_cut=threshold,
      sigma_inflation=sigma_inflation)
    if (original_data is not None):
      self.absences_list = absences_list(obs=original_data,
        was_filtered=False)
    else :
      self.absences_list = absences_list(obs=self.miller_array,
        was_filtered=True)

    self.sg_iterator = sgi_iterator(chiral = True,
      intensity_symmetry = \
        space_group.build_derived_reflection_intensity_group(False) )

    self.sg_choices  = []
    self.mean_i      = []
    self.mean_isigi  = []
    self.n           = []
    self.violations  = []
    self.abs_types   = []
    self.tuple_score = []

    score = []

    for sg in self.sg_iterator.list():
      xs = crystal.symmetry(
        unit_cell = self.miller_array.unit_cell(),
        space_group = sg.group())
      tmp_miller = self.miller_array.customized_copy( crystal_symmetry = xs )
      these_absent_millers = tmp_miller.select(
        tmp_miller.sys_absent_flags().data() )

      if these_absent_millers.data().size() > 0:
        tmp_mean_i = flex.mean( these_absent_millers.data() )
        zero_sel = these_absent_millers.sigmas()==0
        these_absent_millers = these_absent_millers.select(~zero_sel)
        #print sg, list(these_absent_millers.indices()), list(these_absent_millers.data())
        tmp_mean_isigi = flex.mean(
          these_absent_millers.data() / these_absent_millers.sigmas() )
        tmp_n = these_absent_millers.data().size()
        tmp_violations = flex.bool( these_absent_millers.data() /
          these_absent_millers.sigmas() > self.threshold ).count( True )
      else:
        tmp_mean_i = 0
        tmp_mean_isigi = 0
        tmp_n = 0
        tmp_violations = 0

      to_be_checked = []
      for s in sg.group():
        #check if this is an operator that causes absences
        tmp =  conditions_for_operator( s )
        if tmp.absence_type() != "None":
          if tmp.absence_type() in self.absences_table.op_name:
            ii = self.absences_table.op_name.index( tmp.absence_type() )
            if tmp.absence_type() not in to_be_checked:
              if tmp.absence_type() in equivs:
                if equivs[ tmp.absence_type() ] not in to_be_checked:
                  to_be_checked.append( tmp.absence_type() )
                  tmp_score = self.absences_table.score[ ii ]
              else:
                  to_be_checked.append( tmp.absence_type() )
                  tmp_score = self.absences_table.score[ ii ]

      self.abs_types.append( to_be_checked )
      tuple_score =  self.absences_table.propose( to_be_checked )
      self.tuple_score.append( tuple_score )

      self.sg_choices.append(  sg )
      self.mean_i.append( tmp_mean_i )
      self.mean_isigi.append( tmp_mean_isigi )
      self.n.append( tmp_n )
      self.violations.append( tmp_violations )
    tmp_rows = self.suggest_likely_candidates()
    self.sorted_table = table_utils.simple_table(
      column_headers=['space group', '#  absent', '<Z>_absent',
                      '<Z/sigZ>_absent', '+++', '---', 'score'],
      table_rows=tmp_rows)

  def _show_impl(self, out):
    out.show_header("Systematic absences")
    self.absences_table.show(out)
    out.show_sub_header("Space group identification")
    out.show_text("""\
Analyses of the absences table indicates a number of likely space group
candidates, which are listed below. For each space group, the number of
systematic absence violations are listed under the '+++' column. The number of
non-absence violations (weak reflections) are listed under '---'. The last
column is a likelihood based score for the particular space group.  Note that
enantiomorphic spacegroups will have equal scores. Also, if absences were
removed while processing the data, they will be regarded as missing
information, rather then as enforcing that absence in the space group choices.
""")
    out.show_table(self.sorted_table)
    if (getattr(self, "absences_list", None) is not None) : # backwards compat.
      if (self.absences_list.n_possible_max > 0):
        self.absences_list.show(out)

  def suggest_likely_candidates( self, acceptable_violations = 1e+90 ):
    used = flex.bool( len(self.sg_choices), False )
    order = []

    all_done = False
    count = -1
    if (len(self.tuple_score) == 0):
      return []
    tmp_scores = []
    for tt in self.tuple_score:
      tmp_scores.append( tt[0] )
    order = flex.sort_permutation( flex.double( tmp_scores ), False  )


    sorted_rows = []
    max_score = flex.min( flex.double( tmp_scores ) )
    for ii in order:
      sg             = self.sg_choices[ii]
      tmp_n          = self.n[ii]
      tmp_violations = self.violations[ii]
      tmp_mean_i     = self.mean_i[ii]
      tmp_mean_isigi = self.mean_isigi[ii]
      tuple_score    = self.tuple_score[ii]

      sorted_rows.append( [str(sg), '%i'%(tmp_n),
                           '%8.2f  '%(tmp_mean_i),
                           '%8.2f  '%(tmp_mean_isigi),
                           ' %i '%(tuple_score[1]),
                           ' %i '%(tuple_score[2]),
                           ' %8.3e '%((tuple_score[0]-max_score))
                          ])

    return sorted_rows

class absences_list(mmtbx.scaling.xtriage_analysis,
                     miller.systematic_absences_info):
  """
  Container for lists of systematic absences.  This subclass simply overrides
  the default output of the base class in cctbx.miller to be consistent with
  the rest of Xtriage.
  """
  def show(self, *args, **kwds):
    mmtbx.scaling.xtriage_analysis.show(self, *args, **kwds)

  def _show_impl(self, out):
    """
    For each possible space group, show a list of possible systematically
    absent reflections and corresponding I/sigmaI.
    """
    out.show_sub_header("List of individual systematic absences")
    out.show_text("""\
 Note: this analysis uses the original input data rather than the filtered data
 used for twinning detection; therefore, the results shown here may include
 more reflections than shown above.
""")
    if (self.input_amplitudes):
      out.show_text("""\
 Also note that the input data were amplitudes, which means that weaker
 reflections may have been modified by French-Wilson treatment or discarded
 altogether, and the original intensities will not be recovered.
""")
    for group_info, absences in self.space_group_symbols_and_selections :
      group_note = ""
      if (str(group_info) == str(self.space_group_info)):
        group_note = " (input space group)"
      if (absences == False):
        out.show_paragraph_header("%s%s: no systematic absences possible" % \
          (group_info, group_note))
      elif (absences is None):
        out.show_paragraph_header("%s%s: no absences found" % \
          (group_info, group_note))
      else :
        out.show_paragraph_header("%s%s" % (group_info, group_note))
        lines = []
        for i_hkl, hkl in enumerate(absences.indices()):
          intensity = absences.data()[i_hkl]
          sigma = absences.sigmas()[i_hkl]
          indices_fmt = "(%4d, %4d, %4d)" % hkl
          if (sigma == 0):
            lines.append("  %s: i/sigi = undefined" % indices_fmt)
          else :
            lines.append("  %s: i/sigi = %6.1f" % (indices_fmt,
              intensity/sigma))
        out.show_preformatted_text("\n".join(lines))

def test():
  tmp = absences()
  assert tmp.check( "2_1 (c)", (0,0,1),True ) == (True, False)
  assert tmp.check( "2_1 (c)", (0,0,4),True ) == (True, True)
  assert tmp.check( "4_1 (a)", (4,0,0),True ) == (True, True)
  assert tmp.check( "3_1 (c)", (0,0,3),True ) == (True, True)

  tmp = sgi_iterator(chiral = True,
    crystal_system = None,
    intensity_symmetry = sgtbx.space_group_info( "P222").group().build_derived_reflection_intensity_group(False)  )
  sg_list = []
  abs_list = []
  for sg in tmp.list():
    sg_list.append( str(sg) )
    if str(sg)== "P 21 21 21":
      for s in sg.group():
        abs_list.append( conditions_for_operator( s ).absence_type() )
  assert "2_1 (a)" in abs_list
  assert "2_1 (b)" in abs_list
  assert "2_1 (c)" in abs_list


  assert "P 2 2 2" in sg_list
  assert "P 21 2 2" in sg_list
  assert "P 2 21 2" in sg_list
  assert "P 2 2 21" in sg_list
  assert "P 21 21 2" in sg_list
  assert "P 21 2 21" in sg_list
  assert "P 2 21 21" in sg_list
  assert "P 21 21 21" in sg_list

  """
  tmp = sgi_iterator(chiral=None)
  pg  = []
  asu = []
  for sgi in tmp.list():
    a = str( sgtbx.space_group_info( group = sgi.group().build_derived_reflection_intensity_group(False) ) )
    b = str( sgi.reciprocal_space_asu().reference_as_string() )
    if a not in pg:
      pg.append( a )
      asu.append( b )
  for i,j in zip(pg,asu):
    print i, "   ", j
  """

  print("OK")

if __name__ == "__main__":
  test()
