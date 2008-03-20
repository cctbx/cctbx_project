from cctbx import sgtbx, uctbx
from cctbx import miller
from cctbx.array_family import flex
from cctbx import sgtbx
from libtbx import table_utils
import math,sys



#                           name     hkl selector  condition
absence_and_conditions = { "2_1 (a)" : [(None,0,0),  (1.0/2.0,0,0)],
                           "2_1 (b)" : [(0,None,0),  (0,1.0/2.0,0)],
                           "2_1 (c)" : [(0,0,None),  (0,0,1.0/2.0)],

                           "3_1 (c)" : [(0,0,None),  (0,0,1.0/3.0)],
                           "3_2 (c)" : [(0,0,None),  (0,0,1.0/3.0)],

                           "4_1 (c)" : [(0,0,None),  (0,0,1.0/4.0)],
                           "4_2 (c)" : [(0,0,None),  (0,0,2.0/4.0)],
                           "4_3 (c)" : [(0,0,None),  (0,0,3.0/4.0)],

                           "4_1 (a)" : [(None,0,0),  (1.0/4.0,0,0)],
                           "4_2 (a)" : [(None,0,0),  (2.0/4.0,0,0)],
                           "4_3 (a)" : [(None,0,0),  (3.0/4.0,0,0)],

                           "4_1 (b)" : [(0,None,0),  (0,1.0/4.0,0)],
                           "4_2 (b)" : [(0,None,0),  (0,2.0/4.0,0)],
                           "4_3 (b)" : [(0,None,0),  (0,3.0/4.0,0)],

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

absence_classes = { "along a" : ["2_1 (a)", "b (a)", "c (a)", "n (a)", "d (a)"],
                    "along b" : ["2_1 (b)", "a (b)", "c (b)", "n (b)", "d (b)"],
                    "along c" : ["2_1 (c)", "3_1 (c)", "3_2 (c)",
                                 "4_1 (c)", "4_2 (c)", "4_3 (c)",
                                 "6_1 (c)", "6_2 (c)", "6_3 (c)", "6_4 (c)", "6_5 (c)",
                                 "a (c)", "b (c)", "n (c)", "d (c)"] }


absence_classes_via_system_protein = { "Triclinic"   : [],
                                       "Monoclinic"  : ["2_1 (b)"],
                                       "Orthorohmbic": ["2_1 (a)","2_1 (b)","2_1 (c)", ],
                                       "Tetragonal"  : ["2_1 (b)","4_1 (c)","4_2 (c)","4_3 (c)"],
                                       "Trigonal"    : ["3_1 (c)", "3_2 (c)"],
                                       "Hexagonal"   : ["6_1 (c)", "6_2 (c)", "6_3 (c)", "6_4 (c)", "6_5 (c)"],
                                       "Cubic"       : ["2_1 (a)", "4_1 (a)", "4_2 (a)", "4_3(a)"]
                                      }


class absences(object):
  def __init__(self, mult=2.0, threshold=0.95,protein=True):
    self.lib = absence_and_conditions
    self.absence_classes = absence_classes
    if protein:
      self.absence_classes = absence_classes_via_system_protein
    self.mult = mult
    self.threshold = threshold

  def check(self, abs_type, hkl, return_bool=False ):
    # check if it is there
    message = "The reflection with index %s is "%(str(hkl))
    if self.lib.has_key( abs_type ):
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



class analyze_absences(object):
  def __init__(self, miller_array, isigi_cut=3, out=None):
    if out is None:
      out = sys.stdout
    self.out = out
    self.cut = isigi_cut
    self.miller_array = miller_array.deep_copy()
    assert self.miller_array.sigmas() is not None
    #we need to have this in the standard setting please
    self.sg = sgtbx.space_group_info( group = self.miller_array.space_group() )
    self.cb_op = self.sg.change_of_basis_op_to_reference_setting()
    self.miller_array = self.miller_array.change_basis( self.cb_op )
    self.abs_check = absences()
    self.check_conditions()


  def check_conditions(self):
    table_text = """


Systematic absences
-------------------

The following table gives informaton about systematic absences.
For each operator, the reflections are split in three classes:
  Absent    : Reflections that are absent for this operator.
  Non Absent: Reflection of the same type (i.e. (0,0,l)) as above, but they should be present.
  Complement: All other reflections.
For each class, the <I/sigI> is reported, as well as the number of 'violations'. A 'violation'
is designated as a reflection for which a I/sigI criterion is not met. The criteria are

  Absent violation     : I/sigI > %2.1f
  Non Absent violation : I/sigI < %2.1f
  Complement violation : I/sigI < %2.1f

Operators with low associated violations for *both* absent and non absent reflections, are likely to
be true screw axis or glide planes. Both the number of violations and their percentages are given.
The number of violations within the 'complement' class, can be used as a comparison for the number
of violations in the non-absent class.
"""%(self.cut, self.cut, self.cut)



    table_labels = ('Operator', 'absent under operator\n <I/sigI> (violations)','\nn absent',
                                'not absent under operator \n <I/sigI> (violations)', '\nn not absent',
                                'all other reflections \n <I/sigI> (violations)', '\nn compl' )
    for  item in [0]: # absence_class in self.abs_check.absence_classes[ self.sg.group().crystal_system() ]:
      table_rows = []
      for condition in self.abs_check.absence_classes[  self.sg.group().crystal_system() ]:
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

        for hkl, i, sigi in zip(self.miller_array.indices(), self.miller_array.data(), self.miller_array.sigmas() ):
          mc, cc = self.abs_check.check(condition,hkl, return_bool=True)
          if mc: # mask checks out
            if cc: # not absent
              n_n_abs += 1
              isi_n_abs += i/sigi
              i_n_abs   += i
              # should be present. flag if not significant
              if i/sigi < self.cut:
                n_n_abs_viol += 1
            else: #absent
              n_abs += 1
              isi_abs += i/sigi
              i_abs   += i
              # should be absent: flag if significant
              if i/sigi > self.cut:
                n_abs_viol += 1

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

        table_rows.append( [condition, str("%8.2f  (%i, %4.1f%s)"%(isi_abs,n_abs_viol, 100.0*float(n_abs_viol)/max(1,n_abs),'%' )),
                                       str("%8.0f"%(n_abs)),
                                       str("%8.2f  (%i, %4.1f%s)"%(isi_n_abs,n_n_abs_viol, 100.0*float(n_n_abs_viol)/max(1,n_n_abs),'%' )),
                                       str("%8.0f"%(n_n_abs)),
                                       str("%8.2f  (%i, %4.1f%s)"%(isi_tot,n_tot_viol, 100.0*float(n_tot_viol)/max(1,n_tot),'%' )),
                                       str("%8.0f"%(n_tot))] )
      print >> self.out, table_text
      self.table = table_utils.format([table_labels]+table_rows,
                                       comments=None,
                                       has_header=True,
                                       separate_rows=False,
                                       prefix='| ',
                                       postfix=' |')
      print >> self.out, self.table
      print >> self.out
      print >> self.out



class space_group_choices(object):
  def __init__(self, miller_array):
    print


def test():
  tmp = absences()
  assert tmp.check( "2_1 (c)", (0,0,1),True ) == (True, False)
  assert tmp.check( "2_1 (c)", (0,0,4),True ) == (True, True)
  assert tmp.check( "4_1 (a)", (4,0,0),True ) == (True, True)
  assert tmp.check( "3_1 (c)", (0,0,3),True ) == (True, True)
  print "OK"

if __name__ == "__main__":
  test()
