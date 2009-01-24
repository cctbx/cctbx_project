import cctbx
from cctbx.sgtbx import direct_space_asu
from cctbx.sgtbx.direct_space_asu import reference_table
from boost import rational
from cctbx.sgtbx.direct_space_asu.cut_plane import cut



# replacement for cut.base_symbol
def base_symbol_cpp(thecut):
  from cctbx.sgtbx.direct_space_asu import short_cuts
  n = tuple(thecut.n)
  minus_n = tuple([-e for e in n])
  matching_n = None
  for key,value in short_cuts.__dict__.items():
    if (isinstance(value, cut)):
      if (value.n == n):
        if (value.c == thecut.c):
          return key
        elif (value.c == -thecut.c):
          return "-~"+key
        elif (value.c == 1):
          assert matching_n is None
          if (thecut.c < 0):
            matching_n = "-~"+key
          else:
            matching_n = key
      elif (value.n == minus_n):
        if (value.c == -thecut.c):
          return "-"+key
        elif (value.c == thecut.c):
          return "~"+key
        elif (value.c == 1):
          assert matching_n is None
          if (thecut.c < 0):
            matching_n = "-"+key
          else:
            matching_n = "~"+key
  c = short_cuts.r1 * thecut.c
  num = c.numerator()
  abs_num = abs(num)
  den = c.denominator()
  if (matching_n is None):
    if (num == 0):
      s = "0"
    else:
      s = "r1"
      if (num < 0): s = "-"+s
      if (abs_num != 1): s += "*"+str(abs_num)
      if (den != 1): s += "/"+str(den)
    return "cut(vvv" + str(thecut.n).replace(" ", "") + "," + s + ")"
  s = matching_n
  assert num != 0
  if (abs_num != 1): s += "*"+str(abs_num)
  if (den != 1): s += "/"+str(den)
  return s


cctbx.sgtbx.direct_space_asu.cut_plane.cut.base_symbol = base_symbol_cpp

def out_cpp(asu):
  i=0
  for facet in asu.facets:
    if( i!=0 ):
      print "    &", facet
    else:
      print "     ", facet
    i = i + 1



def show_cpp(sg):
  asu = reference_table.get_asu(sg)
  print "/////   Hall: ", asu.hall_symbol
  func = "asu_%03d"%sg
  print "abstract::ptr ", func, "()"
  print "{"
  print "  return abstract_asu("
  out_cpp(asu)
  print "  );"
  print "}"
  print
  return func


def run():
  print
  print "namespace cctbx { namespace sgtbx { namespace asu {"
  print
  print "typedef ivector3_t vvv;"
  print
  str = "asu_func asu_table[230] = {"
  i = 0
  for sg in xrange(1,231):
    if( i%8 == 0 ):
      str += "\n  "
    func = show_cpp(sg)
    str += func + ", "
    i += 1
  print 
  print ""
  print str
  print "};"
  print
  print "}}}"
  print


if (__name__ == "__main__"):
  run()

