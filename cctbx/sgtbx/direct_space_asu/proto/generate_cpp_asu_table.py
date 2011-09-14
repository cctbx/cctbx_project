from scitbx.source_generators.utils import join_open
import cctbx
from cctbx.sgtbx.direct_space_asu import reference_table
from cctbx.sgtbx.direct_space_asu.cut_plane import cut


head1 = """\
///
/// Generated code. DO NOT EDIT
///
/// Generator: cctbx/sgtbx/direct_space_asu/proto/generate_cpp_asu_table.py
///
/// Dependencies: cctbx/sgtbx/direct_space_asu/reference_table.py
///
"""

head2 = """\
#include "reference_table.h"

namespace cctbx { namespace sgtbx { namespace asu {

namespace {

typedef sg_vec3 vvv;
"""

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
# direct_space_asu.cut_plane.cut.base_symbol = base_symbol_cpp
# cut_plane.cut.base_symbol = base_symbol_cpp

def out_cpp(asu, f):
  i=0
  for cut in asu.cuts:
    if( i!=0 ):
      print >>f, "    &", cut
    else:
      print >>f, "     ", cut
    i = i + 1


def show_cpp(sg, f):
  asu = reference_table.get_asu(sg)
  func = "asu_%03d"%sg
  print >>f, "facet_collection::pointer ", func, "()   //  Hall: ", asu.hall_symbol
  print >>f, "{\n  return facet_collection_asu("
  out_cpp(asu, f)
  print >>f, "  );\n}\n"
  return func


def make_md5(path):
  import md5
  import libtbx
  import libtbx.load_env
  real_path = libtbx.env.under_dist( "cctbx", path )
  m = md5.new()
  m.update("\n".join(open(real_path).read().splitlines()))
  return "// " + path + ' ' + m.hexdigest()

def run(dr):
  f = join_open(dr, "reference_table.cpp", "w")
  s1md5 = make_md5("sgtbx/direct_space_asu/reference_table.py")
  s2md5 = make_md5("sgtbx/direct_space_asu/proto/generate_cpp_asu_table.py")
  s3md5 = make_md5("sgtbx/direct_space_asu/short_cuts.py")
  print s1md5, '\n', s2md5, '\n', s3md5
  print >>f, head1
  # comma in print adds one whitespace between parameters
  # so it is better to use + to concatanete strings
  # to avoid trailing white spaces
  print >>f, "////////////////\n" + s1md5 + "\n" + s2md5 + "\n" + s3md5 + "\n////////////////\n"
  print >>f, head2
  table = "asu_func asu_table[230] = {"
  i = 0
  for sg in xrange(1,231):
    if( i%8 == 0 ):
      table += "\n "
    func = show_cpp(sg, f)
    table += " "
    table += func
    if i<229 :
      table += ","
    i += 1
  print >>f, "} // end of unnamed namespace\n\n" + table + "\n};\n\n}}}\n"


if (__name__ == "__main__"):
  run(".")

