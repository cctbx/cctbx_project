from scitbx.source_generators.utils import join_open
from scitbx.source_generators.utils import write_this_is_auto_generated
import libtbx.load_env
import string
import os

this = "cctbx.source_generators.eltbx.generate_henke_cpp"

reference_tables_directory = libtbx.env.under_dist(
  "cctbx", "reference/henke/tables")

def print_header(f):
  write_this_is_auto_generated(f, this)
  print >> f, """\
#include <cctbx/eltbx/henke.h>

namespace cctbx { namespace eltbx { namespace henke {
"""

def print_ftp_info(f):
  print >> f, """\
/*
  Henke Tables

  The original data can be found at:
    ftp://grace.lbl.gov/pub/sf/

  From ftp://grace.lbl.gov/pub/sf/read.me:

                Low-Energy X-ray Interaction Coefficients:
                Photoabsorption, Scattering, and Reflection
                        E = 30-30,000 eV, Z = 1-92

                B. L. Henke, E. M. Gullikson, and J. C. Davis
                        Center for X-Ray Optics, 2-400
                        Lawrence Berkeley Laboratory
                        Berkeley, California 94720

  Reference: B. L. Henke, E. M. Gullikson, and J. C. Davis,
  Atomic Data and Nuclear Data Tables Vol. 54 No. 2 (July 1993).
 */
"""

def collect_tables():
  nff_files = []
  for file in os.listdir(reference_tables_directory):
    fn = file.lower().capitalize()
    if (fn[-4:] == ".nff"): nff_files.append(file)
  tables = [0] * 120
  for file in nff_files:
    f = join_open(reference_tables_directory, file, "r")
    header = f.readline()
    table = f.readlines()
    f.close()
    Symbol = header[1:3].strip()
    Z = int(header[7:9])
    assert len(Symbol) > 0
    assert Symbol[0] in string.lowercase
    assert Symbol[-1] in string.lowercase
    assert Z > 0 and Z < len(tables)
    assert tables[Z] == 0
    Symbol = Symbol.capitalize()
    tables[Z] = (Symbol, table)
  Z = tables[1:].index(0) + 1
  rest = tables[Z:]
  assert rest == [0] * len(rest)
  tables = tables[:Z]
  return tables

def print_table_block(f, tables, Z_begin, Z_end, define_noval=0):
  print >> f, "namespace table_data {"
  print >> f
  print >> f, "using anomalous::e_fp_fdp;"
  print >> f
  # Visual C++ 7.0 compilation is very slow with define_noval=1
  if (define_noval): print >> f, "#define NOVAL fp_fdp_undefined"
  for Z in xrange(Z_begin, Z_end):
    tab = tables[Z]
    print >> f, "e_fp_fdp " + tab[0].lower() \
      + "[] = { /* Z = " + str(Z) + " */"
    for line in tab[1]:
      flds = line.split()
      assert len(flds) == 3
      if (define_noval and flds[1] == "-9999.00"): flds[1] = "NOVAL"
      print >> f, "{%s, %s, %s}," % tuple(flds)
    print >> f, "{0, 0, 0}"
    print >> f, "};"
    print >> f
  if (define_noval):
    print >> f, "#undef NOVAL"
    print >> f
  print >> f, "} // namespace table_data"
  print >> f
  print >> f, "}}} // namespace cctbx::eltbx::henke"

def print_henke_cpp(f, tables):
  print >> f, "namespace table_data {"
  print >> f
  print >> f, "using anomalous::e_fp_fdp;"
  print >> f
  for tab in tables[1:]:
    print >> f, "extern e_fp_fdp " + tab[0].lower() + "[];"
  print >> f
  print >> f, "static const anomalous::label_z_e_fp_fdp all[] = {"
  i = 0
  for tab in tables[1:]:
    i += 1
    print >> f, \
      "{\"" + tab[0] + "\", " + str(i) + ", " + tab[0].lower() + "},"
  print >> f, "{0, 0, 0}"
  print >> f, "};"
  print >> f, """
  } // namespace table_data

  table::table(
    std::string const& label,
    bool exact,
    bool exception_if_no_match)
  {
    std::string work_label = basic::strip_label(label, exact);
    label_z_e_fp_fdp_ = anomalous::find_entry(
      table_data::all, work_label, exact, exception_if_no_match);
  }

  fp_fdp
  table::at_ev(double energy) const
  {
    fp_fdp raw = anomalous::interpolate(label_z_e_fp_fdp_, energy);
    if (!raw.is_valid_fp()) return raw;
    // subtract the number of electrons
    return fp_fdp(raw.fp() - label_z_e_fp_fdp_->z, raw.fdp());
  }

  table_iterator::table_iterator()
  :
    current_("H", true)
  {}

  table
  table_iterator::next()
  {
    table result = current_;
    if (current_.is_valid()) current_.label_z_e_fp_fdp_++;
    return result;
  }

}}} // namespace cctbx::eltbx::henke"""

def collect_points(lines):
  points = []
  for line in lines:
    points.append(line.split()[0])
  return points

def collect_tab_points(tables):
  tab_points = []
  for tab in tables[1:]:
    tab_points.append(collect_points(tab[1]))
  return tab_points

def compare_points(tables):
  tab_points = collect_tab_points(tables)
  for i in xrange(len(tab_points)-1):
    for j in xrange(i+1, len(tab_points)):
      if (tab_points[i] == tab_points[j]):
        print "points %d==%d" % (i+1,j+1)

def run(target_dir):
  tables = collect_tables()
  compare_points(tables) # establish that each list of data points is unique
  f = join_open(target_dir, "henke.cpp", "w")
  print_header(f)
  print_ftp_info(f)
  print_henke_cpp(f, tables)
  f.close()
  Z_block = 12
  for Z_begin in xrange(1, len(tables), Z_block):
    Z_end = min(len(tables), Z_begin + Z_block)
    f = join_open(
      target_dir, "henke_tables_%02d_%02d.cpp" % (Z_begin, Z_end-1), "w")
    print_header(f)
    print_table_block(f, tables, Z_begin, Z_end)
    f.close()

if (__name__ == "__main__"):
  run(".")
