from __future__ import absolute_import, division, print_function
from scitbx.source_generators.utils import join_open
from scitbx.source_generators.utils import write_this_is_auto_generated
import libtbx.load_env
import string
from six.moves import range

this = "cctbx.source_generators.eltbx.generate_sasaki_cpp"

reference_tables_directory = libtbx.env.under_dist(
  "cctbx", "reference/sasaki")

def print_header(f):
  write_this_is_auto_generated(f, this)
  print("""\
#include <cctbx/eltbx/sasaki.h>
#include <scitbx/constants.h>

namespace cctbx { namespace eltbx { namespace sasaki {

namespace table_data {

using namespace detail;
""", file=f)

def print_ftp_info(f):
  print("""\
/*
  Sasaki Tables

  Scattering factors based on the Cromer and Liberman method.
  Original data can be downloaded from:
  ftp://pfweis.kek.jp/pub/Sasaki-table/
  Any reports or publications of these data will acknowledge
  its use by the citation:
    Anomalous scattering factors
        S.Sasaki (1989) Numerical Tables of Anomalous Scattering Factors
        Calculated by the Cromer and Liberman Method,
        KEK Report, 88-14, 1-136
  Questions about these data should be addressed to Dr.Satoshi Sasaki,
  Tokyo Institute of Technology.  Email: sasaki@nc.titech.ac.jp
 */
""", file=f)

class sasaki_table(object):

  def __init__(self, atomic_symbol, atomic_number,
               edge_label=0, edge_wave_length=0):
    assert edge_label in (0, "K", "L1", "L2", "L3")
    self.atomic_symbol = atomic_symbol.lower().capitalize()
    self.atomic_number = atomic_number
    self.edge_label = edge_label
    self.edge_wave_length = edge_wave_length
    self.fp = []
    self.fdp = []

  def check_first_last(self, i_block, first, last):
    assert abs(
      self.first + (i_block * self.incr1000) / 1000. - first) < 1.e-8
    assert abs(first + 9 * self.incr1000 / 10000. - last) < 1.e-8

  def validate(self):
    assert len(self.fp) == len(self.fdp)
    assert len(self.fp) == 280

def collect_tables(file_object, edge):
  import re
  tables = []
  table = 0
  for line in file_object:
    if (not edge):
      m = re.match(
        r"ATOMIC SYMBOL\s+=\s+(\S+)\s+ATOMIC NUMBER\s+=\s+(\d+)", line)
    else:
      m = re.match(
        r"ATOMIC SYMBOL\s+=\s+(\S+)\s+ATOMIC NUMBER\s+=\s+(\d+)"
        + r"\s+(\S+)\s+ABSORPTION EDGE\s+\(\s*(\S+)\s+A", line)
    if (m):
      if (table != 0):
        table.validate()
        tables.append(table)
      if (not edge):
        table = sasaki_table(m.group(1), m.group(2))
      else:
        table = sasaki_table(m.group(1), m.group(2), m.group(3), m.group(4))
      i_block = 0
      continue
    flds = line.split()
    if (flds[1] == "TO"):
      assert flds[3] == "F'"
      first, last = [float(flds[i]) for i in (0, 2)]
      flds = flds[4:]
      if (i_block == 0):
        table.first = 0.1
        table.incr1000 = 100
      table.check_first_last(i_block, first, last)
      i_block += 1
      table.fp += flds
    elif (flds[0].find(",") >= 0):
      assert flds[1] == "F'"
      first, last = [float(x) for x in flds[0].split(",")]
      flds = flds[2:]
      if (i_block == 0):
        table.first = first
        table.incr1000 = 1
      table.check_first_last(i_block, first, last)
      i_block += 1
      table.fp += flds
    else:
      assert flds[0] == 'F"'
      flds = flds[1:]
      table.fdp += flds
    assert len(flds) == 10
  if (table != 0):
    table.validate()
    tables.append(table)
  return tables

class table_references(object):

  def __init__(self, wide):
    self.wide = wide
    self.k = 0
    self.l1 = 0
    self.l2 = 0
    self.l3 = 0

def combine_tables(tables_wide, tables_k, tables_l):
  ctab_dict = {}
  for w in tables_wide:
    ctab_dict[w.atomic_symbol] = table_references(w)
  for k in tables_k:
    ctab_dict[k.atomic_symbol].k = k
  for l in tables_l:
    if (l.edge_label == "L1"):
      ctab_dict[l.atomic_symbol].l1 = l
    elif (l.edge_label == "L2"):
      ctab_dict[l.atomic_symbol].l2 = l
    else:
      assert l.edge_label == "L3"
      ctab_dict[l.atomic_symbol].l3 = l
  ctab_list = []
  for w in tables_wide:
    ctab_list.append(ctab_dict[w.atomic_symbol])
  return ctab_list

def print_table_block(f, tables_combined, i_begin, i_end):
  for i_table in range(i_begin, i_end):
    ctab = tables_combined[i_table]
    for edge in (ctab.wide, ctab.k, ctab.l1, ctab.l2, ctab.l3):
      if (not edge): continue
      lbl = edge.atomic_symbol.lower()
      ann = "Z = " + str(ctab.wide.atomic_number)
      if (edge.edge_label):
        lbl += "_" + edge.edge_label.lower()
        ann += "; edge at " + edge.edge_wave_length + " A"
      print("raw " + lbl + "[] = { // " + ann, file=f)
      for i in range(len(edge.fp)):
        print("{%s,%s}," % (edge.fp[i], edge.fdp[i]), file=f)
      print("};", file=f)
  print(file=f)
  print("} // namespace table_data", file=f)
  print(file=f)
  print("}}} // namespace cctbx::eltbx::sasaki", file=f)

def print_sasaki_cpp(f, tables_combined):
  for ctab in tables_combined:
    print("extern raw " + ctab.wide.atomic_symbol.lower() + "[];", file=f)
    for edge in (ctab.k, ctab.l1, ctab.l2, ctab.l3):
      if (edge):
        assert edge.atomic_symbol == ctab.wide.atomic_symbol
        print("extern raw " + edge.atomic_symbol.lower() \
          + "_" + edge.edge_label.lower() + "[];".lower(), file=f)
  print(file=f)
  print("static const detail::info all[] = {", file=f)
  i = 0
  for ctab in tables_combined:
    i += 1
    out = '{"' + ctab.wide.atomic_symbol \
        + '", ' + str(ctab.wide.atomic_number)
    out += ", " + ctab.wide.atomic_symbol.lower()
    for edge in (ctab.k, ctab.l1, ctab.l2, ctab.l3):
      if (edge):
        out += ", %.4f" % (edge.first,)
        out += ", " + edge.atomic_symbol.lower()
        out += "_" + edge.edge_label.lower()
      else:
        out += ", 0., 0"
    out += "},"
    print(out, file=f)
  print("{0, 0, 0, 0., 0, 0., 0, 0., 0, 0., 0}", file=f)
  print("};", file=f)
  print("""
  } // namespace table_data

  table::table(
    std::string const& label,
    bool exact,
    bool exception_if_no_match)
  {
    std::string work_label = basic::strip_label(label, exact);
    info_ = anomalous::find_entry(
      table_data::all, work_label, exact, exception_if_no_match);
  }

  namespace detail {
  namespace {

    long
    find_table_interval(double given, double first, double incr,
                        double tolerance = 1.e-8)
    {
      double span = (n_raw - 1) * incr;
      double f = (given - first) / span;
      if (f < -tolerance || f > (1.+tolerance)) return -1;
      long i = static_cast<long>(f * (n_raw - 1));
      if (i == n_raw - 1) i--;
      return i;
    }

    bool
    interpolate(double given, double first, const raw* table, bool edge,
                fp_fdp& result)
    {
      if (!table) return false;
      double incr;
      if (!edge) incr = wide_incr;
      else       incr = edge_incr;
      long i = find_table_interval(given, first, incr);
      if (i < 0) return false;
      double x = (given - first) / incr - double(i);
      double fp  = table[i].fp  + x * (table[i+1].fp  - table[i].fp);
      double fdp = table[i].fdp + x * (table[i+1].fdp - table[i].fdp);
      result = fp_fdp(fp, fdp);
      return true;
    }

  } // namespace <anonymous>
  } // namespace detail

  fp_fdp
  table::at_ev(double energy) const
  {
    using detail::interpolate;
    fp_fdp result(fp_fdp_undefined, fp_fdp_undefined);
    double given = scitbx::constants::factor_ev_angstrom / energy;
    if (interpolate(given, info_->first_k, info_->k, true, result)) {
      return result;
    }
    if (interpolate(given, info_->first_l1, info_->l1, true, result)) {
      return result;
    }
    if (interpolate(given, info_->first_l2, info_->l2, true, result)) {
      return result;
    }
    if (interpolate(given, info_->first_l3, info_->l3, true, result)) {
      return result;
    }
    interpolate(given, detail::first_wide, info_->wide, false, result);
    return result;
  }

  table_iterator::table_iterator()
  :
    current_("Li", true)
  {}

  table
  table_iterator::next()
  {
    table result = current_;
    if (current_.is_valid()) current_.info_++;
    return result;
  }

}}} // namespace cctbx::eltbx::sasaki""", file=f)

def run(target_dir):
  f = join_open(reference_tables_directory, "fpwide.tbl", "r")
  tables_wide = collect_tables(f, 0)
  f.close()
  f = join_open(reference_tables_directory, "fpk.tbl", "r")
  tables_k = collect_tables(f, 1)
  f.close()
  f = join_open(reference_tables_directory, "fpl.tbl", "r")
  tables_l = collect_tables(f, 1)
  f.close()
  tables_combined = combine_tables(tables_wide, tables_k, tables_l)
  f = join_open(target_dir, "sasaki.cpp", "w")
  print_header(f)
  print_ftp_info(f)
  print_sasaki_cpp(f, tables_combined)
  f.close()
  block_size = 12
  for i_begin in range(0, len(tables_combined), block_size):
    i_end = min(len(tables_combined), i_begin + block_size)
    f = join_open(
      target_dir, "sasaki_tables_%02d_%02d.cpp" % (i_begin+1, i_end), "w")
    print_header(f)
    print_table_block(f, tables_combined, i_begin, i_end)
    f.close()

if __name__ == "__main__":
  run(".")
