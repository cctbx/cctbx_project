
from libtbx.utils import Sorry, Usage
from libtbx import str_utils
import cStringIO
import os
import sys

class summary (object) :
  def __init__ (self, pdb_hierarchy=None, pdb_file=None, sites_cart=None,
      keep_hydrogens=False) :
    if (pdb_hierarchy is None) :
      assert (pdb_file is not None)
      from iotbx import file_reader
      pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
      pdb_in.assert_file_type("pdb")
      pdb_hierarchy = pdb_in.file_object.construct_hierarchy()
      pdb_hierarchy.atoms().reset_i_seq()
    if (sites_cart is not None) :
      pdb_hierarchy.atoms().set_xyz(sites_cart)
    from mmtbx.validation import ramalyze, rotalyze, cbetadev, clashscore
    log = cStringIO.StringIO()
    rama = ramalyze.ramalyze()
    rama.analyze_pdb(hierarchy=pdb_hierarchy)
    if (rama.numtotal > 0) :
      rama_out_count, rama_out_percent = rama.get_outliers_count_and_fraction()
      rama_fav_count, rama_fav_percent = rama.get_favored_count_and_fraction()
      self.rama_fav = rama_fav_percent * 100.0
      self.rama_out = rama_out_percent * 100.0
    else :
      self.rama_fav = None
      self.rama_out = None
    rota = rotalyze.rotalyze()
    rota.analyze_pdb(hierarchy=pdb_hierarchy)
    if (rota.numtotal > 0) :
      rota_count, rota_perc = rota.get_outliers_count_and_fraction()
      self.rota_out = rota_perc * 100.0
    else :
      self.rota_out = None
    cs = clashscore.clashscore()
    clash_dict, clash_list = cs.analyze_clashes(hierarchy=pdb_hierarchy,
      keep_hydrogens=keep_hydrogens)
    self.clashscore = clash_dict['']
    cbeta = cbetadev.cbetadev()
    cbeta_txt, cbeta_summ, cbeta_list = cbeta.analyze_pdb(
      hierarchy=pdb_hierarchy,
      outliers_only=True)
    # TODO leave self.cbeta_out as None if not protein
    self.cbeta_out = len(cbeta_list)
    self.r_work = None
    self.r_free = None
    self.rms_bonds = None
    self.rms_angles = None
    self.d_min = None
    if (pdb_file is not None) :
      from iotbx.pdb import extract_rfactors_resolutions_sigma
      published_results = extract_rfactors_resolutions_sigma.extract(
        file_name=pdb_file)
      self.r_work = published_results.r_work
      self.r_free = published_results.r_free
      self.d_min   = published_results.high
      lines = open(pdb_file).readlines()
      for line in lines :
        if (line.startswith("REMARK Final:")) :
          fields = line.strip().split()
          self.rms_bonds = float(fields[-4])
          self.rms_angles = float(fields[-1])
          break

  def show (self, out=sys.stdout, prefix="  ") :
    def fs (format, value) :
      return str_utils.format_value(format, value, replace_none_with=("(none)"))
    print >> out, "%sRamachandran outliers = %s %%" % (prefix,
      fs("%6.2f", self.rama_out))
    print >> out, "%s             favored  = %s %%" % (prefix,
      fs("%6.2f", self.rama_fav))
    print >> out, "%sRotamer outliers      = %s %%" % (prefix,
      fs("%6.2f", self.rota_out))
    print >> out, "%sC-beta deviations     = %s" % (prefix,
      fs("%6d", self.cbeta_out))
    print >> out, "%sClashscore            = %6.2f" % (prefix,self.clashscore)
    if (self.r_work is not None) :
      print >> out, "%sR-work                = %8.4f" % (prefix, self.r_work)
    if (self.r_free is not None) :
      print >> out, "%sR-free                = %8.4f" % (prefix, self.r_free)
    if (self.rms_bonds is not None) :
      print >> out, "%sRMS(bonds)            = %8.4f" % (prefix, self.rms_bonds)
    if (self.rms_angles is not None) :
      print >> out, "%sRMS(angles)           = %6.2f" % (prefix,
        self.rms_angles)
    if (self.d_min is not None) :
      print >> out, "%sHigh resolution       = %6.2f" % (prefix, self.d_min)

def run (args, out=sys.stdout) :
  if (len(args) == 0) :
    raise Usage("""
mmtbx.validation_summary model.pdb

Prints a brief summary of validation criteria, including Ramachandran
statistics, rotamer outliers, clashscore, C-beta deviations, plus R-factors
and RMS(bonds)/RMS(angles) if found in PDB header.  (This is primarily used
for evaluating the output of refinement tests; general users are advised to
run phenix.model_vs_data or the validation GUI.)
""")
  pdb_file = args[0]
  if (not os.path.isfile(pdb_file)) :
    raise Sorry("Not a file: %s" % pdb_file)
  s = summary(pdb_file=pdb_file)
  print >> out, ""
  print >> out, "Validation summary for %s:" % pdb_file
  s.show(out=out)
  print >> out, ""

if (__name__ == "__main__") :
  run(sys.argv[1:])
