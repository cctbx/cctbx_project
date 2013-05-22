from __future__ import division

def write_disable_nomenclature_errors (f) :
  f.write("try :\n")
  f.write("  set_nomenclature_errors_on_read(\"ignore\")\n")
  f.write("except Exception :\n")
  f.write("  pass\n")
