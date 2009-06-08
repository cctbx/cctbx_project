# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_interpretation

from mmtbx.monomer_library import pdb_interpretation
from libtbx.option_parser import libtbx_option_parser
from libtbx.utils import date_and_time
from libtbx.str_utils import show_string
import sys, os
op = os.path

def run(args):
  command_line = (libtbx_option_parser(
    usage="phenix.pdb_interpretation [options] pdb_file_or_directory")
    .option(None, "--write_geo_files",
      action="store_true",
      default=False)
    .option(None, "--write_tardy_geo_files",
      action="store_true",
      default=False)
    .option(None, "--strict_conflict_handling",
      action="store_true",
      default=False)
  ).process(args=args)
  co = command_line.options
  processed_pdb_files = pdb_interpretation.run(
    args=command_line.args,
    strict_conflict_handling=co.strict_conflict_handling,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    return_all_processed_pdb_files=True)
  print
  if (   co.write_geo_files
      or co.write_tardy_geo_files):
    for processed_pdb_file in processed_pdb_files:
      acp = processed_pdb_file.all_chain_proxies
      source = acp.pdb_inp.source_info()
      assert source.startswith("file ")
      pdb_file_name = source[5:]
      sites_cart = acp.sites_cart_exact()
      site_labels = [atom.id_str() for atom in acp.pdb_atoms]
      def write_geo(label, geo, geo_file_name):
        f = open(geo_file_name, "w")
        print >> f, "# %sgeometry restraints for file:" % label
        print >> f, "#   %s" % show_string(pdb_file_name)
        print >> f, "#", date_and_time()
        print >> f
        geo.show_sorted(
          sites_cart=sites_cart, site_labels=site_labels, f=f)
        f.close()
      geo = processed_pdb_file.geometry_restraints_manager()
      if (co.write_geo_files):
        geo_file_name = op.basename(pdb_file_name) + ".geo"
        print "Writing file: %s" % show_string(geo_file_name)
        write_geo(label="", geo=geo, geo_file_name=geo_file_name)
        print
      if (co.write_tardy_geo_files):
        geo_file_name = op.basename(pdb_file_name) + ".tardy.geo"
        print "Writing file: %s" % show_string(geo_file_name)
        tardy_tree = geo.construct_tardy_tree()
        reduced_geo = geo.reduce_for_tardy(tardy_tree=tardy_tree)
        write_geo(label="tardy ", geo=reduced_geo, geo_file_name=geo_file_name)
        print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
