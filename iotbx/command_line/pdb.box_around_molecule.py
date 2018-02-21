from __future__ import division

import sys
from cctbx import uctbx
import iotbx.pdb
from libtbx.option_parser import option_parser
import libtbx.load_env
from libtbx.str_utils import show_string
from libtbx.utils import date_and_time
from libtbx.utils import null_out

def run(args):
  log = sys.stdout
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage="%s [options] pdb_file" % libtbx.env.dispatcher_name)
    .option(None, "--buffer_layer",
      action="store",
      type="float",
      default=5)
  ).process(args=args, nargs=1)
  pdb_inp = iotbx.pdb.input(file_name=command_line.args[0])
  atoms = pdb_inp.atoms()
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=atoms.extract_xyz(),
    buffer_layer=command_line.options.buffer_layer)
  atoms.set_xyz(new_xyz=box.sites_cart)
  print >> log, 'REMARK %s --buffer-layer=%.6g %s' % (
    libtbx.env.dispatcher_name,
    command_line.options.buffer_layer,
    show_string(command_line.args[0]))
  print >> log, 'REMARK %s' % date_and_time()
  iotbx.pdb.write_whole_pdb_file(
      output_file=log,
      pdb_hierarchy=pdb_inp.construct_hierarchy(),
      crystal_symmetry=box.crystal_symmetry(),
      ss_annotation=pdb_inp.extract_secondary_structure(log=null_out()))

if (__name__ == "__main__"):
  run(sys.argv[1:])
