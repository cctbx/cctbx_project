from __future__ import division, print_function

import sys
from cctbx import uctbx
import iotbx.pdb
from libtbx.option_parser import option_parser
import libtbx.load_env
from libtbx.str_utils import show_string
from libtbx.utils import date_and_time
import mmtbx.model

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
  model = mmtbx.model.manager(
      model_input = pdb_inp)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=model.get_sites_cart(),
    buffer_layer=command_line.options.buffer_layer)
  model.set_sites_cart(box.sites_cart)
  # Bad hack, never repeat. In fact, all the boxing functionality should
  # go into mmtbx.model.manager
  model._crystal_symmetry = box.crystal_symmetry()
  print('REMARK %s --buffer-layer=%.6g %s' % (
    libtbx.env.dispatcher_name,
    command_line.options.buffer_layer,
    show_string(command_line.args[0])), file=log)
  print('REMARK %s' % date_and_time(), file=log)
  print(model.model_as_pdb(), file=log)

if (__name__ == "__main__"):
  run(sys.argv[1:])
