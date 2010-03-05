def run(args):
  if (len(args) == 0): args = ["--help"]
  from libtbx.option_parser import option_parser
  import libtbx.load_env
  command_line = (option_parser(
    usage="%s [options] pdb_file" % libtbx.env.dispatcher_name)
    .option(None, "--buffer_layer",
      action="store",
      type="float",
      default=5)
  ).process(args=args, nargs=1)
  import iotbx.pdb
  pdb_inp = iotbx.pdb.input(file_name=command_line.args[0])
  atoms = pdb_inp.atoms()
  from cctbx import uctbx
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=atoms.extract_xyz(),
    buffer_layer=command_line.options.buffer_layer)
  atoms.set_xyz(new_xyz=box.sites_cart)
  from libtbx.str_utils import show_string
  print 'REMARK %s --buffer-layer=%.6g %s' % (
    libtbx.env.dispatcher_name,
    command_line.options.buffer_layer,
    show_string(command_line.args[0]))
  from libtbx.utils import date_and_time
  print 'REMARK %s' % date_and_time()
  print iotbx.pdb.format_cryst1_record(crystal_symmetry=box.crystal_symmetry())
  print pdb_inp.construct_hierarchy().as_pdb_string(append_end=True),

if (__name__ == "__main__"):
  import sys
  run(sys.argv[1:])
