"""Calculate FAB elbow angle"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.fab_elbow_angle
from mmtbx.utils.fab_elbow_angle import fab_elbow_angle
import iotbx.pdb
import sys
import iotbx.phil

def get_master_phil():
  import libtbx.phil
  return libtbx.phil.parse(
    input_string="""
      model_file_name = None
        .type = path
        .multiple = True
        .help = '''Enter a PDB file name'''
      light = 'L'
        .type = str
        .help = '''chain ID of light domain'''
      heavy = 'H'
        .type = str
        .help = '''chain ID of heavy domain'''
      limit_l = 107
        .type = int
        .help = the number of the residue separating between variable \
            and constant domains in the light chain
      limit_h = 113
        .type = int
        .help = the number of the residue separating between variable \
            and constant domains in the heavy chain, by default 113
""")

def usage():
  return """
  Calculating Fragment Antigen-Binding (Fab) elbow angle

  Each FAB is made of two chains, Heavy (H) and Light (L), and two domains,
  Variable (V) and Constant (C or C1, The heavy domain has three parts).
  The angle between the variable and constant domains is the elbow angle.
  It is the angle between rotation axes, the pseudo-dyad axes,
  defined by aligning the light portion of each domain, on-to the heavy one.

  command line options:
  light        chain ID of light domain, by default set to L
  heavy        chain ID of heavy domain, by default set to H
  limit_l      the number of the residue separating between variable
               and constant domains in the light chain, by default 107
  limit_H      the number of the residue separating between variable
               and constant domains in the heavy chain, by default 113
  h            help

  usage:
  phenix.fab_elbow_angle fab_name.pdb [light=L] [heavy=H] [limit_l=107] [limit_h=113]

  examples:
  >>>phenix.fab_elbow_angle 7fab.pdb light=L heavy=H limit_l=104 limit_h=117
  >>>123.00
  >>>phenix.fab_elbow_angle 1bbd.pdb
  >>>126.34
  """

def run(args, log=sys.stdout):
  master_phil = get_master_phil()
  input_objects = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="model_file_name")
  work_params = input_objects.work.extract()
  if len(work_params.model_file_name) != 1:
    print(usage(), file=log)
    return

  ph = iotbx.pdb.input(file_name=work_params.model_file_name[0]).construct_hierarchy()
  elbow_angle = fab_elbow_angle(
      pdb_hierarchy=ph,
      chain_id_light=work_params.light,
      chain_id_heavy=work_params.heavy,
      limit_light=work_params.limit_l,
      limit_heavy=work_params.limit_h).fab_elbow_angle
  print('fab elbow angle (deg): {0:.2f}'.format(elbow_angle), file=log)

if __name__ == "__main__":
  run(args=sys.argv[1:])

