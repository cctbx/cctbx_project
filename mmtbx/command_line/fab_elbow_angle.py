from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.fab_elbow_angle
from mmtbx.utils.fab_elbow_angle import fab_elbow_angle
from libtbx.utils import Sorry
import iotbx.pdb
import sys
import os

def run(args, log=sys.stdout):
  """
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
  >>>phenix.fab_elbow_angle fab_name.pdb [light=L] [heavy=H] [limit_l=107] [limit_h=113]

  examples:
  >>>phenix.fab_elbow_angle 7fab.pdb light=L heavy=H limit_l=104 limit_h=117
  >>>123.00
  >>>phenix.fab_elbow_angle 1bbd.pdb
  >>>126.34
  """
  if(len(args)==0 or (len(args)==1 and
                          args[0].lower() in ("-h", "--h","-help","--help"))):
    print >> log, "-"*79
    print >> log, run.__doc__
    print >> log, "-"*79
    sys.exit(0)
  elif not os.path.isfile(args[0]):
    print >> log, 'There is no file %s in the current directory'%args[0]
    raise Sorry('File name error')
  else:
    ph = iotbx.pdb.input(file_name=args[0]).construct_hierarchy()
    inp_params = dict()
    # fab_elbow_angle.py input parameters has different names
    param_manes = dict(light='chain_id_light',
               heavy='chain_id_heavy',
               limit_l='limit_light',
               limit_h='limit_heavy')
    for x in args[1:]:
      kv = x.split('=')
      if kv[0].lower() in ('light','heavy','limit_l','limit_h') and len(kv)==2:
        try:
          val = int(kv[1])
        except ValueError:
          val = kv[1]
        inp_params[param_manes[kv[0]]] = val
      else:
        raise Sorry('Option {} does not exist'.format(kv[0]))

    elbow_angle = fab_elbow_angle(
      pdb_hierarchy=ph,
      **inp_params).fab_elbow_angle
    print  >> log, 'fab elbow angle (deg): {0:.2f}'.format(elbow_angle)

if __name__ == "__main__":
  run(args=sys.argv[1:])
