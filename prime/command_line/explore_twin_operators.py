from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME prime.explore_twin_operators

from prime.index_ambiguity.mod_indexing_ambiguity import indamb_handler
import argparse, os
from prime.postrefine.mod_input import read_frame, read_pickles
import six

def main(data, only_merohedral):
  indambh = indamb_handler()
  intFileList = read_pickles([data])
  if intFileList:
    obsList = {}
    for intFileName in intFileList:
      intPickle = read_frame(intFileName)
      try:
        obs = intPickle['observations'][0]
        obsList[intFileName] = obs
      except Exception as e:
        print("Warning:", e)
        pass
    for key,value in six.iteritems(obsList):
      if only_merohedral:
        flag_all = False
      else:
        flag_all = True
      ops = indambh.generate_twin_operators(value, flag_all=flag_all)
      if ops:
        print(os.path.basename(key), '%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f'%value.unit_cell().parameters(), ' '.join([op.operator.r().as_hkl() for op in ops]))
      else:
        print(os.path.basename(key), '%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f'%value.unit_cell().parameters(), 'Twining operators not found')


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description='Explore (pseudo)merohedral twinning operators'
  )
  parser.add_argument(
    'integration_pickles',
    metavar='Integration Pickles',
    help='Path to integration pickles (regex. is acceptable).'
  )
  parser.add_argument('--only_merohedral', dest='only_merohedral', action='store_const',
                    const=True, default=False,
                    help='Flag for showing only merohedral twinnings')
  args = parser.parse_args()
  if args.only_merohedral:
    print("Only showing results with merohedral twinning operators")
  else:
    print("Showing all possible pseudo- and true-merohedral twinning operators")
  args = parser.parse_args()

  main(args.integration_pickles, args.only_merohedral)
