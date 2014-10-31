from __future__ import division
import iotbx.ncs
import time
import sys
import os
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.ncs_from_pdb


def ncs_from_pdb(args):
  """
  Search for NCS relations in a pdb file

  Arguments and Options:
    file_name (str): xxxx.pdb (must be the first parameter)
    write_messages (True): Print messages produced while searching for NCS
    max_rmsd (10.0): Max allowed rmsd between matching chains
    check_atom_order (False): Ensure that atoms of matching residues are in
      the same order
    exclude_misaligned_residues (True): Exclude matching residues if they are
      spatially misaligned
    max_dist_diff (5.0): Max allowed spatially misaligned distance difference
    min_percent (0.80): Min percent of
      (number of matching residues / number of residues in longer chain)
    min_contig_length (10): Minimum length of matching chain segments

  Example:
  >>> phenix.pdb.ncs_from_pdb 4c5q.pdb min_percent=0.9

  """
  t0 = time.time()
  if (args == []) or ('-h' in args) or ('--help' in args):
    print ncs_from_pdb.__doc__
  else:
    # set parameters
    # clean and remove '=' signs
    args = [y for x in args for y in x.split('=') if y != '']
    params = ncs_search_params()
    # file name must be the first parameter
    params.file_name = args.pop(0)
    if not os.path.isfile(params.file_name):
      print 'No PDB file found.'
      print 'Make sure the PDB file name is the first parameter:'
      print '>>> phenix.ncs_from_pdb xxxx.pdb\n'
    else:
      bool_param = [
        'write_messages','check_atom_order','exclude_misaligned_residues']
      float_param =['max_rmsd','max_dist_diff','min_percent']
      int_param = ['min_contig_length']
      input_is_ok = True
      while args:
        key = args.pop(0)
        if hasattr(params,key):
          val = args.pop(0)
          if key in bool_param:
            val = (val.lower() == 'true')
          elif key in float_param:
            val = float(val)
          elif key in int_param:
            val = int(val)
          params.__dict__[key] = val
        else:
          print '{} is not a valid option!!!'.format(key)
          input_is_ok = False
          break
      if input_is_ok:
        # correct parameters
        if params.min_percent > 1:
          params.min_percent /= 100
        # call function
        ncs_obj = iotbx.ncs.input(
            file_name=params.file_name,
            use_cctbx_find_ncs_tools=True,
            use_minimal_master_ncs=True,
            min_contig_length=params.min_contig_length,
            write_messages=params.write_messages,
            max_rmsd=params.max_rmsd,
            check_atom_order=params.check_atom_order,
            exclude_misaligned_residues=params.exclude_misaligned_residues,
            max_dist_diff=params.max_dist_diff,
            min_percent=params.min_percent)
        print 'File name: {}'.format(params.file_name)
        print 'Search parameters values'
        print '-'*40
        for k,v in params.__dict__.iteritems():
          print '{0:<30}: {1:<20}'.format(k,str(v))

        if len(ncs_obj.ncs_transform) == 0:
          print 'No NCS relations found'
        ncs_obj.get_ncs_info_as_spec(write=True)
        print '\nTime to process: {0:.2f} sec'.format(time.time()-t0)

class ncs_search_params(object):

  def __init__(self):
    """ set default parameters for ncs_from_pdb """
    self.file_name = ''
    self.write_messages = False
    self.max_rmsd = 10.0
    self.check_atom_order = False
    self.exclude_misaligned_residues = True
    self.max_dist_diff = 5.0
    self.min_percent = 0.80
    self.min_contig_length = 10


if __name__ == '__main__':
  ncs_from_pdb(args=sys.argv[1:])
