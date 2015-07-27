from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.biomt_reconstruction

from iotbx import pdb
from libtbx.utils import Sorry
from  iotbx.pdb.multimer_reconstruction import multimer
import sys, os

def biomt_reconstruction(args, command_name="phenix.pdb.biomt_reconstruction"):
  ''' (string) -> file

  Reconstruction of the biological assembly multimer by applying
  BIOMT transformation, from the REMARK 350 in the pdb file, to all
  chains in file.

  Produce a new pdb file containg the biological assembly multimer information

  Arguments:
  args: a string containing the pdb file name, 'file_name.pdb' and an
        optional output file name. If no file name is given, the output
        file will have the name 'biological_assembly_file_name.pdb'

  >>> phenix.pdb.biomt_reconstruction input_file_name.pdb
  bio-assembly reconstructed pdb file was added to your current directory
  /path/biological_assembly_input_file_name.pdb

   >>> phenix.pdb.biomt_reconstruction input_file_name.pdb output_file_name.pdb
  bio-assembly reconstructed pdb file was added to your current directory
  /path/output_file_name.pdb

  @author: Youval Dar (LBL 2013)
  '''
  if (len(args) == 0):
    raise Sorry('No input filename is given. Please provide a pdb file name. \n' + \
                '>>> phenix.pdb.biomt_reconstruction input_file_name.pdb')
  elif ('--help' in args) or ('-h' in args):
    print help(biomt_reconstruction)
  elif (len(args) > 2):
    raise Sorry('To many input parameters are given')
  elif not os.path.isfile(args[0]):
    print 'There is no file {} in the current directory \n {}'.format(args[0],os.getcwd())
    raise Sorry('File name error')
  else:
    m = multimer(file_name=args[0],reconstruction_type='ba')
    if m.number_of_transforms > 0:
      # write the bio-assembly reconstructed file in current directory
      if len(args) == 2:
        # output file name given
        m.write(args[1])
        print '$s file was added to your current directory'%args[1]
      else:
        # output file name not given
        m.write()
        print 'bio-assembly reconstructed pdb file was added to your current directory'
        print os.getcwd() + '/' + m.pdb_output_file_name
    else:
      print 'No BIOMT information in {}'.format(args[0])

if (__name__ == "__main__"):
  biomt_reconstruction(args=sys.argv[1:])
