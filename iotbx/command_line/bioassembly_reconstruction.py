from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.bioassembly_reconstruction

from iotbx import pdb
from libtbx.utils import Sorry
from  iotbx.pdb.multimer_reconstruction import multimer
import sys, os

def run(args, command_name="phenix.pdb.bioassembly_reconstruction"):
  ''' (string) -> file
  
  Reconstruction of the biological assembly multimer by applying
  BIOMT transformation, from the pdb file, to all chains in file.
  
  Produce a new pdb file containg the biological assembly multimer information
  
  Arguments:
  args: a string containing the pdb file name, 'file_name.pdb' and an optional output file name
        if no file name is given, the output file will have the name 'copy_file_name.pdb'
        
  >>> phenix.pdb.bioassembly_reconstruction 'input_file_name.pdb'
  bio-assembly reconstructed pdb file was added to your current directory
  /path/copy_input_file_name.pdb
  
   >>> phenix.pdb.bioassembly_reconstruction 'input_file_name.pdb' 'output_file_name.pdb'
  bio-assembly reconstructed pdb file was added to your current directory
  /path/output_file_name.pdb
  
  '''
  if (len(args) == 0): 
    raise Sorry('No input filename is given. Please provide a pdb file name. \n' + \
                '>>> phenix.pdb.bioassembly_reconstruction "input_file_name.pdb"')
  elif (len(args) > 2):
    raise Sorry('To many input parameters are given')
  elif args[0] in ['--help','-help','-h']:
    print help(run)
  elif not os.path.isfile(args[0]):
    print 'There is no file {} in the current directory \n {}'.format(args[0],os.getcwd())
    raise Sorry('File name error')
  else:
    m = multimer(args[0])
    if len(args) == 2:
      # output file name given
      output_file_name = args[1]
    else:  
      # output file name not given
      output_file_name = 'copy_' + args[0]
    # write the bio-assembly reconstructed file in current directory
    m.write(output_file_name)
    print 'bio-assembly reconstructed pdb file was added to your current directory'
    print os.getcwd() + '/' + output_file_name

  
  

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
