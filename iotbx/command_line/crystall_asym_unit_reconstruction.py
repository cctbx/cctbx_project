from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.crystall_asym_unit_reconstruction

from iotbx import pdb
from libtbx.utils import Sorry
from  iotbx.pdb.multimer_reconstruction import multimer
import sys, os

def crystall_asym_unit_reconstruction(args, command_name="phenix.pdb.crystall_asym_unit_reconstruction"):
  ''' (string) -> file
  
  Reconstruction of the crystallographic asymmetric unit by applying
  MTRIX transformation, from the pdb file, to all chains that are not already in file.
  
  Produce a new pdb file containg the crystallographic asymmetric unit multimer information
  
  Arguments:
  args: a string containing the pdb file name, 'file_name.pdb' and an optional output file name
        if no file name is given, the output file will have the name 'crystall_asymmetric_unit_file_name.pdb'
        
  >>> phenix.pdb.crystall_asym_unit_reconstruction 'input_file_name.pdb'
  crystall_asymmetric_unit reconstructed pdb file was added to your current directory
  /path/crystall_asymmetric_unit_input_file_name.pdb
  
   >>> phenix.pdb.crystall_asym_unit_reconstruction 'input_file_name.pdb' 'output_file_name.pdb'
  crystall_asymmetric_unit reconstructed pdb file was added to your current directory
  /path/output_file_name.pdb
  
  '''
  if (len(args) == 0): 
    raise Sorry('No input filename is given. Please provide a pdb file name. \n' + \
                '>>> phenix.pdb.crystall_asym_unit_reconstruction "input_file_name.pdb"')
  elif ('--help' in args) or ('-h' in args):
    print help(crystall_asym_unit_reconstruction)
  elif (len(args) > 2):
    raise Sorry('To many input parameters are given')  
  elif not os.path.isfile(args[0]):
    print 'There is no file {} in the current directory \n {}'.format(args[0],os.getcwd())
    raise Sorry('File name error')
  else:
    m = multimer(args[0],'cau')
    # write the crystallographic asymmetric unit reconstructed file in current directory
    if len(args) == 2:
      # output file name given
      m.write(args[1])      
    else:  
      # output file name not given
      m.write()
       
    print 'crystallographic asymmetric unit reconstructed pdb file was added to your current directory'
    print os.getcwd() + '/' + m.pdb_output_file_name

  

if (__name__ == "__main__"):
  crystall_asym_unit_reconstruction(args=sys.argv[1:])
