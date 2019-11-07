
# bkpoon (09/27/2016) changed to use cif_as_pdb and pdb_as_cif for more
# functionality

from __future__ import absolute_import, division, print_function
import iotbx.cif
import iotbx.phil
from iotbx.cli_parser import run_program
from iotbx.command_line import cif_as_pdb
from iotbx.pdb import pdb_input_from_any
from libtbx.utils import Sorry
from mmtbx.programs import pdb_as_cif
import os.path
import sys

master_phil = iotbx.phil.parse("""
convert
  .short_caption = mmCIF/PDB conversion
  .caption = This tool will convert between the two supported model formats: \
    files in PDB format will be output as mmCIF and vice-versa.  Note that \
    many of the tools in Phenix support both formats, and phenix.refine can \
    output mmCIF files directly if desired. The output filename will be the \
    same as the input filename, but the extension is changed (e.g. test.pdb \
    becomes test.cif and vice-versa).
  .style = auto_align box caption_img:icons/custom/phenix.pdbtools.png
{
  file_name = None
    .type = path
    .short_caption = Model file
    .style = file_type:pdb bold input_file
}
""")

def run(args=(), params=None, out=sys.stdout):
  assert (params.convert.file_name is not None)
  file_name = params.convert.file_name
  error_message = 'Unable to parse %s. ' % file_name
  error_message += 'Please make sure that the file exists and is a valid model in PDB or CIF format.'
  if (not os.path.exists(file_name)):
    raise Sorry(error_message)
  try:
    model = pdb_input_from_any(file_name=file_name)
  except Exception:
    raise Sorry(error_message)
  output_file = os.path.splitext(os.path.basename(file_name))[0]
  if (model.file_format == "pdb"):
    print("Converting %s to mmCIF format." % file_name, file=out)
    run_program(program_class=pdb_as_cif.Program, args=[file_name])
    output_file += '.cif'
  elif (model.file_format == 'cif'):
    print("Converting %s to PDB format." % file_name, file=out)
    cif_as_pdb.run([file_name])
    output_file += '.pdb'
  else:
    raise Sorry(error_message)
  if (not os.path.exists(output_file)):
    raise Sorry('Could not find output file. Please check that %s exists'
                % output_file)
  return output_file

def validate_params(params):
  if (params.convert.file_name is None):
    raise Sorry("No model file specified!")
  return True
