
# XXX this largely duplicates cif_as_pdb and pdb_as_cif, but I wanted
# something that would convert in both directions

from __future__ import division
import iotbx.phil
from libtbx.utils import Sorry
import os.path
import sys

master_phil = iotbx.phil.parse("""
convert
  .short_caption = mmCIF/PDB conversion
  .caption = This tool will convert between the two supported model formats: \
    files in PDB format will be output as mmCIF and vice-versa.  Note that \
    many of the tools in Phenix support both formats, and phenix.refine can \
    output mmCIF files directly if desired.
  .style = auto_align box caption_img:icons/custom/phenix.pdbtools.png
{
  file_name = None
    .type = path
    .short_caption = Model file
    .style = file_type:pdb bold input_file
  output_file = None
    .type = path
    .short_caption = Output file
    .style = new_file
}
""")

def run (args=(), params=None, out=sys.stdout) :
  assert (params.convert.file_name is not None)
  file_name = params.convert.file_name
  assert os.path.exists(file_name)
  from iotbx.pdb import pdb_input_from_any
  import iotbx.cif
  model = pdb_input_from_any(file_name=file_name)
  pdb_input = model.file_content()
  hierarchy = pdb_input.construct_hierarchy()
  basename = os.path.splitext(os.path.basename(file_name))[0]
  if (model.file_format == "pdb") :
    print >> out, "Converting %s to mmCIF format." % file_name
    cif_object = iotbx.cif.model.cif()
    if (params.convert.output_file is None) :
      params.convert.output_file = basename + ".cif"
    elif (params.convert.output_file.endswith(".pdb")) :
      raise Sorry("The output format is mmCIF, but the specified output "+
        "file name ends with .pdb.  Please use .cif or similar as the file "+
        "extension.")
    cif_object[basename] = hierarchy.as_cif_block(
      crystal_symmetry=pdb_input.crystal_symmetry())
    f = open(params.convert.output_file, "wb")
    print >> f, cif_object
    f.close()
  else :
    if (params.convert.output_file is None) :
      params.convert.output_file = basename + ".pdb"
    elif (params.convert.output_file.endswith(".cif")) :
      raise Sorry("The output format is PDB, but the specified output "+
        "file name ends with .cif.  Please use .pdb or similar as the file "+
        "extension.")
    f = open(params.convert.output_file, "wb")
    print >> f, hierarchy.as_pdb_string(
      crystal_symmetry=pdb_input.crystal_symmetry())
    f.close()
  return params.convert.output_file

def validate_params (params) :
  if (params.convert.file_name is None) :
    raise Sorry("No model file specified!")
  elif (params.convert.output_file is not None) :
    base, ext = os.path.splitext(params.convert.file_name)
    base2, ext2 = os.path.splitext(params.convert.output_file)
    if (ext == ext2) :
      raise Sorry("The file extensions for the input and output files are "+
        "identical, but the format will change.  Please use .pdb if you are "+
        "converting an mmCIF file, or .cif if converting a PDB file.")
  return True
