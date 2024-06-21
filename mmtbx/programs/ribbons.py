from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
import iotbx.pdb
from pathlib import Path
from mmtbx.kinemage.ribbons import find_contiguous_protein_residues, find_contiguous_nucleic_acid_residues
from mmtbx.kinemage.ribbons import make_protein_guidepoints, make_nucleic_acid_guidepoints
from mmtbx.kinemage.ribbons import untwist_ribbon, swap_edge_and_face

version = "1.0.0"

master_phil_str = '''
do_protein = True
  .type = bool
  .short_caption = Construct protein ribbons
  .help = Construct ribbons for the protein sections of the model
do_nucleic_acid = True
  .type = bool
  .short_caption = Construct nucleic acid ribbons
  .help = Construct ribbons for the nucleic acid sections of the model
untwist_ribbons = True
  .type = bool
  .short_caption = Untwist ribbons
  .help = Remove excess twist from ribbons by making the neighboring dot products positive
DNA_style = False
  .type = bool
  .short_caption = DNA style ribbons
  .help = Use the DNA style ribbons instead of the default RNA style ribbons (rotates them by 90 degrees)
color_by = *rainbow secondary_structure
  .type = choice(multi=False)
  .short_caption = How to color the ribbons
  .help = How to color the ribbons
'''

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
mmtbx.ribbons version {}: Given PDB file produce Kinemage file with ribbons representation.

This program produces a kinemage file with ribbons representation of the input file using the
same algorithm as the Richardson lab's MolProbity server.  The code is based on the Java code
within the javadev repository.

How to run:
  mmtbx.ribbons model.pdb
  
Output:
  If neither output.file_name nor output.filename is specified, it will write
  to a file with the same name as the input model file name but with the
  extension replaced with with '.kin'.

'''.format(version)
  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

    if self.params.output.filename is None:
      # If the output file name is not specified, use the same root as the
      # input file and replace the suffix with .kin.
      suffix = '.kin'
      inName = self.data_manager.get_default_model_name()
      p = Path(inName)
      self.params.output.filename = str(p.with_suffix(suffix))
      print('Setting output.filename Phil parameter to',self.params.output.filename)

# ------------------------------------------------------------------------------

  def run(self):
    self.model = self.data_manager.get_model()
    structure = self.model.get_hierarchy()

    # Create the output string that will be written to the output file.  It will be filled in during the run.
    outString = ""
    
    if self.params.do_protein:
      # Find the contiguous protein residues by CA distance
      contiguous_residues = find_contiguous_protein_residues(structure)
      print('Found {} contiguous protein residue lists'.format(len(contiguous_residues)))

      # @todo

      for contig in contiguous_residues:
        guidepoints = make_protein_guidepoints(contig, structure)
        print(' Made {} protein guidepoints for {} residues'.format(len(guidepoints),len(contig)))
        if self.params.untwist_ribbons:
          print('  Untwisted ribbon')
          untwist_ribbon(guidepoints)

      # @todo
      
    if self.params.do_nucleic_acid:
      # Find the contiguous nucleic acid residues by CA distance
      contiguous_residues = find_contiguous_nucleic_acid_residues(structure)
      print('Found {} contiguous nucleic acid residue lists'.format(len(contiguous_residues)))

      # @todo

      for contig in contiguous_residues:
        guidepoints = make_nucleic_acid_guidepoints(contig, structure)
        print(' Made {} NA guidepoints for {} residues'.format(len(guidepoints),len(contig)))
        if self.params.untwist_ribbons:
          print('  Untwisted ribbon')
          untwist_ribbon(guidepoints)
        if self.params.DNA_style:
          print('  Swapped edge and face (DNA style)')
          swap_edge_and_face(guidepoints)
        else:
          print('  Using RNA style ribbons')

      # @todo

    # Write the output to the specified file.
    self.data_manager._write_text("Text", outString, self.params.output.filename)
