"""Construct ribbons representing a model"""
from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
import iotbx.pdb
import mmtbx.secondary_structure
from pathlib import Path
from mmtbx.kinemage.validation import get_chain_color
from mmtbx.kinemage.ribbons import (
    find_contiguous_protein_residues, find_contiguous_nucleic_acid_residues,
    make_protein_guidepoints, make_nucleic_acid_guidepoints,
    untwist_ribbon, swap_edge_and_face, _IsNucleicAcidResidue,
    chain_has_DNA, chain_has_RNA)
from mmtbx.kinemage.ribbon_rendering import (
    Range, build_secondary_structure_map, consolidate_sheets,
    render_fancy_ribbon, generate_chain_ribbons)

version = "1.2.1"

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
color_by = *rainbow secondary_structure solid
  .type = choice(multi=False)
  .short_caption = How to color the ribbons
  .help = How to color the ribbons. Rainbow adjusts the color smoothly along each chain. Solid make each chain a single color. Secondary structure colors every 7th chain by the secondary structure and makes the others solid colors.
do_plain_coils = False
  .type = bool
  .short_caption = Do plain coils
  .help = Do plain coils (no halo) even for the helix and sheet parts of the protein
coil_width = 1.0
  .type = float
  .short_caption = Coil width
  .help = Width of the coil part of the ribbon
selection = (altloc ' ' or altloc '' or altloc a)
  .type = atom_selection
  .short_caption = Atom selection
  .help = Atom selection description
nucleic_acid_as_helix = True
  .type = bool
  .short_caption = Draw nucleic acids as helix
  .help = If true, draw nucleic acids as helix rather than treating as coil because there are not secondary structure records for them
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
  If output.filename is not specified, it will write
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
      print('Setting output.filename Phil parameter to',self.params.output.filename, file=self.logger)

# ------------------------------------------------------------------------------

  def run(self):
    self.model = self.data_manager.get_model()

    # Apply a selection to the model if one is provided.  The default is to select atoms in the
    # first alternate location.
    hierarchy = self.model.get_hierarchy()
    selection_string = self.params.selection
    if selection_string is None:
      selection_string = "(altloc ' ' or altloc '' or altloc a)"
    selection = hierarchy.atom_selection_cache().selection(selection_string)
    hierarchy = hierarchy.select(selection)

    # See if the model file has secondary structure records.
    # This should return None if there are no secondary structure records in the model.
    sec_str_from_pdb_file = self.model.get_ss_annotation()

    # Build the secondary structure map using the annotation API directly
    print('Finding secondary structure:', file=self.logger)
    params = mmtbx.secondary_structure.manager.get_default_ss_params()
    params.secondary_structure.protein.search_method = "from_ca"
    params = params.secondary_structure
    ss_manager = mmtbx.secondary_structure.manager(hierarchy,
                                                   params=params,
                                                   sec_str_from_pdb_file=sec_str_from_pdb_file,
                                                   log=self.logger)
    annotation = ss_manager.actual_sec_str
    self.secondaryStructure = build_secondary_structure_map(
        hierarchy, annotation=annotation, log=self.logger)
    consolidate_sheets(self.secondaryStructure)

    # If we are treating nucleic acids as helices, change the record of each nucleic acid residue to be a helix.
    if self.params.nucleic_acid_as_helix:
      r = Range('HELIX')
      for res in hierarchy.residue_groups():
        if _IsNucleicAcidResidue(res.unique_resnames()[0]):
          self.secondaryStructure[res.resseq_as_int()] = r

    # The name of the structure, which is the root of the input file name (no path, no extension)
    self.idCode = self.data_manager.get_default_model_name().split('.')[0]
    if self.idCode is None or self.idCode == "":
      self.idCode = "macromol"

    # Create the output string that will be written to the output file.  It will be filled in during the run.
    outString = ""

    # Fill in the header information
    outString += "@kinemage 1\n"
    outString += "@onewidth\n"

    # Handle multiple models
    groupByModel = hierarchy.models_size() > 1
    for model in hierarchy.models():
      modelID = model.id
      if modelID == "":
        modelID = "_"
      print('Processing model', modelID, 'with', len(model.chains()), 'chains', file=self.logger)
      if groupByModel:
        outString += "@group {{{} {}}} animate dominant master= {{all models}}\n".format(self.idCode, str(modelID).strip())

      # Make a list of all the chain names in the model with only one entry per name.
      # Use this to make a dictionary to look up the color that is used for each chain name.
      # This ensures that the chains are colored the same no matter their order or repeat in the file.
      chainNames = set()
      for chain in model.chains():
        chainNames.add(chain.id)
      chainNames = sorted(list(chainNames))
      chainCount = 0
      chainColors = {}
      for name in chainNames:
        # Backbone color by model ID if we have multiple models, or by chain ID if we have multiple chains in a model.
        if groupByModel:
          c = get_chain_color(modelID)
        else:
          c = get_chain_color(chainCount)
        chainColors[name] = c
        chainCount += 1

      # Determine whether DNA, RNA, or both are present in the model
      hasDNA = any(chain_has_DNA(chain) for chain in model.chains())
      hasRNA = any(chain_has_RNA(chain) for chain in model.chains())

      # Cycle over all chains in the model and make a group or subgroup for each chain
      # depending on whether we are grouping by model or not.
      for chain in model.chains():
        print('Processing chain', chain.id, file=self.logger)

        if self.params.do_protein:
          # Find the contiguous protein residues by CA distance
          contiguous_residue_lists = find_contiguous_protein_residues(chain)
          print('Found {} contiguous protein residue lists'.format(len(contiguous_residue_lists)), file=self.logger)

          if len(contiguous_residue_lists) > 0:
            if groupByModel:
              outString += "@subgroup {{chain{}}} dominant master= {{chain {}}}\n".format(chain.id, chain.id)
            else:
              outString += "@group {{{} {}}} dominant\n".format(self.idCode, chain.id)
            outString += "@subgroup {ribbon}\n"

            # Set the basic colors for the parts of the ribbon, which may be overridden by rainbow coloring on a
            # per-residue basis.
            bbColor = chainColors[chain.id]
            if bbColor == "white" and not (self.params.color_by == "solid"):
              outString += "@colorset {{alph{}}} red\n".format(chain.id)
              outString += "@colorset {{beta{}}} lime\n".format(chain.id)
              outString += "@colorset {{coil{}}} white\n".format(chain.id)
            else:
              outString += "@colorset {{alph{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{beta{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{coil{}}} {}\n".format(chain.id, bbColor)

            for contig in contiguous_residue_lists:
              guidepoints = make_protein_guidepoints(contig)
              print(' Made {} protein guidepoints for {} residues'.format(len(guidepoints), len(contig)), file=self.logger)
              if self.params.untwist_ribbons:
                print('  Untwisted ribbon', file=self.logger)
                untwist_ribbon(guidepoints)
              if self.params.do_plain_coils:
                outString += render_fancy_ribbon(guidepoints, self.secondaryStructure,
                      width_alpha=2, width_beta=2.2,
                      list_alpha="color= {alph"+chain.id+"} master= {protein ribbon} master= {alpha}",
                      list_beta="color= {beta"+chain.id+"} master= {protein ribbon} master= {beta}",
                      list_coil="width= 4 color= {coil"+chain.id+"} master= {protein ribbon} master= {coil}",
                      do_rainbow=self.params.color_by == "rainbow")
              else:
                outString += render_fancy_ribbon(guidepoints, self.secondaryStructure,
                      width_alpha=2, width_beta=2.2,
                      list_alpha="color= {alph"+chain.id+"} master= {protein ribbon} master= {alpha}",
                      list_beta="color= {beta"+chain.id+"} master= {protein ribbon} master= {beta}",
                      list_coil="width= 4 fore color= {coil"+chain.id+"} master= {protein ribbon} master= {coil}",
                      list_coil_outline="width= 6 rear color= deadblack master= {protein ribbon} master= {coil}",
                      do_rainbow=self.params.color_by == "rainbow")

        if self.params.do_nucleic_acid:
          # Find the contiguous nucleic acid residues by CA distance
          contiguous_residue_lists = find_contiguous_nucleic_acid_residues(chain)
          print('Found {} contiguous nucleic acid residue lists'.format(len(contiguous_residue_lists)), file=self.logger)

          if len(contiguous_residue_lists) > 0:
            if groupByModel:
              outString += "@subgroup {{chain{}}} dominant master= {chain {}}}\n".format(chain.id, chain.id)
            else:
              outString += "@group {{{} {}}} dominant\n".format(self.idCode, chain.id)
            outString += "@subgroup {ribbon}\n"

            bbColor = chainColors[chain.id]
            if bbColor == "white":
              outString += "@colorset {{nucl{}}} lime\n".format(chain.id)
              outString += "@colorset {{ncoi{}}} white\n".format(chain.id)
            else:
              outString += "@colorset {{nucl{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{ncoi{}}} {}\n".format(chain.id, bbColor)

            for contig in contiguous_residue_lists:
              guidepoints = make_nucleic_acid_guidepoints(contig)
              print(' Made {} NA guidepoints for {} residues'.format(len(guidepoints), len(contig)), file=self.logger)
              if self.params.untwist_ribbons:
                print('  Untwisted ribbon', file=self.logger)
                untwist_ribbon(guidepoints)
              if self.params.DNA_style or (hasDNA and hasRNA and chain_has_DNA(chain)):
                print('  Swapped edge and face (DNA style)', file=self.logger)
                swap_edge_and_face(guidepoints)
              else:
                print('  Using RNA style ribbons', file=self.logger)

              outString += render_fancy_ribbon(guidepoints, self.secondaryStructure,
                    width_alpha=3.0, width_beta=3.0,
                    list_alpha="color= {nucl"+chain.id+"} master= {NA ribbon} master= {RNA helix?}",
                    list_beta="color= {nucl"+chain.id+"} master= {NA ribbon} master= {A-form}",
                    list_coil="width= 4 color= {ncoi"+chain.id+"} master= {NA ribbon} master= {coil}",
                    do_rainbow=self.params.color_by == "rainbow")

    # Write the output to the specified file.
    self.data_manager._write_text("Text", outString, self.params.output.filename)
