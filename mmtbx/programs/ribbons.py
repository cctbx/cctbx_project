from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
import iotbx.pdb
import mmtbx.secondary_structure
from pathlib import Path
from mmtbx.kinemage.ribbons import find_contiguous_protein_residues, find_contiguous_nucleic_acid_residues
from mmtbx.kinemage.ribbons import make_protein_guidepoints, make_nucleic_acid_guidepoints
from mmtbx.kinemage.ribbons import untwist_ribbon, swap_edge_and_face, non_CA_atoms_present

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
do_het_atoms = True
  .type = bool
  .short_caption = Make ribbons for HETATM residues
  .help = Construct ribbons for HETATM residue sections of the model
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
do_plain_coils = False
  .type = bool
  .short_caption = Do plain coils
  .help = Do plain coils (no halo) even for the helix and sheet parts of the protein
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

  def printFancyRibbon(self, guides, widthAlpha, widthBeta, listAlpha, listBeta, listCoil, listCoilOutline, doRainbow, rnaPointIDs):
    ret = ""

    secStruct = self.secondaryStructure
    # @todo
    return ret

# ------------------------------------------------------------------------------

  def run(self):
    self.model = self.data_manager.get_model()

    # Analyze the secondary structure and make a dictionary that maps from residue sequnce number to secondary structure type
    # by filling in 'coil' as a default value for each and then parsing all of the secondary structure records in the
    # model and filling in the relevant values for them.
    # @todo Maybe we only mark the residues just past one end of a helix or sheet as coil, not sure what to mark others.
    # @todo Can pass it an atom selection cache, so maybe we do this after we have selected the atoms we want to use?
    print('Finding secondary structure:')
    ss_manager = mmtbx.secondary_structure.manager(self.model.get_hierarchy())
    self.secondaryStructure = {}
    for model in self.model.get_hierarchy().models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          self.secondaryStructure[residue_group.resseq] = 'coil'
    for line in ss_manager.records_for_pdb_file().splitlines():
      rType = line[0:5].strip()
      if rType == 'HELIX':
        id = line[11:14].strip()
        firstResID = int(line[21:26].strip())
        lastResID = int(line[33:38].strip())
        helixID = int(id) - 1
        helixType = ss_manager.get_helix_types()[helixID]
        for i in range(firstResID,lastResID+1):
          self.secondaryStructure[i] = 'alpha'
      elif rType == 'SHEET':  # Java code marks this internally as STRAND but we leave it as SHEET
        id = line[11:14].strip()
        firstResID = int(line[23:28].strip())
        lastResID = int(line[33:38].strip())
        for i in range(firstResID,lastResID+1):
          self.secondaryStructure[i] = 'beta'
      elif rType == 'TURN':
        # @todo Test this code on a file that actually has turns to see if it works
        id = line[11:14].strip()
        firstResID = int(line[22:26].strip())
        lastResID = int(line[33:37].strip())
        for i in range(firstResID,lastResID+1):
          self.secondaryStructure[i] = 'turn'

    # Strip out the hetatoms if we are not doing them
    # @todo

    # The name of the structure, which is the root of the input file name (no path, no extension)
    self.idCode = self.data_manager.get_default_model_name().split('.')[0]
    if self.idCode is None or self.idCode == "":
      self.idCode = "macromol"

    # Round-robin colors for the backbone of the ribbons for when we are coloring each chain or model.
    self.backBoneColors = [ "white", "yellowtint", "peachtint", "greentint", "bluetint", "lilactint", "pinktint" ]

    # Rainbow colors when we are coloring by rainbow
    self.rainbowColors = [ "blue", "sky", "cyan", "sea", "green", "lime", "yellow" ,"gold" ,"orange" ,"red" ]

    # Create the output string that will be written to the output file.  It will be filled in during the run.
    outString = ""

    # Fill in the header information
    outString += "@kinemage 1\n"
    outString += "@onewidth\n"
    
    # Handle multiple models
    groupByModel = self.model.get_hierarchy().models_size() > 1
    for model in self.model.get_hierarchy().models():
      modelID = model.id
      if modelID == "":
        modelID = "_"
      print('Processing model', modelID, 'with', len(model.chains()), 'chains')
      if groupByModel:
        # @todo Test and figure out what to fill in for the "@todo" name based on what the original code produces
        outString += "@group {{{} {}}} dominant master= {{all models}}\n".format(self.idCode, "@todo")
 
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
          c = self.backBoneColors[modelID % len(self.backBoneColors)]
        else:
          c = self.backBoneColors[chainCount % len(self.backBoneColors)]
        chainColors[name] = c
        chainCount += 1

      # Cycle over all chains in the model and make a group or subgroup for each chain
      # depending on whether we are grouping by model or not.
      for chain in model.chains():
        print('Processing chain',chain.id)

        if self.params.do_protein:
          # Find the contiguous protein residues by CA distance
          contiguous_residues = find_contiguous_protein_residues(chain)
          print('Found {} contiguous protein residue lists'.format(len(contiguous_residues)))

          if len(contiguous_residues) > 0:
            if groupByModel:
              outString += "@subgroup {{chain{}}} dominant master= {chain {}}}\n".format(chain.id, chain.id)
            else:
              outString += "@group {{{} {}}} dominant\n".format(self.idCode, chain.id)
            outString += "@subgroup {ribbon}\n"

            bbColor = chainColors[chain.id]
            if bbColor == "white":
              outString += "@colorset {{alph{}}} red\n".format(chain.id)
              outString += "@colorset {{beta{}}} lime\n".format(chain.id)
              outString += "@colorset {{coil{}}} white\n".format(chain.id)
            else:
              outString += "@colorset {{alph{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{beta{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{coil{}}} {}\n".format(chain.id, bbColor)

            for contig in contiguous_residues:
              guidepoints = make_protein_guidepoints(contig)
              print(' Made {} protein guidepoints for {} residues'.format(len(guidepoints),len(contig)))
              if self.params.untwist_ribbons:
                print('  Untwisted ribbon')
                untwist_ribbon(guidepoints)
              # There is always secondary structure looked up for protein residues, so we skip the missing case from the Java code.
              if self.params.do_plain_coils:
                outString += self.printFancyRibbon(guidepoints, 2, 2.2,
                      "color= {alph"+chain.id+"} master= {protein} master= {ribbon} master= {alpha}",
                      "color= {beta"+chain.id+"} master= {protein} master= {ribbon} master= {beta}",
                      "width= 4 color= {coil"+chain.id+"} master= {protein} master= {ribbon} master= {coil}",
                      self.params.color_by == "rainbow",
                      False);
              else:
                outString += self.printFancyRibbon(guidepoints, 2, 2.2,
                      "color= {alph"+chain.id+"} master= {protein} master= {ribbon} master= {alpha}",
                      "color= {beta"+chain.id+"} master= {protein} master= {ribbon} master= {beta}",
                      "width= 4 fore color= {coil"+chain.id+"} master= {protein} master= {ribbon} master= {coil}",
                      "width= 6 rear color= deadblack master= {protein} master= {ribbon} master= {coil}",
                      self.params.color_by == "rainbow",
                      False)
      
        if self.params.do_nucleic_acid:
          # Find the contiguous nucleic acid residues by CA distance
          contiguous_residues = find_contiguous_nucleic_acid_residues(chain)
          print('Found {} contiguous nucleic acid residue lists'.format(len(contiguous_residues)))

          if len(contiguous_residues) > 0:
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

            for contig in contiguous_residues:
              guidepoints = make_nucleic_acid_guidepoints(contig)
              print(' Made {} NA guidepoints for {} residues'.format(len(guidepoints),len(contig)))
              if self.params.untwist_ribbons:
                print('  Untwisted ribbon')
                untwist_ribbon(guidepoints)
              if self.params.DNA_style:
                print('  Swapped edge and face (DNA style)')
                swap_edge_and_face(guidepoints)
              else:
                print('  Using RNA style ribbons')

              outString += self.printFancyRibbon(guidepoints, 3.0, 3.0,
                    "color= {nucl"+chain.id+"} master= {nucleic acid} master= {ribbon} master= {RNA helix?}",
                    "color= {nucl"+chain.id+"} master= {nucleic acid} master= {ribbon} master= {A-form}",
                    "width= 4 color= {ncoi"+chain.id+"} master= {nucleic acid} master= {ribbon} master= {coil}",
                    None,
                    self.params.color_by == "rainbow",
                    True)

          # @todo

    # Write the output to the specified file.
    self.data_manager._write_text("Text", outString, self.params.output.filename)
