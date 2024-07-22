from __future__ import absolute_import, division, print_function
from turtle import window_width
from libtbx.program_template import ProgramTemplate
import iotbx.pdb
import mmtbx.secondary_structure
from pathlib import Path
import numpy as np
from scipy.interpolate import make_interp_spline
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
coil_width = 1.0
  .type = float
  .short_caption = Coil width
  .help = Width of the coil part of the ribbon
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

  def splineInterplate(self, pts, nIntervals):
    # Returns a list of interpolated points between the points using a spline interpolation.
    # Based on https://docs.scipy.org/doc/scipy/tutorial/interpolate/1D.html#parametric-spline-curves
    # @param pts: The guidepoints to interpolate between (a list of 3D points with at least 6 points in it).
    # The first and last points are duplciates to form knots, so are not interpolated between.
    # @param nIntervals: The number of intervals to interpolate between each pair of guidepoints, skipping the first and last.
    # @return: The list of interpolated points

    # Construct the parametric spline for a curve in 3D space.
    # We stick the end points in like any other points, but we don't use them in the interpolation.
    u = range(len(pts))
    x = [pt[0] for pt in pts]
    y = [pt[1] for pt in pts]
    z = [pt[2] for pt in pts]
    p = np.stack( (x, y, z) )
    spl = make_interp_spline(u, p, axis=1)

    # Compute the interpolated points
    count = (len(pts) - 1) * nIntervals + 1
    uu = np.linspace(0, len(pts) - 1, count)
    xx, yy, zz = spl(uu)

    ret = []
    for i in range(count):
      ret.append( (xx[i], yy[i], zz[i]) )
    return ret

# ------------------------------------------------------------------------------

  # Commented out the fields and did not include methods that we don't need for our implementation
  # Switched getType() to making type public.
  # Its fields should be based on the PDB secondary structure records.
  class Range:
    def __init__(self, pType = 'COIL', pInit = 0, pEnd = 0, pSense = 0, pStrand = 1, pSheet = ' ', pFlipped = False):
      # self._rangeIndex = 0
      self.type = pType
      #self._chainId = ''
      self.initSeqNum = pInit
      self.endSeqNum = pEnd
      #self._initICode = ' '
      #self._endICode = ' '
      self.sense = pSense    # 0 if first strand, 1 if parallel, -1 if anti-parallel
      self.strand = pStrand   # Starts at 1 for each strand within a sheet and increases by 1
      self.sheetID = pSheet
      self.flipped = pFlipped  # Used by printing code to make the tops within a set of ribbons all point the same way
      #self.previous = None
      #self.next = None
      #self.duplicateOf = None

  class RibbonElement:
    def __init__(self, other=None, setRange=None):
      self.start = 0
      self.end = 0
      self.type = 'COIL'
      self.range = setRange
      if other is not None:
        self.start = other.start
        self.end = other.end
        self.type = other.type
        self.range = other.range
      if self.range is not None:
        if self.range.type == 'turn':
          self.type = 'COIL'
        else:
          self.type = self.range.type

    def sameSSE(self, that):
      return self.type == that.type and self.range == that.range

    def like(self, that):
      return self.sameSSE(that) and self.start == that.start and self.end == that.end

    def __lt__(self, that):
      a = 0
      b = 0
      if self.range is not None:
        a = self.range.strand
      if that.range is not None:
        b = that.range.strand
      return a < b

  def printFancyRibbon(self, guides, widthAlpha, widthBeta, listAlpha, listBeta, listCoil, listCoilOutline, doRainbow, rnaPointIDs):
    # Constructs a kinemage string for a ribbon representation of a protein or nucleic acid chain.
    # Makes a triangulated ribbon with arrowheads, etc.
    # @param guides: The guidepoints for the ribbon
    # @param widthAlpha: The width of the alpha helix part of the ribbon
    # @param widthBeta: The width of the beta sheet part of the ribbon
    # @param listAlpha: The kinemage string describing the alpha helix part of the ribbon
    # @param listBeta: The kinemage string describing the beta sheet part of the ribbon
    # @param listCoil: The kinemage string describing the coil part of the ribbon
    # @param listCoilOutline: The kinemage string describing the outline of the coil part of the ribbon
    # @param doRainbow: Whether or not to color the ribbon by rainbow
    # @param rnaPointIDs: If true, residue names in point IDs will be decided based on the RNA/DNA alignment of guidepoints to residues.
    # If false (the default), the protein style will be used instead.
    # @return: The kinemage string for the ribbon representation

    ret = ""

    # Initialize local references
    widthCoil = self.params.coil_width
    secStruct = self.secondaryStructure

    # Seven strands of guidpoints: coil, +/-alpha, +/-beta, +/-beta arrowheads.
    # Each strand is a list of 3D points and there is a strand for each offset from the center spline.
    halfWidths = [0.0, -widthAlpha/2, widthAlpha/2, -widthBeta/2, widthBeta/2, -widthBeta, widthBeta]
    strands = [[] for _ in range(len(halfWidths))]
    for g in guides:
      for i in range(len(halfWidths)):
        strands[i].append(g.pos + g.dvec * halfWidths[i])

    # Seven strands of interpolated points
    splinepts = []
    for i in range(len(strands)):
      splinepts.append(self.splineInterplate(strands[i], self.nIntervals))

    # Discovery of ribbon elements: ribbons, ropes, and arrows.
    # We skip the regions associated with the first and last points, which are present to cause the
    # curve to pass through those points (like knots) but are not interpolated between.
    # We do that by looping through fewer entries and by adding 1 to the index.
    # @todo Figure out what this code is doing with its secondary structure class and reproduce with ours.  
    for i in range(len(guides) - 1 - 2):
      g1 = guides[i + 1]
      g2 = guides[i + 2]


    # @todo

    return ret

# ------------------------------------------------------------------------------

  def run(self):
    self.model = self.data_manager.get_model()
    self.nIntervals = 4

    # Analyze the secondary structure and make a dictionary that maps from residue sequence number to secondary structure type
    # by filling in 'COIL' as a default value for each and then parsing all of the secondary structure records in the
    # model and filling in the relevant values for them.
    # @todo Can pass it an atom selection cache, so maybe we do this after we have selected the atoms we want to use?
    print('Finding secondary structure:')
    ss_manager = mmtbx.secondary_structure.manager(self.model.get_hierarchy())
    self.secondaryStructure = {}
    for model in self.model.get_hierarchy().models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          self.secondaryStructure[residue_group.resseq] = 'COIL'
    for line in ss_manager.records_for_pdb_file().splitlines():
      r = self.Range(line[0:5].strip())
      if r.type == 'HELIX':
        r.sheetID = line[11:14].strip()
        r.initSeqNum = int(line[21:26].strip())
        r.endSeqNum = int(line[33:38].strip())
      elif r.type == 'SHEET':  # Java code marks this internally as STRAND but we leave it as SHEET
        r.sheetID = line[11:14].strip()
        r.initSeqNum = int(line[23:28].strip())
        r.endSeqNum = int(line[33:38].strip())
        r.sense = int(line[38:40].strip())
        r.strand = int(line[7:10].strip())
      elif r.type == 'TURN':
        # @todo Test this code on a file that actually has turns to see if it works
        # In fact, we turn turns int coils, so we really don't care.
        # r.sheetID = line[11:14].strip()
        r.initSeqNum = int(line[22:26].strip())
        r.endSeqNum = int(line[33:37].strip())
      for i in range(r.initSeqNum,r.endSeqNum+1):
        self.secondaryStructure[i] = r

    # The name of the structure, which is the root of the input file name (no path, no extension)
    self.idCode = self.data_manager.get_default_model_name().split('.')[0]
    if self.idCode is None or self.idCode == "":
      self.idCode = "macromol"

    # Round-robin colors for the backbone of the ribbons for when we are coloring each chain or model.
    self.backBoneColors = [ "white", "yellowtint", "peachtint", "pinktint", "lilactint", "bluetint", "greentint" ]

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
          contiguous_residue_lists = find_contiguous_protein_residues(chain)
          print('Found {} contiguous protein residue lists'.format(len(contiguous_residue_lists)))

          if len(contiguous_residue_lists) > 0:
            if groupByModel:
              outString += "@subgroup {{chain{}}} dominant master= {chain {}}}\n".format(chain.id, chain.id)
            else:
              outString += "@group {{{} {}}} dominant\n".format(self.idCode, chain.id)
            outString += "@subgroup {ribbon}\n"

            bbColor = chainColors[chain.id]
            if bbColor == "white":
              # Distinguish between the different types of secondary structure for the first chain
              outString += "@colorset {{alph{}}} red\n".format(chain.id)
              outString += "@colorset {{beta{}}} lime\n".format(chain.id)
              outString += "@colorset {{coil{}}} white\n".format(chain.id)
            else:
              # Do all secondary structure in the same color for all but the first chain to clean up the display
              outString += "@colorset {{alph{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{beta{}}} {}\n".format(chain.id, bbColor)
              outString += "@colorset {{coil{}}} {}\n".format(chain.id, bbColor)

            for contig in contiguous_residue_lists:
              guidepoints = make_protein_guidepoints(contig)
              print(' Made {} protein guidepoints for {} residues'.format(len(guidepoints),len(contig)))
              if self.params.untwist_ribbons:
                print('  Untwisted ribbon')
                untwist_ribbon(guidepoints)
              # There is always secondary structure looked up for protein residues, so we skip the case from the Java code
              # where it can be missing.
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
          contiguous_residue_lists = find_contiguous_nucleic_acid_residues(chain)
          print('Found {} contiguous nucleic acid residue lists'.format(len(contiguous_residue_lists)))

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
