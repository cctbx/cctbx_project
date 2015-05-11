from __future__ import division

import iotbx.phil

# standard input params for find_secondary_structure

alpha_params=iotbx.phil.parse("""
     span = 3.5
       .type = float
       .help = number of residues included in i to i+span vector differences
       .short_caption = span for alpha structure

     buffer_residues = 0
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for alpha structure

     standard_length = 6
       .type = int
       .help = standard length for alpha structure
       .short_caption = Standard length for alpha structure

     minimum_length = 6
       .type = int
       .help = Minimum length for alpha structure
       .short_caption = Minimum length for alpha structure

     rise = 1.54
       .type = float
       .help = Alpha helix rise per reside
       .short_caption = Alpha helix rise

     minimum_overlap = 1
       .type = int
       .help = Minimum overlap at ends for alpha structure
       .short_caption = Minimum overlap for alpha structure

     rise_tolerance = 0.5
       .type = float
       .help = Tolerance in rise of helices in input file
       .short_caption = Tolerance in rise of helices

     dot_min_single = 0.3
       .type = float
       .help = Target dot product of i->i+1 with overall direction
       .short_caption = Target dot i to i+1

     dot_min =  0.9
       .type = float
       .help = minimum dot product of directions of average \
          of i-to-i+3  and i-to i+4 vectors to overall average direction

     allow_insertions = True
       .type = bool
       .help = Allow insertions in helices (adding residues)
       .short_caption= Allow insertions in helices

     allow_deletions = False
       .type = bool
       .help = Allow deletions in helices (adding residues)
       .short_caption= Allow deletions in helices

     base_score = 100.
       .type = float
       .help = Base score for helices
       .short_caption = Base score for helices

""")

beta_params=iotbx.phil.parse("""
     span = 2
       .type = float
       .help = number of residues included in i to i+span differences
       .short_caption = span for beta structure

     buffer_residues = 0
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for beta structure

     standard_length = 4
       .type = int
       .help = standard length for beta structure
       .short_caption = Standard length for beta structure

     minimum_length = 4
       .type = int
       .help = Minimum length for beta structure
       .short_caption = Minimum length for beta structure

     minimum_overlap = 1
       .type = int
       .help = Minimum overlap at ends for beta structure
       .short_caption = Minimum overlap for beta structure

     rise = 3.3
       .type = float
       .help = rise per reside for beta structure (3.2-3.4 is typical)
       .short_caption = rise per residue for beta structure

     rise_tolerance = 0.5
       .type = float
       .help = Tolerance in rise for beta structure in input file
       .short_caption = Tolerance in rise for beta structure

     dot_min_single = 0.5
       .type = float
       .help = Target dot product of i->i+1 with overall direction
       .short_caption = Target dot i to i+1

     dot_min =  0.75
       .type = float
       .help = minimum dot product of directions of \
          i-to-i+2 vectors to overall average direction

     allow_insertions = False
       .type = bool
       .help = Allow insertions in strands (adding residues)
       .short_caption= Allow insertions in strands

     allow_deletions = False
       .type = bool
       .help = Allow deletions in strands (adding residues)
       .short_caption= Allow deletions in strands

     base_score = 10.
       .type = float
       .help = Base score for strands
       .short_caption = Base score for strands

""")

other_params = iotbx.phil.parse("""
     span = None
       .type = float
       .help = None
       .short_caption = None

     buffer_residues = 2
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for other structure

     standard_length = 5
       .type = int
       .help = standard length for other structure
       .short_caption = Standard length for other structure

     minimum_length = 4
       .type = int
       .help = Minimum length for other structure
       .short_caption = Minimum length for other structure

     minimum_overlap = 1
       .type = int
       .help = Minimum overlap at ends for other structure
       .short_caption = Minimum overlap for other structure

     rise = None
       .type = float
       .help = None
       .short_caption = None

     rise_tolerance = None
       .type = float
       .help = None
       .short_caption = None

     dot_min_single = None
       .type = float
       .help = Target dot product of i->i+1 with overall direction
       .short_caption = Target dot i to i+1

     dot_min =  None
       .type = float
       .help = None

     allow_insertions = False
       .type = bool
       .help = Allow insertions in other (adding residues)
       .short_caption= Allow insertions in other

     allow_deletions = False
       .type = bool
       .help = Allow deletions in other (adding residues)
       .short_caption= Allow deletions in other

     base_score = 1.
       .type = float
       .help = Base score for other
       .short_caption = Base score for other
""")


master_phil = iotbx.phil.parse("""

  input_files {
    pdb_in = None
      .type = path
      .help = Input PDB file
      .short_caption = Input PDB file
  }

  find_ss_structure {

     find_alpha = True
       .type = bool
       .help = Find alpha helices
       .short_caption = Find alpha helices

     find_beta = True
       .type = bool
       .help = Find beta structure
       .short_caption = Find beta structure

     find_other = False
       .type = bool
       .help = Find other structure
       .short_caption = Find other structure

     exclude_alpha_in_beta  = True
       .type = bool
       .help = Exclude regions already identified as alpha
       .short_caption = Exclude alpha regions from beta

     exclude_alpha_beta_in_other  = True
       .type = bool
       .help = Exclude regions already identified as alpha or beta
       .short_caption = Exclude alpha and beta regions from other

     make_unique = True
       .type = bool
       .help = Assign each residue to a unique type of structure
       .short_caption = Assign residues to unique structure

  }
""")
