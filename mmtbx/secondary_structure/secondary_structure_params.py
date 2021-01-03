from __future__ import absolute_import, division, print_function

import iotbx.phil

# standard input params for find_secondary_structure

alpha_params=iotbx.phil.parse("""

     name = Alpha helix
       .type = str
       .help = Secondary structure name (helix)

     span = 3.5
       .type = float
       .help = number of residues included in i to i+span vector differences
       .short_caption = span for alpha structure

     buffer_residues = 0
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for alpha structure

     n_link_min = 0
       .type = int
       .help = Minimum number of residues linking adjacent segments \
           in H-bonded structure (i.e., HELIX/SHEET records)
       .short_caption = Minimum link residues between segments

     standard_length = 6
       .type = int
       .help = standard length for alpha structure
       .short_caption = Standard length for alpha structure

     minimum_length = 6
       .type = int
       .help = Minimum length for alpha structure
       .short_caption = Minimum length for alpha structure

     residues_per_turn = 3.6
       .type = float
       .help = Alpha helix residues per turn
       .short_caption = Alpha helix residues per turn

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

     target_i_ip3 = None
       .type = float
       .help = None (target i->i+3 distance)
       .short_caption = none

     tol_i_ip3 = None
       .type = float
       .help = None (tolerance in target i->i+3 distance)
       .short_caption = none

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

three_ten_params=iotbx.phil.parse("""
     name = 3-10 helix
       .type = str
       .help = Secondary structure name (alpha helix)

     span = 3
       .type = float
       .help = number of residues included in i to i+span vector differences.\
               This is for three-ten helices
       .short_caption = span for three-ten structure

     buffer_residues = 0
       .type = int
       .help = number of residues included on ends of segments for three-ten
       .short_caption = buffer residues for three-10 structure

     n_link_min = 0
       .type = int
       .help = Minimum number of residues linking adjacent segments \
           in H-bonded structure (i.e., HELIX/SHEET records)
       .short_caption = Minimum link residues between segments

     standard_length = 6
       .type = int
       .help = standard length for three_ten structure
       .short_caption = Standard length for three_ten structure

     minimum_length = 6
       .type = int
       .help = Minimum length for three_ten structure
       .short_caption = Minimum length for three_ten structure

     residues_per_turn = 3
       .type = float
       .help = Three-ten helix residues per turn
       .short_caption = Three-ten helix residues per turn

     rise = 2.0
       .type = float
       .help = Three ten helix rise per reside
       .short_caption = Three ten helix rise

     minimum_overlap = 1
       .type = int
       .help = Minimum overlap at ends for three_ten structure
       .short_caption = Minimum overlap for three_ten structure

     rise_tolerance = 0.5
       .type = float
       .help = Tolerance in rise of helices in input file
       .short_caption = Tolerance in rise of helices

     target_i_ip3 = None
       .type = float
       .help = None (target i->i+3 distance)
       .short_caption = none

     tol_i_ip3 = None
       .type = float
       .help = None (tolerance in target i->i+3 distance)
       .short_caption = none

     dot_min_single = 0.5
       .type = float
       .help = Target dot product of i->i+1 with overall direction
       .short_caption = Target dot i to i+1

     dot_min =  0.9
       .type = float
       .help = minimum dot product of directions of average \
          of i-to-i+3  and i-to i+4 vectors to overall average direction

     allow_insertions = False
       .type = bool
       .help = Allow insertions in helices (adding residues)
       .short_caption= Allow insertions in helices

     allow_deletions = False
       .type = bool
       .help = Allow deletions in helices (adding residues)
       .short_caption= Allow deletions in helices

     base_score = 1.
       .type = float
       .help = Base score for 3_10 helices
       .short_caption = Base score for 3_10 helices

""")

pi_params=iotbx.phil.parse("""
     name = Pi helix
       .type = str
       .help = Secondary structure name (pi helix)

     span = 4.0
       .type = float
       .help = number of residues included in i to i+span vector differences
       .short_caption = span for pi structure

     buffer_residues = 0
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for pi structure

     n_link_min = 0
       .type = int
       .help = Minimum number of residues linking adjacent segments \
           in H-bonded structure (i.e., HELIX/SHEET records)
       .short_caption = Minimum link residues between segments

     standard_length = 6
       .type = int
       .help = standard length for pi structure
       .short_caption = Standard length for pi structure

     minimum_length = 6
       .type = int
       .help = Minimum length for pi structure
       .short_caption = Minimum length for pi structure

     residues_per_turn = 4.1
       .type = float
       .help = Pi helix residues per turn
       .short_caption = Pi helix residues per turn

     rise = 0.95
       .type = float
       .help = Pi helix rise per reside
       .short_caption = Pi helix rise

     minimum_overlap = 1
       .type = int
       .help = Minimum overlap at ends for pi structure
       .short_caption = Minimum overlap for pi structure

     rise_tolerance = 0.5
       .type = float
       .help = Tolerance in rise of pi helices
       .short_caption = Tolerance in rise of pi helices

     target_i_ip3 = None
       .type = float
       .help = None (target i->i+3 distance)
       .short_caption = none

     tol_i_ip3 = None
       .type = float
       .help = None (tolerance in target i->i+3 distance)
       .short_caption = none

     dot_min_single = 0.1
       .type = float
       .help = Target dot product of i->i+1 with overall direction for pi
       .short_caption = Target dot i to i+1 for pi

     dot_min =  0.9
       .type = float
       .help = minimum dot product of directions of average \
          of i-to-i+3  and i-to i+4 vectors to overall average direction

     allow_insertions = False
       .type = bool
       .help = Allow insertions in helices (adding residues)
       .short_caption= Allow insertions in helices

     allow_deletions = False
       .type = bool
       .help = Allow deletions in helices (adding residues)
       .short_caption= Allow deletions in helices

     base_score = 1.
       .type = float
       .help = Base score for pi helices
       .short_caption = Base score for pi helices

""")

beta_params=iotbx.phil.parse("""

     name = Beta strand
       .type = str
       .help = Secondary structure name (strand)

     span = 2
       .type = float
       .help = number of residues included in i to i+span differences
       .short_caption = span for beta structure

     buffer_residues = 0
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for beta structure

     n_link_min = 3
       .type = int
       .help = Minimum number of residues linking adjacent segments \
           in H-bonded structure (i.e., HELIX/SHEET records)
       .short_caption = Minimum link residues between segments

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

     target_i_ip3 = 10.
       .type = float
       .help = Target i->i+3 distance
       .short_caption = Target i to i+3 distance

     tol_i_ip3 = 1.5
       .type = float
       .help = Tolerance in target i->i+3 distance
       .short_caption = Tolerance in target i to i+3 distance

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

     max_sheet_ca_ca_dist = 6.
       .type = float
       .help = Maximum CA-CA distance between paired strands in sheets
       .short_caption = Max CA-CA distance between strands in sheets

     tolerant_max_sheet_ca_ca_dist = 8.
       .type = float
       .help = Tolerant maximum CA-CA distance between paired strands \
           in sheets
       .short_caption = Tolerant max CA-CA distance between strands \
              in sheets

     min_sheet_length = 4
       .type = int
       .help = Minimum H-bonded segment to include in a sheet
       .short_caption = Min length of segment in a sheet

     tolerant_min_sheet_length = 2
       .type = int
       .help = Tolerant minimum H-bonded segment to include in a sheet
       .short_caption = Tolerant min length of segment in a sheet

     allow_ca_only_model = True
       .type = bool
       .help = Allow the use of CA-only models, defining H-bonding for \
               N and O as if they were present
       .short_caption= Allow CA-only models



""")

other_params = iotbx.phil.parse("""
     name = Other structure
       .type = str
       .help = Secondary structure name (other)

     span = None
       .type = float
       .help = None
       .short_caption = None

     buffer_residues = 2
       .type = int
       .help = number of residues included on ends of segments
       .short_caption = buffer residues for other structure

     n_link_min = 0
       .type = int
       .help = Minimum number of residues linking adjacent segments \
           in H-bonded structure (i.e., HELIX/SHEET records)
       .short_caption = Minimum link residues between segments

     standard_length = 8
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

     target_i_ip3 = None
       .type = float
       .help = None (target i->i+3 distance)
       .short_caption = none

     tol_i_ip3 = None
       .type = float
       .help = None (tolerance in target i->i+3 distance)
       .short_caption = none

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
