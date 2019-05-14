from __future__ import absolute_import, division, print_function
inner_phil_str = """\
  scanbox_windows = 101 51 51
    #.type = ints(size_min=1, size_max=3, value_min=10)
    #  future: variable number of window passes
    .type = ints(size=3, value_min=10)
    .help = "Integer scanbox sizes for calculating background,"
            "for cycles 1,2, and 3, respectively."
            "Program defaults are 101, 51, and 51 pixels."
  peripheral_margin = 20
    .type = int(value_min=0)
    .help = "No spot detection inside margin; width in pixels."
"""

phil_str = """\
spotfinder {
%s
}
""" % inner_phil_str

labelit_related_commands = """\

#Preparation

wedgelimit = 2
  .type=int
  .help="maximum number of images to use for labelit indexing"
  .expert_level=2

goniometer_rotation = ''
  .type=str
  .help="Special types are Pringle-Shen, ..."
  .expert_level=3

#Coordinate systems

convention_override=None
  .type=int
  .help="if defined, override the image-format specific behavior (spot_convention number)"
  .expert_level=3

spot_convention=None
  .type=int
  .help="must be set by the calling program; no default"
  .expert_level=3

#Spotfinder

override_pickled_spotfinders = True
  .type=bool
  .help="automatically erase any existing DISTL_pickle file"
  .expert_level=3

spotfinder_header_tests = True
  .type=bool
  .expert_level=3

spotfinder_mode = 'distl'
  .type=str
  .expert_level=3

spotfinder_verbose = False
  .type=bool
  .expert_level=2

force_method2_resolution_limit = None
  .type=float
  .help="override resolution analysis based on spot count falloff; force spots at least this far out."
  .expert_level=2

distl_lowres_limit = 50.0
  .type=float
  .help="don't pick spots inside this resolution limit"
  .expert_level=3

distl_highres_limit = None
  .type=float
  .help="don't pick spots outside this resolution limit"
  .expert_level=3

distl_binned_image_spot_size = 4
  .type=int
  .expert_level=2

distl_maximum_number_spots_for_indexing = 300
  .type=int
  .expert_level=3

distl_minimum_number_spots_for_indexing = 40
  .type=int
  .expert_level=3

distl_profile_bumpiness = 2
  .type=int
  .help="maximum number of local maxima in good Bragg spots"
  .expert_level=2

distl_report_overloads = True
  .type=bool
  .expert_level=3

distl_keep_Zdata = True
  .type=bool
  .expert_level=3

percent_overlap_forcing_detail = 30.
  .type=float
  .help="detail examination of spots with nearest neighbor analysis and overlap likelihood"
  .expert_level=3

overlapping_spot_criterion = 1.2
  .type=float
  .help="in multiples of the semimajor axis"
  .expert_level=3

spots_pickle = './DISTL_pickle'
  .type=str
  .expert_level=3

distl_spotcenter_algorithm = 'center_of_mass'
  .type=str
  .help="either center_of_mass or maximum_pixel"
  .expert_level=3

distl_permit_binning=True
  .type=bool
  .multiple=False
  .help="Permit binning for large images; set False for Web-Ice since diffimage always renders unbinned."
  .expert_level=2

distl_force_binning=False
  .type=bool
  .multiple=False
  .help="Force binning for all images; only used for development and troubleshooting."
  .expert_level=2

#Data Parameters to Autoindex; Some Affect Spotfinder Also

autoindex_override_beam = None
  .type=floats(size=2)
  .help="x and y coordinates of the direct beam in mm"
  .expert_level=1

autoindex_override_distance = None
  .type=float
  .help="crystal-to-detector distance in mm"
  .expert_level=1

autoindex_override_wavelength = None
  .type=float
  .help="incident wavelength in Angstroms"
  .expert_level=1

autoindex_override_twotheta = None
  .type=float
  .help="detector swing angle in degrees"
  .expert_level=1

autoindex_override_deltaphi = None
  .type=float
  .help="single-shot rotation angle in degrees"
  .expert_level=1

image_specific_osc_start = None
  .type=str
  .help="A lambda x expression giving the rotation in degrees given the image number,
         such as lambda x: x-1.0"
  .expert_level=1

codecamp {
  maxcell = None
    .type=float
    .multiple=False
    .help="Directly specify max unit cell; potentially allow contiguous images"
    .expert_level=2
  minimum_spot_count = None
    .type=int
    .help="For determining spot masks on single images, minimum allowable spot count"
    .expert_level=2
}

pdf_output {

  file=""
    .type=str
    .multiple=False
    .help="If given, specify a file path to output a picture of the sublattice model."
    .expert_level=4
  box_selection="all"
    .type=str
    .multiple=False
    .help="index: show original superlattice | coset: show spots unique to the sublattice | all: default, show both"
    .expert_level=4
  enable_legend=False
    .type=bool
    .multiple=False
    .help="Print the Miller indices, in the triclinic sublattice basis system"
    .expert_level=4
  enable_legend_font_size=10
    .type=float
    .multiple=False
    .help="Print the Miller indices, font size in points"
    .expert_level=4
  enable_legend_ink_color=black
    .type=str
    .multiple=False
    .help="Print the Miller indices, ink color"
    .expert_level=4
  enable_legend_vertical_offset=10
    .type=float
    .multiple=False
    .help="Print the Miller indices, vertical legend offset"
    .expert_level=4
  box_linewidth=0.04
    .type=float
    .multiple=False
    .help="Line width for the rectangular box enclosing the spot profile"
    .expert_level=4
  window_fraction=0.666666
    .type=float
    .multiple=False
    .help="Fractional length of image x,y dimensions rendered to pdf; use fraction for x,y"
    .expert_level=4
  window_offset_x=0.16667
    .type=float
    .multiple=False
    .help="Fractional offset of image x dimension for the window rendered to pdf"
    .expert_level=4
  window_offset_y=0.16667
    .type=float
    .multiple=False
    .help="Fractional offset of image y dimension for the window rendered to pdf"
    .expert_level=4
  markup_inliers=True
    .type=bool
    .multiple=False
    .help="Markup the filtered Bragg candidates, peak and profile center"
    .expert_level=4
  render_all=False
    .type=bool
    .multiple=False
    .help="Show spot predictions for all possible sublattices on separate pages"
    .expert_level=4
  profile_shrink=0
    .type=int
    .multiple=False
    .help="For clarity of view, shrink the profile box by # of pixels"
    .expert_level=4

}

"""
