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
