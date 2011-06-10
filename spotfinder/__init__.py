phil_str = """\
spotfinder {
  scanbox_windows = 101 51 51
    .type = ints(size=3, value_min=10)
    .help = "Integer scanbox sizes for calculating background,"
            " for cycles 1,2, and 3."
  peripheral_margin = 20
    .type = int(value_min=0)
    .help = "No spot detection inside margin; width in pixels"
}
"""
