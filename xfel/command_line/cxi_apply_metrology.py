from __future__ import division
from __future__ import print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cxi.apply_metrology
# $Id
#

import sys, os, pycbf
import libtbx.phil
from libtbx.utils import Usage, Sorry
import tempfile, shutil

master_phil = libtbx.phil.parse("""
source_cbf = None
  .type = str
  .help = cbf file to apply to destination file(s), can be just a header.
  .optional = False
dest_cbf = None
  .type = str
  .help = cbf files on which to apply the metrology measurements from the source.
  .help = Can provide multiple files as once
  .multiple = True
apply_detector_distance = False
  .type = bool
  .help = If True, copy detector distance
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) :
      user_phil.append(libtbx.phil.parse("""dest_cbf=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.source_cbf is None) or not os.path.isfile(params.source_cbf):
    master_phil.show()
    raise Usage("source_cbf must be a file")
  if (params.source_cbf is None):
    master_phil.show()
    raise Usage("dest_cbf must be a file (either dest_cbf=XXX, or the file path(s) alone).")

  print("Source file:", params.source_cbf)
  print("Destination file(s):", end=' ')
  for path in params.dest_cbf:
    print(path, end=' ')
  print()

  # categories required to match between the two files.  Tuples of category names
  # and table column names which are keys, I.E., there should be only one row in the
  # category with a given value in the table
  required_categories = [("diffrn"                      , "id"),
                         ("diffrn_source"               , "diffrn_id"),
                         ("diffrn_detector"             , "id"),
                         ("diffrn_detector_axis"        , "axis_id"),
                         ("diffrn_detector_element"     , "id"),
                         ("diffrn_data_frame"           , "detector_element_id"),
                         ("array_structure_list"        , None),
                         ("array_structure_list_axis"   , "axis_set_id"),
                         ("array_structure_list_section", None)]
  optional_categories = [("diffrn_radiation"            , "diffrn_id"),
                         ("diffrn_radiation_wavelength" , "id"),
                         ("diffrn_measurement"          , "diffrn_id"),
                         ("diffrn_scan"                 , "id"),
                         ("diffrn_scan_frame"           , "frame_id"),
                         ("array_intensities"           , "array_id"),
                         ("array_structure"             , "id"),
                         ("array_data"                  , "array_id")]

  # categories whose data to copy from one cbf to another
  copy_categories =     [("axis"                     , "id"),
                         ("diffrn_scan_axis"         , "axis_id"),
                         ("diffrn_scan_frame_axis"   , "axis_id")]
  req_names = [n[0] for n in required_categories]; req_names.extend([c[0] for c in copy_categories])
  opt_names = [n[0] for n in optional_categories]
  keys      = [n[1] for n in required_categories]; keys .extend([c[1] for c in copy_categories])

  src_cbf = pycbf.cbf_handle_struct()
  src_cbf.read_widefile(params.source_cbf, pycbf.MSG_DIGEST)

  # verify all the categories are present in the source cbf
  print("Testing for required categories in source:")
  src_cbf.select_category(0)
  n_found = 0
  while True:
    if src_cbf.category_name() in req_names:
      print("Found", src_cbf.category_name())
      n_found += 1
    else:
      if src_cbf.category_name() not in opt_names:
        raise Sorry("%s not a recognized category"%src_cbf.category_name())
    try:
      src_cbf.next_category()
    except Exception as e:
      assert "CBF_NOTFOUND" in e.message
      break
  assert n_found == len(req_names)
  print("OK")

  # iterate through the files, validate the required tables and copy the others
  for path in params.dest_cbf:
    print("Validating %s..."%os.path.basename(path), end=' ')

    dst_cbf = pycbf.cbf_handle_struct()
    dst_cbf.read_widefile(path, pycbf.MSG_DIGEST)

    # Validate
    for category, key in required_categories:
      src_cbf.find_category(category)
      dst_cbf.find_category(category)

      assert src_cbf.count_columns() == dst_cbf.count_columns()
      assert src_cbf.count_rows() == dst_cbf.count_rows()

      for j in range(src_cbf.count_rows()):
        src_cbf.rewind_column()
        dst_cbf.rewind_column()

        src_cbf.select_row(j)
        if key is None:
          # for tables with non-unique key column, compare row by row
          dst_cbf.select_row(j)
        else:
          # for tables with unique key column, don't assume the rows are in the
          # same order
          src_cbf.find_column(key)
          dst_cbf.find_column(key)
          dst_cbf.find_row(src_cbf.get_value())

        for i in range(src_cbf.count_columns()):
          src_cbf.select_column(i)
          dst_cbf.find_column(src_cbf.column_name())
          if src_cbf.get_value() != dst_cbf.get_value():
            raise Sorry("Non matching values: table %s, row %d, column %s, %s vs. %s"%(category, i, src_cbf.column_name(),src_cbf.get_value(), dst_cbf.get_value()))

    # Copy
    for category, key in copy_categories:
      src_cbf.find_category(category)
      dst_cbf.find_category(category)

      assert src_cbf.count_columns() == dst_cbf.count_columns()
      assert src_cbf.count_rows() == dst_cbf.count_rows()

      for j in range(src_cbf.count_rows()):
        src_cbf.rewind_column()
        dst_cbf.rewind_column()

        src_cbf.select_row(j)
        src_cbf.find_column(key)

        # don't overwrite detector distance
        if category == "diffrn_scan_frame_axis" and src_cbf.get_value() == "AXIS_D0_Z" and not params.apply_detector_distance:
          continue

        dst_cbf.find_column(key)
        dst_cbf.find_row(src_cbf.get_value())

        for i in range(src_cbf.count_columns()):
          src_cbf.select_column(i)
          dst_cbf.find_column(src_cbf.column_name())

          if src_cbf.get_value() == '.' or src_cbf.get_typeofvalue() == "null":
            dst_cbf.set_typeofvalue("null")
          else:
            dst_cbf.set_value(src_cbf.get_value())

    print("writing cbf...", end=' ')

    t = tempfile.NamedTemporaryFile(delete=False)
    destpath = t.name
    t.close()

    dst_cbf.write_widefile(destpath,pycbf.CBF,\
                           pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)

    del dst_cbf

    shutil.move(destpath, path)

    print("Done")
