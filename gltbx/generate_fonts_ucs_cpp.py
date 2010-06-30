from __future__ import division
import libtbx.load_env
from libtbx.str_utils import show_string
import sys, os
op = os.path

class font_info(object):

  def __init__(O, short_name, file_name, width, height, xorig, yorig):
    O.short_name = short_name
    O.file_name = file_name
    O.width = width
    O.height = height
    O.xorig = xorig
    O.yorig = yorig

font_infos = [
  font_info("8x13", "8x13.bdf", 8, 13, 0, -2),
  font_info("9x15", "9x15.bdf", 9, 15, 0, -3),
  font_info("10x20", "10x20.bdf", 10, 20, 0, -4)
]

class read_bitmap(object):

  __slots__ = ["label", "encoding", "swidth", "dwidth", "bbx", "bitmap"]

  def __init__(O, bdf_file):
    """Converts a block of this form:
STARTCHAR char0
ENCODING 0
SWIDTH 568 0
DWIDTH 8 0
BBX 8 13 0 -2
BITMAP
00
...
00
ENDCHAR
"""
    O.label = None
    O.encoding = None
    O.swidth = None
    O.dwidth = None
    O.bbx = None
    O.bitmap = None
    line = bdf_file.next().strip()
    if (line == "ENDFONT"):
      return
    assert line.startswith("STARTCHAR ")
    O.label = line.split(None, 1)[1]
    while True:
      line = bdf_file.next().strip()
      if (line.startswith("ENCODING ")):
        fields = line.split()
        assert len(fields) == 2
        O.encoding = int(fields[1])
      elif (line.startswith("SWIDTH ")):
        fields = line.split()
        assert len(fields) == 3
        O.swidth = (int(fields[1]), int(fields[2]))
      elif (line.startswith("DWIDTH ")):
        fields = line.split()
        assert len(fields) == 3
        O.dwidth = (int(fields[1]), int(fields[2]))
      elif (line.startswith("BBX ")):
        fields = line.split()
        assert len(fields) == 5
        O.bbx = [int(field) for field in fields[1:]]
      elif (line == "BITMAP"):
        break
      else:
        raise RuntimeError(
          "Font file %s: unknown line in STARTCHAR block: %s" % (
            show_string(bdf_file.name), line))
    assert O.encoding is not None
    O.bitmap = []
    while True:
      line = bdf_file.next().strip()
      if (line == "ENDCHAR"):
        break
      O.bitmap.append(line)
    assert len(O.bitmap) == O.bbx[1]

  def as_glbitmap(O):
    result = []
    w,h = O.bbx[:2]
    n_bytes = w//8
    mask = (256**n_bytes)-1
    remainder = w - n_bytes*8
    if (remainder > 0):
      n_bytes += 1
      mask <<= remainder
      mask |= (2**remainder)-1
      mask <<= 8-remainder
    n_hex = 2*n_bytes
    padding = "0"*n_hex
    rows = list(O.bitmap)
    rows.reverse()
    for row in rows:
      v = int((row+padding)[:n_hex],16) & mask
      bytes = []
      for i in xrange(n_bytes):
        bytes.append(v & 255)
        v >>= 8
      bytes.reverse()
      result.extend(bytes)
    return result

  def format_cpp(O):
    return "/* %5d */ " % O.encoding \
         + ",".join(["0x%.2x" % i for i in O.as_glbitmap()])

class encoding_range(object):

  __slots__ = ["start", "count"]

  def __init__(O, start):
    O.start = start
    O.count = 1

class encoding_ranges(object):

  __slots__ = ["ranges"]

  def __init__(O, bitmaps):
    O.ranges = []
    current_range = None
    previous_encoding = -2
    for bitmap in bitmaps:
      if (previous_encoding+1 == bitmap.encoding):
        current_range.count += 1
      else:
        current_range = encoding_range(start=bitmap.encoding)
        O.ranges.append(current_range)
      previous_encoding = bitmap.encoding

  def format_cpp(O):
    return ",\n".join([
       "%d, %d" % (range.start, range.count)
         for range in O.ranges])

def format_font_ucs_cpp(
      short_name,
      full_name,
      width,
      height,
      xorig,
      yorig,
      char_records,
      encoding_ranges):
  number_of_chars = len(char_records)
  raw_bitmaps = ",\n".join(char_records)
  return """\
#include <gltbx/fonts_ucs.h>

namespace gltbx { namespace fonts { namespace ucs {

namespace {

static const unsigned char raw_bitmaps[] = {
%(raw_bitmaps)s
};

static const encoding_range encoding_ranges[] = {
%(encoding_ranges)s,
0, 0
};

} // namespace <anonymous>

bitmap_font_record bitmap_%(width)dx%(height)d = {
  "%(short_name)s",
  "%(full_name)s",
  %(width)d,
  %(height)d,
  %(xorig)d,
  %(yorig)d,
  %(number_of_chars)d,
  sizeof(raw_bitmaps),
  raw_bitmaps,
  encoding_ranges,
};

}}} // namespace gltbx::fonts::ucs
""" % vars()

def convert(ucs_fonts_dir, target_dir, font_info):
  bdf_file = open(op.join(ucs_fonts_dir, font_info.file_name))
  full_name = None
  number_of_chars = None
  for line in bdf_file:
    if (line.startswith("FONT ")):
      assert full_name is None
      fields = line.split()
      assert len(fields) == 2
      full_name=fields[1]
    if (line.startswith("CHARS ")):
      assert number_of_chars is None
      fields = line.split()
      assert len(fields) == 2
      number_of_chars = int(fields[1])
      break
  assert full_name is not None
  assert number_of_chars is not None
  bitmaps = []
  while True:
    bitmap = read_bitmap(bdf_file=bdf_file)
    if (bitmap.label is None):
      break
    bitmaps.append(bitmap)
  assert len(bitmaps) == number_of_chars
  char_records = [bitmap.format_cpp() for bitmap in bitmaps]
  cpp_file_name = op.join(target_dir, "font_ucs_%s.cpp" % font_info.short_name)
  open(cpp_file_name, "w").write(
    format_font_ucs_cpp(
      short_name=font_info.short_name,
      full_name=full_name,
      width=font_info.width,
      height=font_info.height,
      xorig=font_info.xorig,
      yorig=font_info.yorig,
      char_records=char_records,
      encoding_ranges=encoding_ranges(bitmaps=bitmaps).format_cpp()))

def run(target_dir):
  if (not op.isdir(target_dir)):
    os.makedirs(target_dir)
  for relative_path in ["gui_resources/ucs-fonts", "ucs-fonts"]:
    ucs_fonts_dir = libtbx.env.find_in_repositories(
      relative_path=relative_path, test=op.isdir, optional=True)
    if (ucs_fonts_dir is not None):
      break
  else:
    raise RuntimeError("Cannot find ucs-fonts directory.")
  done_flag_file = op.join(target_dir, "FONTS_UCS_DONE_FLAG_FILE")
  if (op.isfile(done_flag_file)):
    print "      Info: Re-using existing font cpp files."
    print "      Hint: Remove %s" % op.join(
      op.basename(target_dir), op.basename(done_flag_file))
    print "            to force generation of new font files."
    return
  print "      fonts:",
  for font_info in font_infos:
    print font_info.short_name,
    sys.stdout.flush()
    convert(
      ucs_fonts_dir=ucs_fonts_dir, target_dir=target_dir, font_info=font_info)
  open(done_flag_file, "w")
  print

if (__name__ == "__main__"):
  run(sys.argv[1])
