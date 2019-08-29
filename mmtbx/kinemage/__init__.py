from __future__ import absolute_import, division, print_function

def kin_vec(start_key, start_xyz, end_key, end_xyz, width=None):
  start_altloc = start_key[0:1]
  if start_altloc == ' ':
    start_altloc_txt = ""
  else:
    start_altloc_txt = " '%s'" % start_altloc.lower()
  end_altloc = end_key[0:1]
  if end_altloc == ' ':
    end_altloc_txt = ""
  else:
    end_altloc_txt = " '%s'" % end_altloc.lower()
  if width is None:
    return "{%s} P%s %.3f %.3f %.3f {%s} L%s %.3f %.3f %.3f\n" % (
           start_key,
           start_altloc_txt,
           start_xyz[0],
           start_xyz[1],
           start_xyz[2],
           end_key,
           end_altloc_txt,
           end_xyz[0],
           end_xyz[1],
           end_xyz[2])
  else:
    return "{%s} P%s %.3f %.3f %.3f {%s} L%s width%d %.3f %.3f %.3f\n" % (
           start_key,
           start_altloc_txt,
           start_xyz[0],
           start_xyz[1],
           start_xyz[2],
           end_key,
           end_altloc_txt,
           width,
           end_xyz[0],
           end_xyz[1],
           end_xyz[2])
