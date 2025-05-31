"""Obsolete method to store and output custom links"""
# MARKED_FOR_DELETION_OLEG
# Reason: sub-optimal approach to store and output custom links.
# Proposal: use origin_id in bond proxies. Write a class capable of
# taking GRM and hierarchy/labels to produce links in various formats.
#
from __future__ import absolute_import, division, print_function
def link_record_output(all_chain_proxies):
  outl = ""
  for key, item in all_chain_proxies.pdb_link_records.items():
    if key=="SSBOND":
      for ssi, (atom1, atom2, sym_op) in enumerate(item):
        #
        def _format_ssbond_atom(atom):
            return "%3s %s %4s%s" % (
              atom.parent().resname,
              atom.parent().parent().parent().id,
              atom.parent().parent().resseq,
              atom.parent().parent().icode,
              )
        #
        if str(sym_op)!="x,y,z":
          continue
        outl += "SSBOND %3s %s   %s" % (
          ssi+1,
          _format_ssbond_atom(atom1),
          _format_ssbond_atom(atom2),
          )
        assert str(sym_op)=="x,y,z"
        outl += "\n"
    elif key=="LINK":
      for ssi, (atom1, atom2, sym_op) in enumerate(item):
        def _format_link_atom(atom):
          altloc = atom.parent().altloc
          if not altloc: altloc=" "
          return "%4s%s%3s%2s%4s%s" % (
            atom.name,
            altloc,
            atom.parent().resname,
            atom.parent().parent().parent().id,
            atom.parent().parent().resseq,
            atom.parent().parent().icode,
            )
        if str(sym_op)!="x,y,z":
          continue
        outl += "LINK%s%s%s%s" % (
          ' '*8,
          _format_link_atom(atom1),
          ' '*15,
          _format_link_atom(atom2),
          )
        outl += "\n"
#          print """         1         2         3         4         5         6#         7
#123456789012345678901234567890123456789012345678901234567890123456789012
#LINK         O1  DDA     1                 C3  DDL     2
#LINK        MN    MN   391                 OE2 GLU   217            2565"""
#          print outl
    else:
      raise Sorry("PDB record type %s unknown" % key)
  return outl[:-1]
# END_MARKED_FOR_DELETION_OLEG
