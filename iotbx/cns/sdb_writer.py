from __future__ import absolute_import, division, print_function
from iotbx.cns import space_group_symbols
from cctbx import adptbx
from six.moves import cStringIO as StringIO

def write_header(s, file=None, description=None, comment=None,
                    space_group_info=None, n_rows=None):
  assert n_rows is not None
  if (file is not None):
    print("{+ file: %s +}" % str(file), file=s)
  if (description is not None):
    print("{+ description: %s +}" % str(description), file=s)
  if (comment is not None):
    print("{+ comment:", file=s)
    for line in comment: print(line, file=s)
    print("+}", file=s)
  print(file=s)
  print("{- begin block parameter definition -} define(", file=s)
  print(file=s)
  if (space_group_info is not None):
    print("""\
{============================ space group ============================}

{* space group *}
{* use International Table conventions with subscripts substituted
   by parenthesis *}
{===>} sg="%s";
""" % space_group_symbols.cns_format(space_group_info), file=s)
  print("""\
{==================== derivative/MAD coordinates =====================}

{+ list: for each site define:
         - whether the site is to be refined, fixed or ignored
         - derivative name (must be the same for all sites in a derivative)
         - chemical type (note: wavelength-dependent form factors
                                are specified in mad_refine.inp)
         - coordinates (x, y, z)
         - B-value (b)
         - occupancy (q)
         - group name (g) +}

{+ list: the optional group name (g) is a string of upto 4 characters.
         If a group is left blank, each site is refined individually.
         If a group is specified, all sites with the same group name
         and the same derivative name are treated as a rigid body, and their
         occupancies, B-values, and form factors are refined as a group. +}

{+ table: rows=%d numbered
          cols=9 "action" "derivative name" "chemical type"
                 "x coordinate" "y coordinate" "z coordinate"
                 "B-value" "occupancy" "group" +}
""" % n_rows, file=s)

def write_tail(s):
  print("""\
{* to appended new entries or merge this file with other
   site database files use sdb_manipulate.inp *}

{* to delete sites from this file either set the derivative
   name to be blank or use delete_sites.inp *}

{===========================================================================}
{         things below this line do not normally need to be changed         }
{===========================================================================}

) {- end block parameter definition -}""", file=s)

def write_scatterer(s, running_index, scatterer,
                       action=None, segid=None, group=None):
  if (action is None): action = "refine"
  if (segid is None): segid = "SITE"
  if (group is None): group = ""
  assert running_index > 0
  assert scatterer.flags.use_u_iso_only()
  assert action in ("refine", "fix", "ignore")
  i = running_index
  print("""\
{+ choice: "refine" "fix" "ignore" +}
{===>} site.action_%d="%s";
{===>} site.segid_%d="%s"; site.type_%d="%s";
{===>} site.x_%d=%.6g; site.y_%d=%.6g; site.z_%d=%.6g;
{===>} site.b_%d=%.6g; site.q_%d=%.6g; site.g_%d="%s";
""" % (i,action,
       i,segid, i,scatterer.scattering_type,
       i,scatterer.site[0],
       i,scatterer.site[1],
       i,scatterer.site[2],
       i,adptbx.u_as_b(scatterer.u_iso), i,scatterer.occupancy, i,group), file=s)

def xray_structure_as_cns_sdb_file(self, file=None,
                                         description=None,
                                         comment=None,
                                         action=None,
                                         segid=None,
                                         group=None):
  s = StringIO()
  write_header(s, file, description, comment,
                  self.space_group_info(),
                  self.scatterers().size())
  for running_index,scatterer in enumerate(self.scatterers()):
    write_scatterer(s,
      running_index+1,
      scatterer.customized_copy(
        site=self.unit_cell().orthogonalize(scatterer.site)),
      action=action,
      segid=segid,
      group=group)
  write_tail(s)
  return s.getvalue()
