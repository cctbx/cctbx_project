from iotbx.cns import space_group_symbols
from cctbx import adptbx
from cStringIO import StringIO

def write_header(s, file=None, description=None, comment=None,
                    space_group_info=None, n_rows=None):
  assert n_rows != None
  if (file != None):
    print >> s, "{+ file: %s +}" % str(file)
  if (description != None):
    print >> s, "{+ description: %s +}" % str(description)
  if (comment != None):
    print >> s, "{+ comment:"
    for line in comment: print >> s, line
    print >> s, "+}"
  print >> s
  print >> s, "{- begin block parameter definition -} define("
  print >> s
  if (space_group_info != None):
    print >> s, """\
{============================ space group ============================}

{* space group *}
{* use International Table conventions with subscripts substituted
   by parenthesis *}
{===>} sg="%s";
""" % space_group_symbols.cns_format(space_group_info)
  print >> s, """\
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
""" % n_rows

def write_tail(s):
  print >> s, """\
{* to appended new entries or merge this file with other
   site database files use sdb_manipulate.inp *}

{* to delete sites from this file either set the derivative
   name to be blank or use delete_sites.inp *}

{===========================================================================}
{         things below this line do not normally need to be changed         }
{===========================================================================}

) {- end block parameter definition -}"""

def write_scatterer(s, running_index, scatterer,
                       action="fix", segid="SITE", group=""):
  assert running_index > 0
  assert not scatterer.anisotropic_flag
  assert action in ("refine", "fix", "ignore")
  i = running_index
  print >> s, """\
{+ choice: "refine" "fix" "ignore" +}
{===>} site.action_%d="%s";
{===>} site.segid_%d="%s"; site.type_%d="%s";
{===>} site.x_%d=%.6g; site.y_%d=%.6g; site.z_%d=%.6g;
{===>} site.b_%d=%.6g; site.q_%d=%.6g; site.g_%d="%s";
""" % (i,action,
       i,segid, i,scatterer.caasf.label(),
       i,scatterer.site[0],
       i,scatterer.site[1],
       i,scatterer.site[2],
       i,adptbx.u_as_b(scatterer.u_iso), i,scatterer.occupancy, i,group)

def xray_structure_as_cns_sdb_file(self, file=None,
                                         description=None,
                                         comment=None,
                                         group=""):
  s = StringIO()
  write_header(s, file, description, comment,
                  self.space_group_info(),
                  self.scatterers().size())
  for running_index,scatterer in self.scatterers().items():
    write_scatterer(s,
      running_index+1,
      scatterer.copy(site=self.unit_cell().orthogonalize(scatterer.site)),
      group=group)
  write_tail(s)
  return s.getvalue()
