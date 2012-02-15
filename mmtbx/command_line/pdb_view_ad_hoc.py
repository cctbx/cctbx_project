# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from gltbx import wx_viewer
from libtbx.option_parser import libtbx_option_parser
import wx
import sys

class viewer(wx_viewer.show_points_and_lines_mixin):

  def set_points_and_lines(
        O,
        app,
        minimum_covering_sphere_view_scale=1.3):
    pdb_atoms = app.pdb_atoms
    pdb_atoms.set_chemical_element_simple_if_necessary()
    atom_tmp_sentinel = pdb_atoms.reset_tmp(first_value=0, increment=0)
    for i,rg in enumerate(app.pdb_hierarchy.residue_groups()):
      for atom in rg.atoms():
        atom.tmp = i
    O.points = pdb_atoms.extract_xyz()
    if (app.co.serial_labels):
      m = ""
      need_m = (app.pdb_hierarchy.models_size() != 1)
      for mdl in app.pdb_hierarchy.models():
        if (need_m): m = mdl.id.strip() + ":"
        for i in xrange(mdl.atoms_size()):
          O.labels.append(m+str(i))
      assert len(O.labels) == O.points.size()
    else:
      rg_done = set()
      for atom in pdb_atoms:
        i = atom.tmp
        if (i not in rg_done):
          rg_done.add(i)
          l = atom.id_str()
        else:
          l = atom.name
        O.labels.append(l)
    from cctbx.crystal.distance_based_connectivity import \
      build_simple_two_way_bond_sets
    bond_sets = build_simple_two_way_bond_sets(
      sites_cart=O.points,
      elements=pdb_atoms.extract_element())
    for i,bond_set in enumerate(bond_sets):
      for j in bond_set:
        if (i < j):
          line = (i,j)
          ai, aj = [pdb_atoms[_] for _ in line]
          if (ai.is_in_same_conformer_as(aj)):
            O.line_i_seqs.append(line)
            if (ai.tmp == aj.tmp):
              if (ai.tmp % 2 == 0):
                color = (0,0,1)
              else:
                color = (0,1,0)
            else:
              color = (1,0,0)
            O.line_colors[line] = color
    del atom_tmp_sentinel
    from scitbx.math import minimum_covering_sphere, sphere_3d
    mcs = minimum_covering_sphere(O.points, epsilon=1.e-2)
    O.minimum_covering_sphere = sphere_3d(
      center=mcs.center(),
      radius=mcs.radius()*minimum_covering_sphere_view_scale)
    O.flag_show_minimum_covering_sphere = False
    O.flag_show_rotation_center = False
    _ = app.co.labels_threshold
    O.flag_show_labels = (_ == 0 or len(O.points) <= _)
    O.labels_display_list = None
    O.lines_display_list = None
    O.points_display_list = None

class App(wx_viewer.App):

  def __init__(O, args):
    import libtbx.load_env
    command_line = (libtbx_option_parser(
      usage="%s [options] pdb_file" % libtbx.env.dispatcher_name)
      .option(None, "--labels_threshold",
        action="store",
        type="int",
        default=20,
        help="do not show atom labels if more than given number of atoms",
        metavar="INT")
      .option(None, "--serial_labels",
        action="store_true",
        default=False)
    ).process(args=args, nargs=1)
    O.co = command_line.options
    file_name = command_line.args[0]
    import iotbx.pdb
    O.pdb_inp = iotbx.pdb.input(file_name=file_name)
    O.pdb_hierarchy = O.pdb_inp.construct_hierarchy()
    O.pdb_atoms = O.pdb_hierarchy.atoms()
    print file_name
    print "  number of models:", O.pdb_hierarchy.models_size()
    print "  number of atoms:", O.pdb_atoms.size()
    super(App, O).__init__(title=libtbx.env.dispatcher_name)

  def init_view_objects(O):
    box = wx.BoxSizer(wx.VERTICAL)
    O.view_objects = viewer(O.frame, size=(600,600))
    O.view_objects.set_points_and_lines(O)
    box.Add(O.view_objects, wx.EXPAND, wx.EXPAND)
    O.frame.SetSizer(box)
    box.SetSizeHints(O.frame)

def run(args):
  App(args).MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
