from PyQt4 import QtGui
from PyQt4.QtCore import Qt
from crys3d import qttbx
import crys3d.qttbx.xray_structure_viewer_controls
from gltbx import quadrics, gl_managed
from gltbx.gl import *
from gltbx.glu import *
from cctbx import crystal
from cctbx.array_family import flex
from cctbx.eltbx import covalent_radii
from scitbx import matrix as mat
import itertools


def display(**kwds):
  app = QtGui.QApplication([])
  view = xray_structure_viewer(**kwds)
  view.resize(800, 800)
  ctrls = xray_structure_viewer_controls(view)
  ctrls.show()
  view.show()
  app.exec_()


class xray_structure_viewer_controls(
  qttbx.widget_control_mixin,
  qttbx.xray_structure_viewer_controls.Ui_Form):

  def __init__(self, xray_structure_viewer):
    qttbx.widget_control_mixin.__init__(self, xray_structure_viewer,
                                        None, Qt.Tool)
    self.move(0,50)
    self.view.set_labels_type(self.labels.currentText())
    self.labels.currentIndexChanged[str].connect(self.view.set_labels_type)
    self.labelFontSize.setValue(self.view.label_font.pointSize())
    self.labelFontSize.valueChanged[int].connect(
      self.view.set_font_size)


class xray_structure_viewer(qttbx.widget):

  material_for = ([
    ('Br', (0.97, 0.86, 0.03), (1   , 0.5 , 0   ) ),
    ( 'C', (0.32,)*3         , (0.75,)*3          ),
    ( 'N', (0   , 0   , 1   ), (0.37, 0.37, 0.63) ),
    ( 'O', (0.91, 0   , 0   ), (0.63, 0.37, 0.37) ),
    ( 'F', (0   , 1   , 0   ), (0.07, 0.35, 0.07) ),
    ('Al', (0   , 0.5 , 0.5 ), (0.37, 0.87, 0.87) ),
    ('Si', (0.98, 0.42, 0.01), (0.62, 0.62, 0.37) ),
    ( 'P', (0.5 , 0   , 0.5 ), (0.5,)*3           ),
    ( 'S', (0.97, 0.85, 0.03), (0.5,)*3           ),
    ('Cl', (0   , 0.25, 0   ), (0.13, 0.65, 0.13) ),
    ('Br', (0.51, 0   , 0   ), (0.72, 0.53, 0.53) ),
    ( 'I', (0.27, 0   , 0.27), (0.87, 0.37, 0.87) ),
    ]

    # Metals (1st row)
    + [
    ( elt, (0   , 0   , 0.49), (0.37, 0.37, 0.62) )
    for elt in ('Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')
    ]
    )

  default_material = gl_managed.material_model(
    ambient_front_colour=(0.2,)*3,
    diffuse_front_colour=(0.1,)*3)

  material_for = dict([
    (elt, gl_managed.material_model(ambient_front_colour=a,
                                    diffuse_front_colour=b))
    for (elt, a, b) in material_for
    ])

  covalent_bond_tolerance = 0.5 # Angstrom
  distance_cutoff = None
  bonding="covalent"

  def __init__(self, xray_structure, name='??',
               **kwds):
    super(xray_structure_viewer, self).__init__(
      unit_cell=xray_structure.unit_cell(),
      orthographic=True,
      light_position=(-1, 1, 1, 0),
      **kwds)
    assert self.bonding in ("covalent", "all")
    assert self.bonding != "all" or self.distance_cutoff is not None
    self.xray_structure = xs = xray_structure
    self.setWindowTitle("%s in %s" % (name,
                                      xs.space_group().type().hall_symbol()))
    sites_frac = xs.sites_frac()
    self.set_extent(sites_frac.min(), sites_frac.max())
    self.is_unit_cell_shown = False

    sites_cart = self.sites_cart = xs.sites_cart()
    thermal_tensors = xs.extract_u_cart_plus_u_iso()
    self.ellipsoid_to_sphere_transforms = {}
    self.scatterer_indices = flex.std_string()
    self.scatterer_labels = flex.std_string()
    for i, (sc, site, u_cart) in enumerate(itertools.izip(xs.scatterers(),
                                                          sites_cart,
                                                          thermal_tensors)):
      t = quadrics.ellipsoid_to_sphere_transform(site, u_cart)
      self.ellipsoid_to_sphere_transforms.setdefault(
        sc.element_symbol(),
        quadrics.shared_ellipsoid_to_sphere_transforms()).append(t)
      self.scatterer_indices.append("# %i" % i)
      self.scatterer_labels.append(sc.label)
    self.labels = None
    self.label_font = QtGui.QFont("Arial Black", pointSize=18)

    if self.bonding == "covalent":
      radii = [
        covalent_radii.table(elt).radius()
        for elt in xs.scattering_type_registry().type_index_pairs_as_dict() ]
      buffer_thickness = 2*max(radii) + self.covalent_bond_tolerance
      asu_mappings = xs.asu_mappings(buffer_thickness=buffer_thickness)
      bond_table = crystal.pair_asu_table(asu_mappings)
      bond_table.add_covalent_pairs(xs.scattering_types(),
                                    tolerance=self.covalent_bond_tolerance)
    elif self.bonding == "all":
      asu_mappings = xs.asu_mappings(buffer_thickness=self.distance_cutoff)
      bond_table = crystal.pair_asu_table(asu_mappings)
      bond_table.add_all_pairs(self.distance_cutoff)

    pair_sym_table = bond_table.extract_pair_sym_table(
      all_interactions_from_inside_asu=True)
    self.bonds = flex.vec3_double()
    self.bonds.reserve(len(xs.scatterers()))
    uc = self.xray_structure.unit_cell()
    frac = mat.rec(uc.fractionalization_matrix(), (3,3))
    inv_frac = frac.inverse()
    site_symms = xs.site_symmetry_table()
    scatt = self.xray_structure.scatterers()
    for i, neighbours in enumerate(pair_sym_table):
      x0 = sites_cart[i]
      sc0 = scatt[i]
      for j, ops in neighbours.items():
        sc1 = scatt[j]
        if sc0.scattering_type == 'H' and sc1.scattering_type == 'H':
          continue
        for op in ops:
          if op.is_unit_mx():
            x1 = sites_cart[j]
          else:
            x1 = uc.orthogonalize(op*sites_frac[j])
            op_cart = inv_frac*mat.rec(op.r().as_double(), (3,3))*frac
            u1 = (op_cart
                  *mat.sym(sym_mat3=thermal_tensors[j])
                  *op_cart.transpose())
            t = quadrics.ellipsoid_to_sphere_transform(x1, u1.as_sym_mat3())
            self.ellipsoid_to_sphere_transforms[sc1.element_symbol()].append(t)
            self.sites_cart.append(x1)
            op_lbl = (" [%s]" % op).lower()
            self.scatterer_indices.append("# %i%s" % (j, op_lbl))
            self.scatterer_labels.append("%s%s" % (sc1.label, op_lbl))
          self.bonds.append(x0)
          self.bonds.append(x1)

  def initialise_opengl(self):
    self.ellipsoid_proto = quadrics.proto_ellipsoid(
      slices=32, stacks=32)
    self.principal_ellipses_tex = \
        quadrics.ellipsoid_principal_sections_texture(darkening=0.75,
                                                      n_s=64, n_t=64)
    self.cylindre_proto = quadrics.proto_cylinder(slices=16)
    self.bond_material = gl_managed.material_model(
      ambient_front_colour=(0.2,)*3,
      diffuse_front_colour=(0.1,)*3)
    glEnable(GL_TEXTURE_2D)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)

  def draw_object_in_cartesian_coordinates(self):
    self.principal_ellipses_tex.bind()
    for element, transforms in self.ellipsoid_to_sphere_transforms.iteritems():
      material = self.material_for.get(element, self.default_material)
      material.execute()
      transforms.draw(self.ellipsoid_proto)
    self.principal_ellipses_tex.unbind()

    self.bond_material.execute()
    for i in xrange(0, len(self.bonds), 2):
      start, end = self.bonds[i], self.bonds[i+1]
      self.cylindre_proto.draw(start, end, base_radius=0.05)

    if self.labels is not None:
      glPushAttrib(GL_LIGHTING_BIT)
      glDisable(GL_LIGHTING)
      glPushAttrib(GL_DEPTH_BUFFER_BIT)
      glDisable(GL_DEPTH_TEST)
      glColor3f(1, 1, 1)
      e = 0.1
      for x, lbl in itertools.izip(self.sites_cart, self.labels):
        self.renderText(x[0]-e, x[1]+e, x[2]-e,
                        lbl,
                        self.label_font)
      glPopAttrib()
      glPopAttrib()

  def set_labels_type(self, kind):
    self.labels = getattr(self, str(kind).lower().replace(' ', '_'), None)
    self.updateGL()

  def set_font_size(self, s):
    self.label_font.setPointSize(s)
    self.updateGL()
