
# TODO:
#  - cached scenes

from __future__ import division
from cctbx.miller import display

class hklview_3d () :
  def __init__ (self, mysettings) :
    # FIXME orthographic is definitely best for this application, but it isn't
    # working properly right now
    #self.orthographic = True
    self.settings = mysettings
    self.buffer_factor = 2.0
    self.min_slab = 4
    self.min_viewport_use_fraction = 0.1
    self.min_dist = 4.0
    self.flag_show_fog = True
    self.flag_use_lights = True
    self.flag_use_quadrics = False
    self.minimum_covering_sphere = None
    self.spheres_display_list = None
    self.points_display_list = None
    self.labels_display_list = None
    self.miller_array = None
    self.d_min = None
    self.scene = None
    self.animation_time = 0
    #self.fps = gltbx.viewer_utils.fps_monitor()
    # XXX prevent exception when no data are loaded
    from scitbx.math import minimum_covering_sphere
    from scitbx.array_family import flex
    points = flex.vec3_double([(0.0,0.0,0.0),(1.0,1.0,1.0)])
    mcs = minimum_covering_sphere(points=points, epsilon=0.1)
    self.minimum_covering_sphere = mcs

  def set_miller_array (self, miller_array, zoom=False, merge=None) :
    if (miller_array is None) : return
    self.miller_array = miller_array
    self.merge = merge
    self.d_min = miller_array.d_min()
    self.construct_reciprocal_space(merge=merge)



  def construct_reciprocal_space (self, merge=None) :
    self.scene = display.scene(miller_array=self.miller_array,
      merge=merge,
      settings=self.settings)
    from scitbx.math import minimum_covering_sphere
    mcs = minimum_covering_sphere(points=self.scene.points,
                                  epsilon=0.1)
    self.minimum_covering_sphere = mcs
    self.spheres_display_list = None
    self.points_display_list = None
    self.labels_display_list = None
    self.rotation_center = (0,0,0)


  def DrawNGLJavaScript(self):
    if self.miller_array is None :
      return

    h_axis = self.scene.axes[0]
    k_axis = self.scene.axes[1]
    l_axis = self.scene.axes[2]

    Hstararrowstart = [-h_axis[0]*100, -h_axis[1]*100, -h_axis[2]*100]
    Hstararrowend = [h_axis[0]*100, h_axis[1]*100, h_axis[2]*100]
    Hstararrowtxt  = [h_axis[0]*102, h_axis[1]*102, h_axis[2]*102]
    Kstararrowstart = [-k_axis[0]*100, -k_axis[1]*100, -k_axis[2]*100]
    Kstararrowend = [k_axis[0]*100, k_axis[1]*100, k_axis[2]*100]
    Kstararrowtxt  = [k_axis[0]*102, k_axis[1]*102, k_axis[2]*102]
    Lstararrowstart = [-l_axis[0]*100, -l_axis[1]*100, -l_axis[2]*100]
    Lstararrowend = [l_axis[0]*100, l_axis[1]*100, l_axis[2]*100]
    Lstararrowtxt  = [l_axis[0]*102, l_axis[1]*102, l_axis[2]*102]

    arrowstr = """
    // xyz arrows
    var a=[0,0,0];
    shape.addSphere( a , [ 1, 1, 1 ], 0.3, 'Origo');
    //blue-x
    shape.addArrow( %s, %s , [ 0, 0, 1 ], 0.1);
    //green-y
    shape.addArrow( %s, %s , [ 0, 1, 0 ], 0.1);
    //red-z
    shape.addArrow( %s, %s , [ 1, 0, 0 ], 0.1);

    shape.addText( %s, [ 0, 0, 1 ], 5, 'H');
    shape.addText( %s, [ 0, 1, 0 ], 5, 'K');
    shape.addText( %s, [ 1, 0, 0 ], 5, 'L');

    """ %(str(Hstararrowstart), str(Hstararrowend), str(Kstararrowstart), str(Kstararrowend),
          str(Lstararrowstart), str(Lstararrowend), Hstararrowtxt, Kstararrowtxt, Lstararrowtxt)

    colors = self.scene.colors
    radii = self.scene.radii * self.settings.scale
    points = self.scene.points
    assert (colors.size() == radii.size() == self.scene.points.size())
    shapespherestr = ""
    for i, hkl in enumerate(points) :
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      shapespherestr += "shape.addSphere( %s , %s, %s, %s);\n" \
         %(str(list(hkl)), str(list(colors[i])), str(radii[i]), "' '" )
    documentstr = """
    window.addEventListener( 'resize', function( event ){
        stage.handleResize();
    }, false );

        var stage;

        document.addEventListener('DOMContentLoaded', function () {

          stage = new NGL.Stage('viewport', { backgroundColor: "grey" });
          var shape = new NGL.Shape('shape');

    %s

    %s

    var shapeComp = stage.addComponentFromObject(shape);
    shapeComp.addRepresentation('buffer');
    shapeComp.autoView();

     });

    """ % (arrowstr, shapespherestr)
    with open( r"C:\Users\oeffner\Buser\NGL_HKLviewer\myjstr3.js", "w") as f:
      f.write(documentstr)

  #--- user input and settings
  def update_settings (self) :
    self.construct_reciprocal_space(merge=self.merge)
    self.DrawNGLJavaScript()

  def process_pick_points (self) :
    self.closest_point_i_seq = None
    if (self.pick_points is not None) and (self.scene is not None) :
      closest_point_i_seq = gltbx.viewer_utils.closest_visible_point(
        points=self.scene.points,
        atoms_visible=self.scene.visible_points,
        point0=self.pick_points[0],
        point1=self.pick_points[1])
      if (closest_point_i_seq is not None) :
        self.closest_point_i_seq = closest_point_i_seq
    if (self.closest_point_i_seq is not None) :
      self.scene.label_points.add(self.closest_point_i_seq)
      self.labels_display_list = None
      self.GetParent().update_clicked(index=self.closest_point_i_seq)
      #hkl, d_min, value = self.scene.get_reflection_info(
      #  self.closest_point_i_seq)
      #self.GetParent().update_clicked(hkl, d_min, value)
    else :
      self.GetParent().update_clicked(index=None)



