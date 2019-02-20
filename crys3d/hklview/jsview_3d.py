
# TODO:
#  - cached scenes

from __future__ import division
from libtbx.math_utils import roundoff
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
    self.miller_array = None
    self.d_min = None
    self.scene = None
    self.animation_time = 0
    self.cameratype = "orthographic"

  def set_miller_array (self, miller_array, merge=None, details="") :
    if (miller_array is None) : return
    self.miller_array = miller_array
    self.merge = merge
    self.d_min = miller_array.d_min()
    array_info = miller_array.info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % miller_array.unit_cell().parameters()
    print "Data: %s %s, %d reflections in space group: %s, unit Cell: %s" \
      % (array_info.label_string(), details, miller_array.indices().size(),
          miller_array.space_group_info(), uc)

    self.construct_reciprocal_space(merge=merge)


  def construct_reciprocal_space (self, merge=None) :
    self.scene = display.scene(miller_array=self.miller_array,
      merge=merge,
      settings=self.settings)
    self.rotation_center = (0,0,0)


  def DrawNGLJavaScript(self):
    if self.miller_array is None :
      print "A miller array must be selected for drawing"
      return

    h_axis = self.scene.axes[0]
    k_axis = self.scene.axes[1]
    l_axis = self.scene.axes[2]

    Hstararrowstart = roundoff( [-h_axis[0]*100, -h_axis[1]*100, -h_axis[2]*100] )
    Hstararrowend = roundoff( [h_axis[0]*100, h_axis[1]*100, h_axis[2]*100] )
    Hstararrowtxt  = roundoff( [h_axis[0]*102, h_axis[1]*102, h_axis[2]*102] )
    Kstararrowstart = roundoff( [-k_axis[0]*100, -k_axis[1]*100, -k_axis[2]*100] )
    Kstararrowend = roundoff( [k_axis[0]*100, k_axis[1]*100, k_axis[2]*100] )
    Kstararrowtxt  = roundoff( [k_axis[0]*102, k_axis[1]*102, k_axis[2]*102] )
    Lstararrowstart = roundoff( [-l_axis[0]*100, -l_axis[1]*100, -l_axis[2]*100] )
    Lstararrowend = roundoff( [l_axis[0]*100, l_axis[1]*100, l_axis[2]*100] )
    Lstararrowtxt  = roundoff( [l_axis[0]*102, l_axis[1]*102, l_axis[2]*102] )

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
    data = self.scene.data
    hkls = self.scene.indices
    dres = self.scene.work_array.d_spacings().data()
    colstr = self.scene.miller_array.info().label_string()
    assert (colors.size() == radii.size() == self.scene.points.size())
    shapespherestr = ""
    for i, hklstars in enumerate(points) :
      tooltip = "'H,K,L: %s, %s, %s" %(hkls[i][0], hkls[i][1], hkls[i][2])
      tooltip += "\\ndres: %s" %str(roundoff(dres[i])  )
      tooltip += "\\n%s: %s" %(colstr, str(roundoff(data[i]) ) )
      tooltip += "'"
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      shapespherestr += "shape.addSphere( %s, %s, %s, %s);\n" \
         %(str(roundoff(list(hklstars))), str(roundoff(list(colors[i]), 2)),
            str(roundoff(radii[i], 2)), tooltip )
    documentstr = """
    window.addEventListener( 'resize', function( event ){
        stage.handleResize();
    }, false );

        var stage;

        document.addEventListener('DOMContentLoaded', function () {

          stage = new NGL.Stage('viewport', { backgroundColor: "grey" });
          var shape = new NGL.Shape('shape');
          stage.setParameters( { cameraType: "%s" } );


    %s

    %s

    var shapeComp = stage.addComponentFromObject(shape);
    shapeComp.addRepresentation('buffer');
    shapeComp.autoView();

     });

    """ % (self.cameratype, arrowstr, shapespherestr)
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
      self.GetParent().update_clicked(index=self.closest_point_i_seq)
      #hkl, d_min, value = self.scene.get_reflection_info(
      #  self.closest_point_i_seq)
      #self.GetParent().update_clicked(hkl, d_min, value)
    else :
      self.GetParent().update_clicked(index=None)



