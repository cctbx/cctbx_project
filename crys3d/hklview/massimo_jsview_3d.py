
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
    self.jscriptfname = ""
    self.NGLscriptstr = ""
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
    nrefls = self.scene.points.size()

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
    assert (colors.size() == radii.size() == nrefls)
    shapespherestr = ""
    shapespherestr2 = ""
    shapespherestr3 = ""
    shapespherestr4 = ""
    shapespherestr5 = ""
    colours = []
    positions = []
    radii2 = []
    tooltips =[]
    spbufttips = []
    for i, hklstars in enumerate(points) :
      tooltip = "'H,K,L: %s, %s, %s" %(hkls[i][0], hkls[i][1], hkls[i][2])
      tooltip += "\\ndres: %s" %str(roundoff(dres[i])  )
      tooltip += "\\n%s: %s" %(colstr, str(roundoff(data[i]) ) )
      tooltip += "'"
      spbufttip = "H,K,L: %s, %s, %s" %(hkls[i][0], hkls[i][1], hkls[i][2])
      spbufttip += "\ndres: %s" %str(roundoff(dres[i])  )
      spbufttip += "\n%s: %s" %(colstr, str(roundoff(data[i]) ) )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )

      shapespherestr += "shape.addSphere( %s, %s, %s, %s);\n" \
         %(str(roundoff(list(hklstars))), str(roundoff(list(colors[i]), 2)),
            str(roundoff(radii[i], 2)), tooltip )


      positions.extend( roundoff(list(hklstars)) )
      colours.extend( roundoff(list(colors[i]), 2) )
      radii2.append( roundoff(radii[i], 2) )
      spbufttips.append(spbufttip)

    # spherebuffer is likely more efficient but how to change alpha values on individual spheres?
    spherebufferstr = """
  ttips = %s
  positions = new Float32Array( %s )
  colours = new Float32Array( %s )
  radii = new Float32Array( %s )
  sphereBuffer = new NGL.SphereBuffer({
    position: positions,
    color: colours,
    radius: radii,
    picking: ttips,
  })
  shape.addBuffer(sphereBuffer)

// create tooltip element and add to the viewer canvas
  tooltip = document.createElement("div");
  Object.assign(tooltip.style, {
    display: "none",
    position: "absolute",
    zIndex: 10,
    pointerEvents: "none",
    backgroundColor: "rgba(0, 0, 0, 0.6)",
    color: "lightgrey",
    padding: "0.5em",
    fontFamily: "sans-serif"
  });

  stage.viewer.container.appendChild(tooltip);
  // listen to `hovered` signal to move tooltip around and change its text
  stage.signals.hovered.add(function (pickingProxy) {
    if (pickingProxy && (pickingProxy.picker == ttips )){
      var sphere = pickingProxy.sphere;
      var cp = pickingProxy.canvasPosition;
      tooltip.innerText = ttips[pickingProxy.pid];
      tooltip.style.bottom = cp.y + 3 + "px";
      tooltip.style.left = cp.x + 3 + "px";
      tooltip.style.display = "block";
    }else{
      tooltip.style.display = "none";
    }
  });

    """ %(str(spbufttips), str(positions), str(colours), str(radii2) )




    self.NGLscriptstr = """


var connection = new WebSocket('ws://127.0.0.1:7894/');

connection.onopen = function () {
  connection.send('Waffle'); // Send the message 'Waffle' to the server
};

// Log errors
connection.onerror = function (error) {
  console.log('WebSocket Error ' + error);
};


window.addEventListener( 'resize', function( event ){
    stage.handleResize();
}, false );


var a=[0,0,0];
var sphereBuffer;
var stage;
var shape;
var shapeComp;
var repr;
var colours = new Float32Array(3)
var sphererepr = new Array( %d )
var spherebufs = new Array( %d )

var hklscene = function (b) {
  shape = new NGL.Shape('shape');
  stage = new NGL.Stage('viewport', { backgroundColor: "grey", tooltip:true });
  stage.setParameters( { cameraType: "%s" } );

  %s

  %s

  shapeComp = stage.addComponentFromObject(shape);
  repr = shapeComp.addRepresentation('buffer');
  shapeComp.autoView();
  repr.update()
}

document.addEventListener('DOMContentLoaded', function() { hklscene(a) }, false );


connection.onmessage = function (e) {
  var retval //= [1,2,-3]
  //alert('Server: ' + e.data);
  val = e.data.split(",")
  //alert("cordinates from server= " + val)
  //document.removeEventListener('DOMContentLoaded', function() { hklscene(a) }, false );
  a = [ val[0], val[1], val[2] ]
  var si = parseInt( val[3] )
  //hklscene(a)
  //shapeComp.removeRepresentation('buffer');

  //shapeComp.dispose()
  //repr.clear()
  //repr.dispose()

  //stage.removeAllComponents();
  //buf = shape.getBufferList()
  //shape._primitiveData.sphereRadius[si] = 2
  //shape._primitiveData.spherePosition[0] = val[0]
  //shape._primitiveData.spherePosition[1] = val[1]
  //shape._primitiveData.spherePosition[2] = val[2]
  //shape._primitiveData.sphereColor[3*si] = (Math.cos(val[0]) +1.0 )/2.0
  //shape._primitiveData.sphereColor[3*si + 1] = 1.0
  //shape._primitiveData.sphereColor[3*si + 2] = 1.0
  //vc1 = shape._primitiveData.sphereColor[3*si]
  //shapeComp.setPosition( val )

  //shapeComp = stage.addComponentFromObject(shape);
  //repr = shapeComp.addRepresentation('buffer');
  //repr.setColor()
  //repr.updateVisibility()
  //repr.update({ color: true })
  //repr.update({ radius: true })
  //repr.setParameters({ color: true })

  var opi = (Math.sin(val[0]) +1.0 )/2.0
  colours[0] = (Math.cos(val[0]) +1.0 )/2.0
  colours[1] = 1.0
  colours[2] = 1.0
  var vc1 = colours[0]

  spherebufs[si].setAttributes({ color: colours })
  //sphererepr[si].setParameters({opacity: opi})
  stage.viewer.requestRender()


  //document.addEventListener('DOMContentLoaded', function() { hklscene(a) }, false );
  connection.send('sphere ' + si + '. Red colour now: ' + vc1); // Send the message 'Waffle' to the server
};




    """ % (nrefls, nrefls, self.cameratype, arrowstr, shapespherestr)
    #""" % (nrefls, nrefls, self.cameratype, arrowstr, spherebufferstr)
    #""" % (nrefls, nrefls, self.cameratype, arrowstr, shapespherestr3)
    if self.jscriptfname:
      with open( self.jscriptfname, "w") as f:
        f.write( self.NGLscriptstr )

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


"""
# python2 code

from websocket_server import WebsocketServer
import time

nc = None
def new_client(client, server):
  #server.send_message(client, u"Oink!")
  server.send_message(client, u"1,2,3" )
  nc = client

def new_message(client, server, message):
  print message

server = WebsocketServer(7894, host='127.0.0.1')
server.set_fn_new_client(new_client)
time.sleep(10)
server.set_fn_message_received(new_message)
server.send_message_to_all(u"Oink!")
server.run_forever()

# python3 code

import asyncio
import datetime
import random
import websockets

async def time(websocket, path):
  x = 0
  for i in range(100):
    x += 0.2
    msg = str(x) + u", 2, 3, " + str(i)
    await websocket.send( msg )
    message = await websocket.recv()
    print( message)
    await asyncio.sleep(1)

start_server = websockets.serve(time, '127.0.0.1', 7894)
asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()



"""





