
// 'use strict';

// Microsoft Edge users follow instructions on
// https://stackoverflow.com/questions/31772564/websocket-to-localhost-not-working-on-microsoft-edge
// to enable websocket connection

var pagename = location.pathname.substring(1);

if ((typeof websocket_portnumber) != "number")
{
// get portnumber for websocket from the number embedded in the filename
// with a regular expression
  var websocket_portnumber = parseInt( pagename.match("hkl_([0-9]+).htm")[1] );
}

if ((typeof websocket_portnumber) != "number")
    alert("Specify port number for the websocket in your html file either like: \n \
<script> var websocket_portnumber = 42673; </script> \n \
or embedded in the filename such as:\n \
C:/Users/Oeffner/AppData/Local/Temp/hkl_42673.htm"
   );

var mysocket;
var socket_intentionally_closed = false;

//var stage = null;
//var shape = null;

var shape = new NGL.Shape('shape');
var stage = new NGL.Stage('viewport', {  backgroundColor: "rgb(128, 128, 128)",
                                    tooltip:false, // create our own tooltip from a div element
                                    fogNear: 100, fogFar: 100 });

var shapeComp = null;
var vectorshape = null;
var repr = null;
var AA = String.fromCharCode(197); // short for angstrom
var DGR = String.fromCharCode(176); // short for degree symbol
var current_ttip = "";
var ttips = [];
var vectorreprs = [];
var vectorshapeComps = [];
var positions = [];
var br_positions = [];
var br_colours = [];
var br_radii = [];
var br_ttips = [];
var colours = [];
var alphas = [];
var radii = [];
var shapebufs = [];
var br_shapebufs = [];
var nrots = 0;
var fontsize = 9;
var postrotmxflag = false;
var cvorient = new NGL.Matrix4();
var oldmsg = "";
var clipFixToCamPosZ = false;
var origclipnear;
var origclipfar;
var origcameraZpos;
var nbins = 0;
var rerendered = false;
var expstate = "";
var current_ttip_ids;
var isdebug = false;
var tdelay = 100;
var displaytooltips = true;
var colourchart = null;
var sockwaitcount = 0;
var ready_for_closing = false;
var columnSelect = null;

var Hstarstart = null;
var Hstarend = null;
var Kstarstart = null;
var Kstarend = null;
var Lstarstart = null;
var Lstarend = null;
var Hlabelpos = null;
var Klabelpos = null;
var Llabelpos = null;


function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}


function createElement(name, properties, style, fsize=10)
{
// utility function used in for loop over colourgradvalarray
  var el = document.createElement(name);
  Object.assign(el, properties);
  Object.assign(el.style, style);
  Object.assign(el.style,
  {
    display: "block",
    position: "absolute",
    fontFamily: "sans-serif",
    fontSize: fsize.toString() + "pt",
  }
  );
  return el;
}


function createSelect(options, properties, style)
{
  var select = createElement("select", properties, style);
  options.forEach(function (d)
  {
    select.add(createElement("option", { value: d[0], text: d[1] }));
  })
  return select;
}


function addElement(el)
{
// utility function used in for loop over colourgradvalarray
  Object.assign(el.style,
  {
    position: "absolute",
    zIndex: 10
  }
  );
  stage.viewer.container.appendChild(el);
}


function addDivBox(name, t, l, w, h, bgcolour="rgba(255, 255, 255, 0.0)", fsize=10)
{
  if (name != null && name != "null")
    txt = name.toString();
  else
    txt = "";

  divbox = createElement("div",
  {
    innerText: txt
  },
  {
    backgroundColor: bgcolour,
    color:  "rgba(0, 0, 0, 1.0)",
    top: t.toString() + "px",
    left: l.toString() + "px",
    width: w.toString() + "px",
    height: h.toString() + "px",
  },
  fsize
  );
  addElement(divbox);
  return divbox;
}


function addDiv2Container(container, name, t, l, w, h, bgcolour="rgba(255, 255, 255, 0.0)", fsize=10)
{
  divbox = addDivBox(name, t, l, w, h, bgcolour, fsize)
  container.append( divbox );
}


function CreateWebSocket()
{
  try
  {
    mysocket = new WebSocket('ws://127.0.0.1:' + websocket_portnumber);
    mysocket.bufferType = "arraybuffer"; // "blob";
    //if (mysocket.readyState !== mysocket.OPEN)
    //  alert('Cannot connect to websocket server! \nAre the firewall permissions or browser security too strict?');
    //  socket_intentionally_closed = false;
    mysocket.onerror = function(e) { onError(e)  };
    mysocket.onopen = function(e) { onOpen(e)  };
    mysocket.onclose = function(e) { onClose(e)  };
    mysocket.onmessage = function(e) { onMessage(e)  };

  }
  catch(err)
  {
    alert('JavaScriptError: ' + err.stack );
    //addDivBox("Error!", window.innerHeight - 50, 20, 40, 20, "rgba(100, 100, 100, 0.0)");
  }
}

CreateWebSocket();


function RemoveStageObjects()
{
  // delete the shapebufs[] that holds the positions[] arrays
  if (shapeComp != null)
  {
    shapeComp.removeRepresentation(repr);
    // remove shapecomp from stage first
    stage.removeAllComponents();
  }
  ttips = [];
  vectorreprs = [];
  vectorshapeComps = [];
  positions = [];
  br_positions = [];
  br_colours = [];
  br_radii = [];
  br_ttips = [];
  expstate = "";
  colours = [];
  alphas = [];
  radii = [];
  shapebufs = [];
  br_shapebufs = [];
  shapeComp = null;
  vectorshape = null;
  repr = null;
  nbins = 0;
}


function WebsockSendMsg(msg, message_is_complete = true)
{
  try
  {
    if (socket_intentionally_closed == true)
      return;
    // Avoid "WebSocket is already in CLOSING or CLOSED state" errors when using QWebEngineView
    // See https://stackoverflow.com/questions/48472977/how-to-catch-and-deal-with-websocket-is-already-in-closing-or-closed-state-in

    if (mysocket.readyState === mysocket.CONNECTING )
    {
      sleep(50).then(()=> {
         WebsockSendMsg(msg);
          return;
        }
      );
    }

    if (mysocket.readyState === mysocket.OPEN)
    {
      mysocket.send(msg);
      if (message_is_complete == true)
        mysocket.send( 'Ready ' + pagename + '\n' );
    }
    else
      if (mysocket.readyState !== mysocket.CONNECTING)
      {
        sleep(200).then(()=> {
            if (mysocket.readyState !== mysocket.OPEN )
            {
              //alert('Closing socket');
              mysocket.close(4242, 'Refreshing ' + pagename); // not sure this is ever received by server
              //window.location.reload(true);
              //alert('Creating socket');
              CreateWebSocket();
              WebsockSendMsg('Connection lost and reestablished')
              WebsockSendMsg(msg);
            }
            return;
          }
        );
        //alert('Cannot send data! \nAre the firewall permissions or browser security too strict?');
      }
  }
  catch(err)
  {
    alert('JavaScriptError: ' + err.stack );
    addDivBox("Error!", window.innerHeight - 50, 20, 40, 20, "rgba(100, 100, 100, 0.0)");
  }
}


var dbgmsg = "";
// debug message window
var debugmessage = document.createElement("div");
Object.assign(debugmessage.style, {
  position: "absolute",
  zIndex: 10,
  pointerEvents: "none",
  backgroundColor: "rgba(255, 255, 255, 0.8 )",
  color: "black",
  padding: "0.1em",
  fontFamily: "sans-serif",
  bottom: "10px",
  left: "10px",
  fontSize: "smaller",
  display: "block"
});



function ReturnClipPlaneDistances()
{
  if (stage.viewer.parameters.clipScale == 'relative')
    cameradist = stage.viewer.cDist;
  if (stage.viewer.parameters.clipScale == 'absolute')
    if (stage.viewer.cDist != 0
         && stage.viewer.parameters.clipFar > stage.viewer.cDist
         && stage.viewer.cDist > stage.viewer.parameters.clipNear)
      cameradist = stage.viewer.cDist;
    else if (stage.viewer.camera.position.z != 0
         && stage.viewer.parameters.clipFar > -stage.viewer.camera.position.z
         && -stage.viewer.camera.position.z > stage.viewer.parameters.clipNear)
      cameradist = stage.viewer.camera.position.z;
    else if (stage.viewer.camera.position.z == -stage.viewer.cDist)
      cameradist = stage.viewer.cDist;
    else
      return;

  msg = String( [stage.viewer.parameters.clipNear,
                  stage.viewer.parameters.clipFar,
                  cameradist ] )
  WebsockSendMsg('ReturnClipPlaneDistances:\n' + msg );
}


async function RenderRequest()
{
  await sleep(100);
  stage.viewer.requestRender();
  WebsockSendMsg( 'RenderRequest ' + pagename );
}

// Log errors to debugger of your browser
function onError(e)
{
  msg = 'WebSocket Error ' + e;
  console.log(msg);
  dbgmsg =msg;
};


function onOpen(e)
{
  msg = 'Now connected via websocket to ' + pagename + '\n';
  WebsockSendMsg(msg);
  dbgmsg =msg;
  rerendered = false;
};


function onClose(e)
{
  msg = 'Now disconnecting from websocket ' + pagename + '\n';
  console.log(msg);
  dbgmsg =msg;
};


function onMessage(e)
{
  var c,
  si;
  WebsockSendMsg('\n    Browser: Got ' + e.data ); // tell server what it sent us
  try
  {
    var datval = e.data.split(":\n");
    var msgtype = datval[0];
    var val = datval[1].split(",");

    if (msgtype === "Reload")
    {
    // refresh browser with the javascript file
      if (stage != null)
      {
        msg = getOrientMsg();
        WebsockSendMsg('OrientationBeforeReload:\n' + msg );
      }
      WebsockSendMsg( 'Refreshing ' + pagename );

      sleep(200).then(()=> {
          socket_intentionally_closed = true;
          mysocket.close(4242, 'Refreshing ' + pagename);
          ready_for_closing = true;
          window.location.reload(true);
          // In 200ms we are gone. A new javascript file will be loaded in the browser
        }
      );
    }

    if (stage == null) // everything below assumes stage!=null
      return;

    if (msgtype === "alpha")
    {
      bin = parseInt(val[0]);
      if (bin < shapebufs.length)
      {
        alphas[bin] = parseFloat(val[1]);
        shapebufs[bin].setParameters({opacity: alphas[bin]});
        for (var g=0; g < nrots; g++ )
          br_shapebufs[bin][g].setParameters({opacity: alphas[bin]});
        //stage.viewer.requestRender();
        RenderRequest();
      }
    }

    if (msgtype === "colour")
    {
      bin = parseInt(val[0]);
      if (bin < shapebufs.length)
      {
        si =  parseInt(val[1]);
        colours[bin][3*si] = parseFloat(val[2]);
        colours[bin][3*si+1] = parseFloat(val[3]);
        colours[bin][3*si+2] = parseFloat(val[4]);
        shapebufs[bin].setAttributes({ color: colours[bin] });

        for (var g=0; g < nrots; g++ )
        {
          br_colours[bin][3*si] = parseFloat(val[2]);
          br_colours[bin][3*si+1] = parseFloat(val[3]);
          br_colours[bin][3*si+2] = parseFloat(val[4]);
          br_shapebufs[bin][g].setAttributes({ color: br_colours[bin] });
        }
        //stage.viewer.requestRender();
        RenderRequest();
      }
    }

    if (msgtype === "DisplayTooltips")
    {
      displaytooltips = val[0];
      stage.signals.hovered.removeAll();
      stage.signals.clicked.removeAll();
      if (displaytooltips == "click")
        stage.signals.clicked.add( PickingProxyfunc );
      if (displaytooltips == "hover")
        stage.signals.hovered.add( PickingProxyfunc );
    }

    if (msgtype === "ShowThisTooltip")
    {
      current_ttip = eval(datval[1]).split("\n\n")[0];
      current_ttip_ids = eval(datval[1]).split("\n\n")[1];
    }

    if (msgtype === "TooltipOpacity")
    {
      Object.assign(tooltip.style, {
        backgroundColor: "rgba(255, 255, 255, " + val[0] + " )",
      });
    }

    if (msgtype === "Redraw")
    {
      RenderRequest();
      WebsockSendMsg( 'Redrawing ' + pagename );
    }

    if (msgtype === "ReOrient")
    {
      WebsockSendMsg( 'Reorienting ' + pagename );
      sm = new Float32Array(16);
      for (j=0; j<16; j++)
      {
        sm[j] = parseFloat(val[j]);
        if (isNaN( sm[j] ))
          return; // do nothing just in case
      }

      var m = new NGL.Matrix4();
      m.fromArray(sm);
      stage.viewerControls.orient(m);
      //stage.viewer.renderer.setClearColor( 0xffffff, 0.01);
      //stage.viewer.requestRender();
      RenderRequest();
      msg = getOrientMsg();
      WebsockSendMsg('CurrentViewOrientation:\n' + msg );
    }

    if (msgtype.includes("Expand") )
    {
      WebsockSendMsg( 'Expanding data...' );

      if (msgtype == "Expand" && expstate == "")
        return;

      if (msgtype == "ExpandP1" && expstate == "isP1Expanded")
        return;

      if (msgtype == "ExpandFriedel" && expstate == "isFriedelExpanded")
        return;

      if (msgtype == "ExpandP1Friedel" && expstate == "isP1FriedelExpanded")
        return;

      // delete the shapebufs[] that holds the positions[] arrays
      shapeComp.removeRepresentation(repr);
      // remove shapecomp from stage first
      stage.removeComponent(shapeComp);

      br_positions = [];
      br_colours = [];
      br_radii = [];
      br_ttips = [];
      br_shapebufs = [];
      var nexpandrefls = 0;

      //alert('rotations:\n' + val);
      // Rotation matrices are concatenated to a string of floats
      // separated by line breaks between each roation matrix
      rotationstrs = datval[1].split("\n");
      var Rotmats = [];
      var r = new NGL.Vector3();

      for (var rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
      {
        Rotmats.push( new NGL.Matrix3() );
        // convert string of rotation matrix elements into a Matrix3
        var elmstrs = rotationstrs[rotmxidx].split(",");
        for (j=0; j<9; j++)
          Rotmats[rotmxidx].elements[j] = parseFloat(elmstrs[j]);
      }

      var Imx = new NGL.Matrix3();
      Imx.identity(); // for testing
      if ( !(msgtype.includes("P1")) && rotationstrs.length == 1 && Rotmats[0].equals(Imx) )
        throw "Only the identity matrix is provided. That means no P1 expansion of reflections!";

      for (var bin=0; bin<nbins; bin++)
      {
        var nsize = positions[bin].length/3; // number of reflections in each bin
        var csize = nsize*3;
        var nsize3 = nsize*3;
        var anoexp = false;

        if (msgtype.includes("Friedel") )
        {
          anoexp = true;
          csize = nsize*6;
        }
        br_positions.push( [] );
        br_shapebufs.push( [] );
        br_colours.push( [] );
        br_radii.push( [] );
        br_ttips.push( [] );

        br_colours[bin] = colours[bin];
        br_radii[bin] = radii[bin];
        if (anoexp)
        {
          var colarr = [];
          var cl = colours[bin].length;
          for (var i=0; i<cl; i++)
          {
            colarr[i] = colours[bin][i];
            colarr[i+cl] = colours[bin][i];
          }
          br_colours[bin] = new Float32Array(colarr);

          var radiiarr = [];
          var rl = radii[bin].length;
          for (var i=0; i<rl; i++)
          {
            radiiarr[i] = radii[bin][i];
            radiiarr[i+rl] = radii[bin][i];
          }
          br_radii[bin] = new Float32Array(radiiarr);
        }

        nrots = 0;
        nexpandrefls = 0;
        for (var rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
        {
          if (rotationstrs[rotmxidx].length < 1 )
            continue;
          nrots++;

          br_positions[bin].push( [] );
          br_shapebufs[bin].push( [] );
          br_ttips[bin].push( [] );
          Object.assign(br_ttips[bin][rotmxidx], ttips[bin]); // deep copy the ttips[bin] object
          br_ttips[bin][rotmxidx][0] = rotmxidx;
          br_positions[bin][rotmxidx] = new Float32Array( csize );
          nexpandrefls += csize;

          for (var i=0; i<nsize; i++)
          {
            idx= i*3;
            r.x = positions[bin][idx];
            r.y = positions[bin][idx+1];
            r.z = positions[bin][idx+2];

            r.applyMatrix3(Rotmats[rotmxidx])

            br_positions[bin][rotmxidx][idx] = r.x;
            br_positions[bin][rotmxidx][idx + 1] = r.y;
            br_positions[bin][rotmxidx][idx + 2] = r.z;

            if (anoexp)
            {
              r.negate(); // inversion for anomalous pair
              br_positions[bin][rotmxidx][nsize3 + idx] = r.x;
              br_positions[bin][rotmxidx][nsize3 + idx + 1] = r.y;
              br_positions[bin][rotmxidx][nsize3 + idx + 2] = r.z;
            }
          }

          br_shapebufs[bin][rotmxidx] = new NGL.SphereBuffer({
              position: br_positions[bin][rotmxidx],
              color: br_colours[bin],
              radius: br_radii[bin],
              // rotmxidx works as the id of the rotation of applied symmetry operator when creating tooltip for an hkl
              picking: br_ttips[bin][rotmxidx],
              } );
          shape.addBuffer(br_shapebufs[bin][rotmxidx]);
          //WebsockSendMsg( 'Memory usage: ' + String(window.performance.memory.totalJSHeapSize) +
          //        ', ' + String(window.performance.memory.totalJSHeapSize) );
        }
        if (nexpandrefls == nsize*3)
          expstate = "";
        if (nexpandrefls == nsize*6)
          expstate = "isFriedelExpanded";
        if (nexpandrefls == nsize*3*nrots && nrots > 1)
          expstate = "isP1Expanded";
        if (nexpandrefls == nsize*6*nrots && nrots > 1)
          expstate = "isP1FriedelExpanded";
      }
      MakeHKL_Axis(shape);

      shapeComp = stage.addComponentFromObject(shape);
      repr = shapeComp.addRepresentation('buffer');

      for (var bin=0; bin<nbins; bin++)
      {
        for (var rotmxidx=0; rotmxidx < nrots; rotmxidx++ )
        {
          br_shapebufs[bin][rotmxidx].setParameters({opacity: alphas[bin]});
        }
      }

      //stage.viewer.requestRender();
      RenderRequest();
      WebsockSendMsg( 'Done ' + msgtype );
    }

    if (msgtype === "DisableMouseRotation")
    {
      WebsockSendMsg( 'Fix mouse rotation' + pagename );
      stage.mouseControls.remove("drag-left");
      stage.mouseControls.remove("scroll-ctrl");
      stage.mouseControls.remove("scroll-shift");
    }

    if (msgtype === "EnableMouseRotation")
    {
      WebsockSendMsg( 'Can mouse rotate ' + pagename );
      stage.mouseControls.add("drag-left", NGL.MouseActions.rotateDrag);
      stage.mouseControls.add("scroll-ctrl", NGL.MouseActions.scrollCtrl);
      stage.mouseControls.add("scroll-shift", NGL.MouseActions.scrollShift);
    }

    if (msgtype === "RotateStage")
    {
      WebsockSendMsg( 'Rotating stage ' + pagename );

      var sm = new Float32Array(9);
      var m4 = new NGL.Matrix4();

      for (j=0; j<9; j++)
        sm[j] = parseFloat(val[j]);

      // GL matrices are the transpose of conventional rotation matrices
      m4.set( sm[0], sm[3], sm[6], 0.0,
              sm[1], sm[4], sm[7], 0.0,
              sm[2], sm[5], sm[8], 0.0,
              0.0,   0.0,   0.0,   1.0
      );
      stage.viewerControls.orient(m4);
      if (val[9]=="verbose")
        postrotmxflag = true;
      ReturnClipPlaneDistances();
      RenderRequest();
      sleep(100).then(()=> {
          msg = getOrientMsg();
          WebsockSendMsg('CurrentViewOrientation:\n' + msg );
        }
      );
    }

    if (msgtype === "SpinAnimate")
    {
      WebsockSendMsg( 'SpinAnimating ' + pagename );
      //strs = datval[1].split("\n");
      var r = new Float32Array(3);
      //var elmstrs = strs[0].split(",");
      for (j=0; j<3; j++)
        r[j] = parseFloat(val[j]);
      if (r[0] == 0.0 && r[1] == 0.0 && r[2] == 0.0)
      {
        // default bindings as per ngl\src\controls\mouse-actions.ts
        stage.mouseControls.add("drag-ctrl-left", NGL.MouseActions.panDrag);
        stage.mouseControls.add("drag-ctrl-right", NGL.MouseActions.focusScroll);
        stage.mouseControls.add("drag-shift-left", NGL.MouseActions.zoomDrag);
        stage.mouseControls.add("drag-shift-right", NGL.MouseActions.zoomDrag);
        stage.mouseControls.add("drag-middle", NGL.MouseActions.zoomDrag);
        stage.mouseControls.add("drag-right", NGL.MouseActions.panDrag);
        stage.mouseControls.add("drag-left", NGL.MouseActions.rotateDrag);
        stage.mouseControls.add("scroll-ctrl", NGL.MouseActions.scrollCtrl);
        stage.mouseControls.add("scroll-shift", NGL.MouseActions.scrollShift);
        stage.setSpin(false);
      }
      else
      {
        stage.spinAnimation.axis.set(r[0], r[1], r[2]);

        stage.mouseControls.remove("drag-ctrl-left");
        stage.mouseControls.remove("drag-ctrl-right");
        stage.mouseControls.remove("drag-shift-left");
        stage.mouseControls.remove("drag-shift-right");
        stage.mouseControls.remove("drag-middle");
        stage.mouseControls.remove("drag-right");
        stage.mouseControls.remove("drag-left");
        stage.mouseControls.remove("scroll-ctrl");
        stage.mouseControls.remove("scroll-shift");
        stage.setSpin(true);
      }
    }

    if (msgtype === "TranslateHKLpoints")
    {
      WebsockSendMsg( 'Translating HKLs ' + pagename );
      strs = datval[1].split("\n");
      var sm = new Float32Array(3);
      var elmstrs = strs[0].split(",");
      for (j=0; j<3; j++)
        sm[j] = parseFloat(elmstrs[j]);
      shapeComp.setPosition([ sm[0], sm[1], sm[2] ]);
      //stage.viewer.requestRender();
      RenderRequest();
      sleep(100).then(()=> {
          msg = getOrientMsg();
          WebsockSendMsg('CurrentViewOrientation:\n' + msg );
        }
      );
    }

    function DeleteVectors(reprname)
    {
      thisrepr = stage.getRepresentationsByName(reprname);
      var wasremoved = false;
      for (i=0; i<stage.compList.length; i++)
        if (stage.compList[i].reprList[0].name == reprname)
        {
          thiscomp = stage.compList[i];
          thiscomp.removeRepresentation(thisrepr);
          stage.removeComponent(thiscomp);
          wasremoved = true;
        }
      return wasremoved;
    };

    if (msgtype === "AddVector")
    {
      var r1 = new Float32Array(3);
      var r2 = new Float32Array(3);
      var rgb = new Float32Array(3);
      for (j=0; j<3; j++)
      {
        r1[j] = parseFloat(val[j]);
        r2[j] = parseFloat(val[j+3]);
        rgb[j]= parseFloat(val[j+6]);
      }
      radius = parseFloat(val[11]);

      if (vectorshape == null)
        vectorshape = new NGL.Shape('vectorshape');

      vectorshape.addArrow( r1, r2 , [rgb[0], rgb[1], rgb[2]], radius);
      if (val[6] !== "") {
        var txtR = [ (r1[0] + r2[0])/2.0, (r1[1] + r2[1])/2.0, (r1[2] + r2[2])/2.0 ];
        vectorshape.addText( txtR, [rgb[0], rgb[1], rgb[2]], fontsize, val[9] );
      }
      // if reprname is supplied with a vector then make a representation named reprname
      // of this and all pending vectors stored in vectorshape and render them.
      // Otherwise just accummulate the new vector
      var reprname = val[10].trim();
      if (reprname != "")
      {
        DeleteVectors(reprname); // delete any existing vectors with the same name
        vectorshapeComps.push( stage.addComponentFromObject(vectorshape) );
        vectorreprs.push(
          vectorshapeComps[vectorshapeComps.length-1].addRepresentation('vecbuf',
                                                                      { name: reprname} )
        );
        vectorshape = null;
        //stage.viewer.requestRender();
        RenderRequest();
      }
    }

    if (msgtype === "RemoveVectors")
    {
      var reprname = val[0].trim(); // elmstrs[0].trim();
      // if reprname is supplied only remove vectors with that name
      var reprnamegone = false;
      var clipvecgone = false;
      var unitcellgone = false;
      var reciprocunitcellgone = false;
      if (reprname != "")
        reprnamegone = DeleteVectors(reprname);
      else // otherwise remove all vectors
      {
        clipvecgone = DeleteVectors("clip_vector");
        unitcellgone = DeleteVectors("unitcell");
        reciprocunitcellgone = DeleteVectors("reciprocal_unitcell");
      }
      if (reprnamegone || clipvecgone || unitcellgone || reciprocunitcellgone)
        //stage.viewer.requestRender();
        RenderRequest();
    }

    if (msgtype === "DefineHKL_Axes")
    {
      strarrs = datval[1].split("\n\n");
      hstart = eval(strarrs[0]);
      hend = eval(strarrs[1]);
      kstart = eval(strarrs[2]);
      kend = eval(strarrs[3]);
      lstart = eval(strarrs[4]);
      lend = eval(strarrs[5]);
      hlabelpos = eval(strarrs[6]);
      klabelpos = eval(strarrs[7]);
      llabelpos = eval(strarrs[8]);

      DefineHKL_Axes(hstart, hend, kstart, kend, 
                 lstart, lend, hlabelpos, klabelpos, llabelpos)
    }

    if (msgtype === "SetFontSize")
    {
      fontsize = parseFloat(val[0]);
      RenderRequest();
    }

    if (msgtype === "SetMouseSpeed")
    {
      stage.trackballControls.rotateSpeed = parseFloat(val[0]);
    }

    if (msgtype === "SetCameraType")
    {
      camtype = val[0];
      stage.setParameters( { cameraType: camtype } );
      RenderRequest();
    }

    if (msgtype === "GetMouseSpeed")
    {
      msg = String( [stage.trackballControls.rotateSpeed] )
      WebsockSendMsg('ReturnMouseSpeed:\n' + msg );
    }

    if (msgtype === "SetClipPlaneDistances")
    {
      var near = parseFloat(val[0]);
      var far = parseFloat(val[1]);
      origcameraZpos = parseFloat(val[2]);
      stage.viewer.parameters.clipMode =  'camera';
      // clipScale = 'absolute' means clip planes are using scene dimensions
      stage.viewer.parameters.clipScale = 'absolute';
      clipFixToCamPosZ = true;

      if (near >= far )
      { // default to no clipping if near >= far
        stage.viewer.parameters.clipMode = 'scene';
      // clipScale = 'relative' means clip planes are in percentage
        stage.viewer.parameters.clipScale = 'relative';
        clipFixToCamPosZ = false;
        near = 0;
        far = 100;
      }
      else
        stage.viewer.camera.position.z = origcameraZpos;
      stage.viewer.parameters.clipNear = near;
      stage.viewer.parameters.clipFar = far;
      origclipnear = near;
      origclipfar = far;
      //stage.viewer.requestRender();
      RenderRequest();
    }

    if (msgtype === "GetClipPlaneDistances")
      ReturnClipPlaneDistances();

    if (msgtype === "GetBoundingBox")
    {
      msg = String( [stage.viewer.boundingBoxSize.x,
                     stage.viewer.boundingBoxSize.y,
                     stage.viewer.boundingBoxSize.z]
                  )
      WebsockSendMsg('ReturnBoundingBox:\n' + msg );
    }

    if (msgtype ==="JavaScriptCleanUp")
    {
      RemoveStageObjects();
      stage.mouseObserver.dispose();
      stage.dispose();
      stage = null;
      ready_for_closing = true;
      WebsockSendMsg('JavaScriptCleanUpDone:\nDestroying JavaScript objects');
      socket_intentionally_closed = true;
      mysocket.close(4241, 'Cleanup done');
      document = null;
    }

    if (msgtype ==="SetBrowserDebug")
    {
      isdebug = (val[0] === "true" );
    }

    if (msgtype ==="RemoveStageObjects")
    {
      RemoveStageObjects();
    }

    if (msgtype === "AddSpheresBin2ShapeBuffer")
    {
      strarrs = datval[1].split("\n\n");
      coordarray = eval(strarrs[0]);
      colourarray = eval(strarrs[1]);
      radiiarray = eval(strarrs[2]);
      ttipids = eval(strarrs[3]);
      AddSpheresBin2ShapeBuffer(coordarray, colourarray, radiiarray, ttipids);
    }

    if (msgtype === "MakeColourChart")
    {
      msg = datval[1].split("\n\n");
      ctop = eval(msg[0]);
      cleft = eval(msg[1]);
      label = msg[2];
      fomlabel = msg[3];
      colourgradvalarrays = eval(msg[4]);
      MakeColourChart(ctop, cleft, label, fomlabel, colourgradvalarrays);
      RenderRequest();
    }

    if (msgtype ==="RenderStageObjects")
    {
      //HKLscene();
      MakeHKL_Axis(shape);
      shapeComp = stage.addComponentFromObject(shape);
      repr = shapeComp.addRepresentation('buffer');
      RenderRequest();
      WebsockSendMsg('Drawing new reflections');
    }

    if (msgtype === "SetAutoView")
    {
      if (shapeComp != null) // workaround for QTWebEngine bug sometimes failing to render scene
        shapeComp.autoView();
      WebsockSendMsg('AutoViewSet ' + pagename);
    }

    if (msgtype === "MakeImage")
    {
      filename = val[0];
      stage.viewer.makeImage( {
                factor: 1,
                antialias: true,
                trim: false,
                transparent: false
            } ).then( function( blob ){
              if (parseInt(val[1]) < 3)
              {
// Using websocket_server in python2 which doesn't allow streaming large compressed data
// So use NGL's download image function
                NGL.download( blob, filename );
              }
              else
              { // websockets in python3 which supports streaming large blobs
                WebsockSendMsg('Imageblob', false);
                WebsockSendMsg( blob );
              }

              WebsockSendMsg('ImageWritten ' + pagename);
        } );
    }

    if (msgtype === "MakeBrowserDataColumnComboBox")
    {
      if (columnSelect != null)
        columnSelect.remove(); 

      msg = datval[1].split("\n\n");
      var columnSelect = createElement("select", {
        onchange: function (e)
        {
          WebsockSendMsg('SelectedBrowserDataColumnComboBox: ' + e.target.value);
        }
      }, { top: "25px", right: "10px", width: "130px", position: "absolute" }, fsize = fontsize);

      for (i = 0; i < msg.length; i++)
      {
        labelval = msg[i].split("\n");
        columnSelect.add(createElement("option", { text: labelval[0], value: labelval[1] }, fsize = fontsize));
      }
      addElement(columnSelect);

      // 
      divlabel = createElement("div",
        {
          innerText: "Select Data"
        },
        {
          backgroundColor: "rgba(255, 255, 255, 1.0)",
          color: "rgba(0, 0, 0, 1.0)",
          top: "10px", right: "10px", width: "130px",
          position: "absolute"
        },
        fsize = fontsize
      );
      addElement(divlabel);

    }

    if (msgtype === "Testing")
    {
      // test something new
      /*
      var newradii = radii[0].map(function(element) {
        return element*1.5;
      });
      shapebufs[0].setAttributes({
          radius: newradii
      })
      repr = shapeComp.addRepresentation('buffer');
      //stage.viewer.requestRender();
      RenderRequest();
      */
    }
    WebsockSendMsg('Received message: ' + msgtype );
    if (isdebug)
      debugmessage.innerText = dbgmsg;
  }

  catch(err)
  {
    WebsockSendMsg('JavaScriptError: ' + err.stack );
  }

};


var ttipalpha = 0.7;
var camtype = "orthographic";
var negativeradiistr

function timefunc() {
  var d = new Date();
  var now = d.getTime();
  return now
}

var timenow = timefunc();
var rightnow = timefunc();


window.addEventListener( 'resize',
  function( event )
  {
    stage.handleResize();
  },
  false
);


window.onbeforeunload = function(event) 
{
  if (!ready_for_closing)
    WebsockSendMsg('Warning!: Web browser closed unexpectedly perhaps by an external process. Call JavaScriptCleanUp() or Reload() instead.')
};


if (isdebug)
{
  var script=document.createElement('script');
  script.src='https://rawgit.com/paulirish/memory-stats.js/master/bookmarklet.js';
  document.head.appendChild(script);
}



// define tooltip element
var tooltip = document.createElement("div");
Object.assign(tooltip.style, {
  display: "none",
  position: "absolute",
  zIndex: 10,
  pointerEvents: "none",
  backgroundColor: "rgba(255, 255, 255, ttipalpha )",
  color: "black",
  padding: "0.1em",
  fontFamily: "sans-serif"
});


function DefineHKL_Axes(hstart, hend, kstart, kend, 
                 lstart, lend, hlabelpos, klabelpos, llabelpos)
{
  Hstarstart = hstart;
  Hstarend = hend;
  Kstarstart = kstart;
  Kstarend = kend;
  Lstarstart = lstart;
  Lstarend = lend;
  Hlabelpos = hlabelpos;
  Klabelpos = klabelpos;
  Llabelpos = llabelpos;
};


function MakeHKL_Axis()
{
  // xyz arrows
  // shape.addSphere( [0,0,0] , [ 1, 1, 1 ], 0.3, 'Origin');
  //blue-x
  shape.addArrow( Hstarstart, Hstarend , [ 0, 0, 1 ], 0.1);
  //green-y
  shape.addArrow( Kstarstart, Kstarend, [ 0, 1, 0 ], 0.1);
  //red-z
  shape.addArrow( Lstarstart, Lstarend, [ 1, 0, 0 ], 0.1);

  shape.addText( Hlabelpos, [ 0, 0, 1 ], fontsize, 'h');
  shape.addText( Klabelpos, [ 0, 1, 0 ], fontsize, 'k');
  shape.addText( Llabelpos, [ 1, 0, 0 ], fontsize, 'l');
};


function getOrientMsg()
{
  cvorientmx = stage.viewerControls.getOrientation();
  if (cvorientmx.determinant() == 0)
      return oldmsg; // don't return invalid matrix

  cvorient = cvorientmx.elements;
  for (j=0; j<16; j++)
  {
    if (Number.isNaN( cvorient[j]) )
      return oldmsg; // don't return invalid matrix
  }

  if (stage.viewer.cDist != 0
        && stage.viewer.parameters.clipFar > stage.viewer.cDist
        && stage.viewer.cDist > stage.viewer.parameters.clipNear)
    cameradist = stage.viewer.cDist;
  else if (stage.viewer.camera.position.z != 0
        && stage.viewer.parameters.clipFar > -stage.viewer.camera.position.z
        && -stage.viewer.camera.position.z > stage.viewer.parameters.clipNear)
    cameradist = -stage.viewer.camera.position.z;
  else
    cameradist = cvorient[14]; // fall back if stage.viewer.camera.position.z is corrupted
  cvorient.push( cameradist );
  msg = String(cvorient);
  oldmsg = msg;
  return msg;
}


// listen to hover or click signal to move tooltip around and change its text
function PickingProxyfunc(pickingProxy)
{
// adapted from http://nglviewer.org/ngl/api/manual/interaction-controls.html#clicked
  if (pickingProxy
        && (Object.prototype.toString.call(pickingProxy.picker["ids"]) === '[object Array]' )
        && displaytooltips )
  {
    var cp = pickingProxy.canvasPosition;
    var sym_id = -1;
    var hkl_id = -1;
    var ttipid = "";
    if (pickingProxy.picker["ids"].length > 0)
    { // get stored id number of symmetry operator applied to this hkl
      sym_id = pickingProxy.picker["ids"][0];
      var ids = pickingProxy.picker["ids"].slice(1);
      var is_friedel_mate = 0;
      hkl_id = ids[ pickingProxy.pid % ids.length ];
      if (pickingProxy.pid >= ids.length)
        is_friedel_mate = 1;
    }
    // tell python the id of the hkl and id number of the symmetry operator
    rightnow = timefunc();
    if (rightnow - timenow > tdelay)
    { // only post every 50 milli second as not to overwhelm python
      ttipid = String([hkl_id, sym_id, is_friedel_mate]);
      WebsockSendMsg( 'tooltip_id: [' + ttipid + ']' );
      timenow = timefunc();
    }

    if (isdebug)
      console.log( "current_ttip_ids: " + String(current_ttip_ids) + ", ttipid: " + String(ttipid) );
    if (current_ttip !== "" && current_ttip_ids == ttipid )
    {
      tooltip.innerText = current_ttip;
      tooltip.style.bottom = cp.y + 7 + "px";
      tooltip.style.left = cp.x + 8 + "px";
      tooltip.style.fontSize = fontsize.toString() + "pt";
      tooltip.style.display = "block";
    }
  }
  else
  {
    tooltip.style.display = "none";
    current_ttip = "";
  }
};


function getTextWidth(text, fsize=8)
{
  // re-use canvas object for better performance
  var canvas = getTextWidth.canvas || (getTextWidth.canvas = document.createElement("canvas"));
  var context = canvas.getContext("2d");
  context.font = fsize.toString() + "pt sans-serif";
  var metrics = context.measureText(text);
  return metrics.width;
}


function MakeColourChart(ctop, cleft, millerlabel, fomlabel, colourgradvalarrays)
{
  /* colourgradvalarrays is a list of colour charts. If only one list then it's one colour chart.
  Otherwise it's usually a list of colour charts that constitute a gradient across colours,
  typically used for illustrating figure of merits attenuating phase values in map coefficients
  */
  var ih = 3,
  topr = 25,
  topr2 = 0,
  lp = 10;

  var maxnumberwidth = 0;
  for (j = 0; j < colourgradvalarrays[0].length; j++)
  {
    val = colourgradvalarrays[0][j][0];
    maxnumberwidth = Math.max( getTextWidth(val, fontsize), maxnumberwidth );
  }
  wp = maxnumberwidth + 5,
  //wp = 60,
  lp2 = lp + wp,
  gl = 3,
  wp2 = gl,
  fomlabelheight = 25;

  if (colourgradvalarrays.length === 1)
  {
    wp2 = 15;
    fomlabelheight = 0;
  }
  var wp3 = wp + colourgradvalarrays.length * wp2 + 2;

  totalheight = ih*colourgradvalarrays[0].length + 35 + fomlabelheight;

  if (colourchart != null)
    colourchart.remove(); // delete previous colour chart if any
  colourchart = addDivBox(null, ctop, cleft, wp3, totalheight, bgcolour="rgba(255, 255, 255, 1.0)");

  // make a white box on top of which boxes with transparent background are placed
  // containing the colour values at regular intervals as well as label legend of
  // the displayed miller array
  addDiv2Container(colourchart, null, topr2, lp, wp3, totalheight, 'rgba(255, 255, 255, 1.0)');

  // print label of the miller array used for colouring
  lblwidth = getTextWidth(millerlabel, fontsize);
  addDiv2Container(colourchart, millerlabel, topr2, lp, lblwidth + 5, 20, 'rgba(255, 255, 255, 1.0)', fsize=fontsize);

  if (fomlabel != "" )
  {
    // print FOM label, 1, 0.5 and 0.0 values below colour chart
    fomtop = topr2 + totalheight - 18;
    fomlp = lp + wp;
    fomwp = wp3;
    fomtop2 = fomtop - 13;
    // print the 1 number
    addDiv2Container(colourchart, 1, fomtop2, fomlp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fsize=fontsize);
    // print the 0.5 number
    leftp = fomlp + 0.48 * gl * colourgradvalarrays.length;
    addDiv2Container(colourchart, 0.5, fomtop2, leftp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fsize=fontsize);
    // print the FOM label
    addDiv2Container(colourchart, fomlabel, fomtop, fomlp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fsize=fontsize);
    // print the 0 number
    leftp = fomlp + 0.96 * gl * colourgradvalarrays.length;
    addDiv2Container(colourchart, 0, fomtop2, leftp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fsize=fontsize);
  }

  for (j = 0; j < colourgradvalarrays[0].length; j++)
  {
    val = colourgradvalarrays[0][j][0];
    topv = j*ih + topr;
    toptxt = topv - 5;
    // print value of miller array if present in colourgradvalarrays[0][j][0]
    addDiv2Container(colourchart,val, toptxt, lp, wp, ih, 'rgba(255, 255, 255, 0.0)', fsize=fontsize);
  }

  // if colourgradvalarrays is an array of arrays then draw each array next to the previous one
  for (g = 0; g < colourgradvalarrays.length; g++)
  {
    leftp = g*gl + lp + wp;
    for (j = 0; j < colourgradvalarrays[g].length; j++)
    {
      R = colourgradvalarrays[g][j][1];
      G = colourgradvalarrays[g][j][2];
      B = colourgradvalarrays[g][j][3];
      rgbcol = 'rgba(' + R.toString() + ',' + G.toString() + ',' + B.toString() + ', 1.0)'
      topv = j*ih + topr;
      addDiv2Container(colourchart, null, topv, leftp, wp2, ih, rgbcol);
    }
  }

  colourchart.onselectstart = function ()
  { // don't select numbers or labels on chart when double clicking the coulour chart
    return false;
  }

  colourchart.oncontextmenu = function (e) 
  { // oncontextmenu captures right clicks
    e.preventDefault()
    alert("in oncontextmenu")
    return false;
  };

  colourchart.ondblclick = function (e)
  {
    var sel = window.getSelection();
    sel.removeAllRanges(); // don't select numbers or labels on chart when double clicking the coulour chart
    WebsockSendMsg('doubleclick colour chart');
  };
}


function AddSpheresBin2ShapeBuffer(coordarray, colourarray, radiiarray, ttipids) 
{
  ttiplst = [-1].concat(ttipids);
  ttips.push( { ids: ttiplst,
       getPosition: function() { return { x:0, y:0 }; } // dummy function to avoid crash
  }  );
  positions.push( new Float32Array( coordarray ) );
  colours.push( new Float32Array( colourarray ) );
  radii.push( new Float32Array( radiiarray ) );
  curridx = positions.length -1;
  shapebufs.push( new NGL.SphereBuffer({
    position: positions[curridx],
    color: colours[curridx], 
    radius: radii[curridx],
    picking: ttips[curridx],
    })
  );
  shape.addBuffer(shapebufs[curridx]);
  alphas.push(1.0);
  nbins = nbins + 1;
}


function HKLscene()
{
  shape = new NGL.Shape('shape');
  stage = new NGL.Stage('viewport', {  backgroundColor: "rgb(128, 128, 128)",
                                    tooltip:false, // create our own tooltip from a div element
                                    fogNear: 100, fogFar: 100 });

  stage.setParameters( { cameraType: camtype } );

  //MakeHKL_Axis(shape);

//  placeholder for spherebufferstr
  

// create tooltip element and add to the viewer canvas
  stage.viewer.container.appendChild(tooltip);
  stage.signals.clicked.add( PickingProxyfunc );


  stage.mouseObserver.signals.dragged.add(
    function ( deltaX, deltaY)
    {
      if (clipFixToCamPosZ === true)
      {
        stage.viewer.parameters.clipNear = origclipnear + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = origclipfar + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.requestRender();
      }
      msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        postrotmxflag = true;
        WebsockSendMsg('CurrentViewOrientation:\n' + msg );
        timenow = timefunc();
      }
    }
  );


  stage.mouseObserver.signals.clicked.add(
    function (x, y)
    {
      msg = getOrientMsg();
      WebsockSendMsg('CurrentViewOrientation:\n' + msg );
    }
  );


  stage.mouseObserver.signals.scrolled.add(
    function (delta)
    {
      if (clipFixToCamPosZ === true)
      {
        stage.viewer.parameters.clipNear = origclipnear + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = origclipfar + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.requestRender();
      }
      msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        WebsockSendMsg('CurrentViewOrientation:\n' + msg );
        timenow = timefunc();
      }
    }
  );


  stage.viewer.signals.rendered.add(
    function()
    {
      if (postrotmxflag === true)
      {
        msg = getOrientMsg();
        WebsockSendMsg('CurrentViewOrientation:\n' + msg );
        postrotmxflag = false;
      }
    }
  );


  stage.viewerControls.signals.changed.add(
    function()
    {
      msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        WebsockSendMsg('CurrentViewOrientation:\n' + msg );
        //ReturnClipPlaneDistances();
        sleep(250).then(()=> {
            ReturnClipPlaneDistances();
          }
        );
        timenow = timefunc();
      }
    }
  );

  shapeComp = stage.addComponentFromObject(shape);
  repr = shapeComp.addRepresentation('buffer');
  shapeComp.autoView();
  repr.update();

  if (isdebug)
    stage.viewer.container.appendChild(debugmessage);

  // avoid NGL zoomFocus messing up clipplanes positions. So reassign those signals to zoomDrag
  stage.mouseControls.remove("drag-shift-right");
  stage.mouseControls.add("drag-shift-right", NGL.MouseActions.zoomDrag);
  stage.mouseControls.remove("drag-middle");
  stage.mouseControls.add("drag-middle", NGL.MouseActions.zoomDrag);
  stage.mouseControls.remove('clickPick-left'); // avoid undefined move-pick when clicking on a sphere

  stage.viewer.requestRender();
  if (isdebug)
    debugmessage.innerText = dbgmsg;
}


function OnUpdateOrientation()
{
  msg = getOrientMsg();
  WebsockSendMsg('MouseMovedOrientation:\n' + msg );
}


function PageLoad()
{
  try
  {
    document.addEventListener('DOMContentLoaded', function() { HKLscene() }, false );
    document.addEventListener('mouseup', function() { OnUpdateOrientation() }, false );
    document.addEventListener('wheel', function(e) { OnUpdateOrientation() }, false );
    document.addEventListener('scroll', function(e) { OnUpdateOrientation() }, false );
    // mitigate flickering on some PCs when resizing
    document.addEventListener('resize', function () { RenderRequest() }, false);


  }
  catch(err)
  {
    WebsockSendMsg('JavaScriptError: ' + err.stack );
  }
}


PageLoad();
