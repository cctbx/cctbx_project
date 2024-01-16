// HKLviewer driver file leveraging NGL.js accessed from python via a websocket 

'use strict';

// Microsoft Edge users follow instructions on
// https://stackoverflow.com/questions/31772564/websocket-to-localhost-not-working-on-microsoft-edge
// to enable websocket connection

if ((typeof isHKLviewer) != "boolean")
  var isHKLviewer = false;

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

var shape = new NGL.Shape('shape');
var stage = new NGL.Stage('viewport', {  backgroundColor: "rgb(128, 128, 128)",
                                    tooltip:false, // create our own tooltip from a div element
                                    fogNear: 100, fogFar: 100 });

var shapeComp = null;
var vectorshape = null;
var repr = null;
const AA = String.fromCharCode(197); // short for angstrom
const DGR = String.fromCharCode(176); // short for degree symbol
var current_ttip = "";
var ttips = [];
var filename = ""
var vectorreprs = [];
var vectorshapeComps = [];
var positions = [];
var expansion_positions = [];
var expansion_colours = [];
var expansion_radii = [];
var expansion_ttips = [];
var colours = [];
var alphas = [];
var radii = [];
var shapebufs = [];
var expansion_shapebufs = [];
var nrots = 0;
var fontsize = 9;
var vectorwidth = 2
var postrotmxflag = false;
var cvorient = new NGL.Matrix4();
var oldmsg = "";
var binmsgtype = "";
var clipFixToCamPosZ = false;
var nbins = 0;
var rerendered = false;
var expstate = "";
var current_ttip_ids;
var isdebug = false;
var tdelay = 200;
var ttipid = "";
var displaytooltips = true;
var colourchart = null;
var millerlabel = null;
var fomlabel = null;
var colourgradvalarrays = null;
var infobanner = null;
var ResetViewBtn = null;
var PlusBtn = null;
var MinusBtn = null;
var pmleft = null;
var pmbottom = null;
var btnwidth = null;
var StopAnimateBtn = null;
var maxReflInFrustum = 0;
var hklequationmsg = "";
var animatheta = 0.0;
var animaaxis = new NGL.Vector3();
var componenttheta = 0.0;
var componentaxis = new NGL.Vector3();
var sockwaitcount = 0;
var ready_for_closing = false;
var columnSelect = null;
var animationspeed = 0.0;
var XYZaxes = null;
var Helm = null;
var Kelm = null;
var Lelm = null;
var Hstarstart = null;
var Hstarend = null;
var Kstarstart = null;
var Kstarend = null;
var Lstarstart = null;
var Lstarend = null;
var Hlabelpos = null;
var Klabelpos = null;
var Llabelpos = null;
var Hlabelvec = new NGL.Vector3();
var Klabelvec = new NGL.Vector3();
var Llabelvec = new NGL.Vector3();
var annodivs = [];
var rotationdisabled = false;
var div_annotation_opacity = 0.7;
var camtype = "orthographic";
var canvaspos = null;
var requestedby = "";
var isAutoviewing = false;
var zoomanis = null;


const simulate_click_evt = new MouseEvent("simulate_click", {
  view: window,
  bubbles: true,
  cancelable: true,
  clientX: 20,
});


function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}


function createElement(name, properties, style, fsize=10)
{
// utility function used in for loop over colourgradvalarray
  let el = document.createElement(name);
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
  let select = createElement("select", properties, style);
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


function addDivBox(name, t, l, w, h, bgcolour = "rgba(255, 255, 255, 0.0)", fsize = 10) {
  let txt = "";
  if (name != null && name != "null")
    txt = name.toString();
  else
    txt = "";

  let divbox = createElement("div",
    {
      innerText: txt
    },
    {
      backgroundColor: bgcolour,
      color: "rgba(0, 0, 0, 1.0)",
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


function addBottomDivBox(name, b, l, w, h, bgcolour = "rgba(255, 255, 255, 0.0)", fsize = 10) {
  let txt = "";
  if (name != null && name != "null")
    txt = name.toString();
  else
    txt = "";

  let divbox = createElement("div",
    {
      innerText: txt
    },
    {
      backgroundColor: bgcolour,
      color: "rgba(0, 0, 0, 1.0)",
      bottom: b.toString() + "px",
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
  let divbox = addDivBox(name, t, l, w, h, bgcolour, fsize)
  container.append( divbox );
}


function createDivElement(label, rgb)
{
  let elm = createElement("div", { innerText: label },
    {
      color: "rgba(255, 255, 255, 255, 1.0)",    
      backgroundColor: "rgba(" + rgb[0] + ", " + rgb[1] + ", " + rgb[2] + ", " + div_annotation_opacity + ")",
      padding: "4px"
    }, fontsize
  );
  
  let pointelm = createElement("pointdiv", // make a small white square to indicate where this elm is pointing
    {  innerText: ""  },
    {
      backgroundColor: "rgba(255, 255, 255, 0.0)",
      color: "rgba(255, 255, 255, 1.0)",
      top: "0px",
      left: "0px",
      width: "2px",
      height: "2px",
    },
    fontsize
  );
  
  elm.append(pointelm);
  return elm;
}


function CreateWebSocket()
{
  try
  {
    mysocket = new WebSocket('ws://localhost:' + websocket_portnumber + '/');
    mysocket.binaryType = "arraybuffer"; // "blob";
    mysocket.onerror = function(e) { onError(e)  };
    mysocket.onopen = function(e) { onOpen(e)  };
    mysocket.onclose = function(e) { onClose(e)  };
    mysocket.onmessage = function(e) { onMessage(e)  };
  }
  catch(err)
  {
    alert('JavaScriptError: ' + err.stack );
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
  if (colourchart != null) {
    colourchart.remove(); // delete previous colour chart if any
    colourchart = null;
  }
  ttips = [];
  vectorreprs = [];
  vectorshapeComps = [];
  positions = [];
  expansion_positions = [];
  expansion_colours = [];
  expansion_radii = [];
  expansion_ttips = [];
  expstate = "";
  colours = [];
  alphas = [];
  radii = [];
  shapebufs = [];
  expansion_shapebufs = [];
  shapeComp = null;
  // delete shape to ensure shape.boundingbox will equal viewer.boundingbox of the currently loaded reflections
  shape = null;
  vectorshape = null;
  repr = null;
  nbins = 0;
  tooltip.style.display = "none";
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
              mysocket.close(4242, 'Refreshing ' + pagename); // not sure this is ever received by server
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


function ReturnClipPlaneDistances(calledby = "")
{ // If onrequest=true then jsview.py will increase the semaphore 
  // count by calling clipplane_msg_sem.release().
  // Only do this in response to a explicit message like "GetClipPlaneDistances"
  // from where the sempahore has been acquired
  let cameradist = stage.viewer.camera.position.z;
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

  let msg = String([stage.viewer.parameters.clipNear, stage.viewer.parameters.clipFar,
    cameradist, stage.viewer.camera.zoom, calledby]);
  WebsockSendMsg('ReturnClipPlaneDistances:\n' + msg );
}


function DeletePrimitives(reprname)
{
  let thisrepr = stage.getRepresentationsByName(reprname);
  let wasremoved = false;
  for (let i = 0; i < stage.compList.length; i++)
    if (stage.compList[i].reprList.length > 0 && stage.compList[i].reprList[0].name == reprname) {
      let thiscomp = stage.compList[i];
      thiscomp.removeAllAnnotations();
      thiscomp.removeRepresentation(thisrepr);
      stage.removeComponent(thiscomp);
      wasremoved = true;
    }
  return wasremoved;
};

function RemovePrimitives(reprname)
{
  // if reprname is supplied only remove vectors with that name
  let reprnamegone = false;
  let clipvecgone = false;
  let unitcellgone = false;
  let symHKLsgone = false;
  let reciprocunitcellgone = false;
  if (reprname != "")
    reprnamegone = DeletePrimitives(reprname);
  else // otherwise remove all vectors
  {
    symHKLsgone = DeletePrimitives("sym_HKLs");
    clipvecgone = DeletePrimitives("clip_vector");
    unitcellgone = DeletePrimitives("unitcell");
    reciprocunitcellgone = DeletePrimitives("reciprocal_unitcell");
  }
  if (reprnamegone || clipvecgone || unitcellgone || reciprocunitcellgone || symHKLsgone)
    RenderRequest();
}


function CameraZoom(t, deltaX, deltaY) {
  let dx = 0;
  if (Number.isInteger(deltaX))
    dx = deltaX;
  let dy = 0;
  if (Number.isInteger(deltaY))
    dy = deltaY;
  let z = stage.viewer.camera.zoom;
  z -= (dx + dy) * stage.trackballControls.rotateSpeed; // which is adjusted by SetMouseSpeed message
  if (z > 0.0) {// use positive zoom values to avoid mirroring the stage
    stage.viewer.camera.zoom = z;
    stage.viewer.requestRender();
    ReturnClipPlaneDistances("CameraZoom");
  }
  let msg = "dx: " + dx.toString() + ", dy: " + dy.toString()
    + ", zoom:" + stage.viewer.camera.zoom.toString();
  if (isdebug)
    console.log(msg);
};


function AnimateRotation(axis, animatheta) {
  let m4 = new NGL.Matrix4();
  let then = 0;

  function render(now) {
// as in https://developer.mozilla.org/en-US/docs/Web/API/WebGL_API/Tutorial/Animating_objects_with_WebGL
    now *= 0.001;
    const deltaTime = now - then;
    then = now;

    if (animationspeed > 0.0)
      animatheta = (animatheta + deltaTime * animationspeed) % 360;
    if (shapeComp == null)
      return;

    RotateAxisComponents(axis, animatheta);
    stage.viewer.requestRender();

    if (animationspeed > 0.0)
      requestAnimationFrame(render);
  }
  if (animationspeed > 0.0)
    requestAnimationFrame(render);
}


function RotateAxisComponents(axis, theta)
{
  if (shapeComp == null)
    return;
  let m4 = new NGL.Matrix4();
  m4.makeRotationAxis(axis, theta);

  shapeComp.setTransform(m4);
  for (let i = 0; i < vectorshapeComps.length; i++) {
    if (typeof vectorshapeComps[i].reprList != "undefined")
      vectorshapeComps[i].setTransform(m4);
  }
}


async function RenderRequest(note = "")
{
  requestedby = note;
  if (requestedby != "")
    WebsockSendMsg(requestedby + '_BeforeRendering');
    // requestedby + '_AfterRendering' message is then sent by stage.viewer.signals.rendered
    // triggered by stage.viewer.requestRender
  await sleep(100).then(()=> { 
    stage.viewer.requestRender();
  });
  if (isdebug)
    WebsockSendMsg('RenderRequest ');
};


function getRotatedZoutMatrix()
{
  // Default in WebGL is for x-axis to point left and z-axis to point into the screen.
  // But we want x-axis pointing right and z-axis pointing out of the screen. 
  // Rotate coordinate system to that effect
  let m4 = new NGL.Matrix4();
  // GL matrices are the transpose of conventional rotation matrices
  m4.set(-1.0,  0.0,  0.0,  0.0,
          0.0,  1.0,  0.0,  0.0,
          0.0,  0.0, -1.0,  0.0,
          0.0,  0.0,  0.0,  1.0
  );
  return m4;
}


function SetDefaultOrientation(release_python_semaphore=true) {
  if (!rotationdisabled && animationspeed == 0.0)
    stage.viewerControls.orient(getRotatedZoutMatrix());
  SetAutoviewTimeout(shapeComp, 500, release_python_semaphore);
}


function SetAutoviewNoAnim(mycomponent)
{
// position component explicitly with SetAutoviewNoAnim() if autoview animation in 
// ResolveAutoview() is stalling which could happen on VMs in regression tests
  if (mycomponent == null || isAutoviewing==false)
    return;
  if (zoomanis != null)
  {
    mycomponent.stage.animationControls.pause();
    mycomponent.stage.animationControls.remove(zoomanis[0]);
    mycomponent.stage.animationControls.remove(zoomanis[1]);
    zoomanis = null;
  }
  let zaim = mycomponent.getZoom();
  let m = getRotatedZoutMatrix();
  m.multiplyScalar(-zaim);
  stage.viewerControls.orient(m);
  // ensure next call to ReturnClipPlaneDistances posts correct values to python
  // so python computes correct clipplane values
  stage.viewer.cDist = -stage.viewer.camera.position.z; 
  isAutoviewing = false;
  stage.viewer.requestRender();
  WebsockSendMsg('SetAutoView camera.z = ' + stage.viewer.camera.position.z.toString()); 
  WebsockSendMsg('FinishedSetAutoViewNoAnim forced (zaim= ' + zaim.toString() + ')'); // equivalent of the signal function
  return true;
};


async function ResolveAutoview(mycomponent, t)
{
  try {
    let zaim = mycomponent.getZoom();
    //zoomanis = mycomponent.stage.animationControls.zoomMove(mycomponent.getCenter(), zaim, t);
    let center = new NGL.Vector3();
    center.x=0;
    center.y=0;
    center.z=0;
    //if (camz===null || camz === -1)
    zaim = zaim*1.5;
    //else zaim = camz;
    // with HKL axes now part of vectorshape and not shape getCenter no longer returns (0,0,0)
    // and zaim often is too small. Work around this with explicit center at 0,0,0 and doubling zaim
    zoomanis = mycomponent.stage.animationControls.zoomMove(center, zaim, t);
    let dt = 50;
    let sumt = 0;
    while (isAutoviewing) 
    {
      // A workaround for lack of a signal function fired when autoView() has finished. autoView() runs 
      // asynchroneously in the background. Its completion time is at least t milliseconds and depends on the 
      // data size of mycomponent. It will have completed once the condition 
      // stage.viewer.camera.position.z == mycomponent.getZoom() is true. So fire our own signal 
      // at that point in time
      if (stage.viewer.camera.position.z == zaim && sumt > 0) 
      {
        if (isAutoviewing==true) 
        {
          let m = stage.viewerControls.getOrientation();
          let det = Math.pow(m.determinant(), 1/3);
          m.multiplyScalar(-zaim/det);
          stage.viewerControls.orient(m);
          WebsockSendMsg('FinishedSetAutoView');
          isAutoviewing = false;   
          return true;
        }
        else // set by SetAutoviewNoAnim()
          return;
      }
      await sleep(dt).then(()=> { 
        sumt += dt; 
        WebsockSendMsg('SetAutoView camera.z = ' + stage.viewer.camera.position.z.toString()); 
      } );   
    }
  }
  catch(err)
  {
    WebsockSendMsg('JavaScriptError: ' + err.stack );
  }
};


async function AutoViewPromiseRace(mycomponent, t)
{
  function onTimeoutResolveDefaultThreshold() 
  {// position component explicitly with SetAutoviewNoAnim() if autoview animation in 
   // ResolveAutoview() is stalling which could happen on VMs in regression tests
    return new Promise(async (resolve) => {
      setTimeout(() => {
        resolve(SetAutoviewNoAnim(mycomponent));
      }, t * 10) // start if 10 times t ms has lapsed
    });
  }
  return await Promise.race([ResolveAutoview(mycomponent, t),
                  onTimeoutResolveDefaultThreshold() ]);
}


async function SetAutoviewTimeout(mycomponent, t, release_python_semaphore)
{
  WebsockSendMsg('StartSetAutoView ');
  WebsockSendMsg('SetAutoView camera.z = ' + stage.viewer.camera.position.z.toString()); 
  isAutoviewing = true;

  let currentmode = stage.viewer.parameters.clipMode;
  let currentscale = stage.viewer.parameters.clipScale
  let currentnear = stage.viewer.parameters.clipNear;
  let currentfar = stage.viewer.parameters.clipFar;
  let currentcameraZ = stage.viewer.camera.position.z;

  AutoViewPromiseRace(mycomponent, t).then(()=>{
    requestedby = ""
    if (release_python_semaphore==true)
      requestedby = "AutoViewFinished"; // posts AutoViewFinished_AfterRendering in stage.viewer.signals.rendered.add()

    if (rotationdisabled || animationspeed > 0.0)
    {
      stage.viewer.parameters.clipMode = currentmode; 
      stage.viewer.parameters.clipScale = currentscale;
      stage.viewer.parameters.clipNear = currentnear; 
      stage.viewer.parameters.clipFar = currentfar;
      stage.viewer.camera.position.z = currentcameraZ;
    }
    stage.viewer.requestRender();
    ReturnClipPlaneDistances(); // updates zoom value in python */
  });
}


async function SendComponentRotationMatrixMsg() {
  if (shapeComp === null)
    return;
  await sleep(100);
  try {
    let msg = String(shapeComp.matrix.elements);
    WebsockSendMsg('CurrentComponentRotation:\n' + msg);
  }
  catch (err) {
    WebsockSendMsg('JavaScriptError: ' + err.stack);
  }
};

//async 
function SendOrientationMsg(funcname) {
  try {
    let msg = getOrientMsg();
    WebsockSendMsg(funcname + '_Orientation:\n' + msg);
  }
  catch (err) {
    WebsockSendMsg('JavaScriptError: ' + err.stack);
  }
};


// Log errors to debugger of your browser
function onError(e)
{
  let msg = 'WebSocket Error ' + e;
  console.log(msg);
  dbgmsg =msg;
};


function onOpen(e)
{
  let msg = 'Now connected via websocket to ' + pagename + '\n';
  WebsockSendMsg(msg);
  dbgmsg =msg;
  rerendered = false;
};


function onClose(e)
{
  let msg = 'Now disconnecting from websocket ' + pagename + '\n';
  console.log(msg);
  dbgmsg =msg;
};

var coordarray; // global for the binary data that are sent after receiving the AddHKLCoordinates message
var colourarray;// global for the binary data that are sent after receiving the AddHKLColours message
var radiiarray; // global for the binary data that are sent after receiving the AddHKLRadii message
var ttipids; // global for the binary data that are sent after receiving the AddHKLTTipIds message
var iswireframe = false;

function onMessage(e)
{
  let c,
    si;
  let val = null;
  let val2 = null;
  let datval = null;
  
  try
  {
    let msgtype = "";
    if (e.data instanceof ArrayBuffer == false) { // plain string data
      let showdata = e.data;
      if (showdata.length > 400)
        showdata = e.data.slice(0, 200) + '\n...\n' + e.data.slice(e.data.length - 200, -1);

      datval = e.data.split(":\n");
      msgtype = datval[0];
      if (datval.length == 1 && typeof msgtype === 'string') {// if message is empty we expect the next message to be an ArrayBuffer, i.e. a byte array sent from python
        binmsgtype = msgtype; // store the msgtype of the next message
        WebsockSendMsg('Waiting for ArrayBuffer for ' + binmsgtype);
        return;
      }
      else
        binmsgtype = "";
      val = datval[1].split(","); // assuming no commas in the received strings
      val2 = datval[1].split(";;"); // in case the received strings contain intended commas
    }

    if (e.data instanceof ArrayBuffer) { // binary data is received. 
      // Binary data is sent as a pair of messages by send_msg_to_browser(). First a plain string 
      // containing the message type is sent.
      // Then the actual data array is sent in binary format as a bytearray.
      // When HKLjavascripts notes the second message is in binary format it pairs it up with the 
      // previous message which is the message type so it can be processed
      msgtype = binmsgtype;
    }
    if (isdebug)
      WebsockSendMsg('Browser.JavaScript Got: ' + msgtype); // tell server what it sent us

    if (msgtype === "Reload")
    {
    // refresh browser with the javascript file
      if (stage != null)
      {
        let msg = getOrientMsg();
        WebsockSendMsg('OrientationBeforeReload:\n' + msg );
      }
      WebsockSendMsg( 'Refreshing ');

      sleep(200).then(()=> {
          socket_intentionally_closed = true;
          mysocket.close(4242, 'Refreshing ');
          ready_for_closing = true;
          window.location.reload(true);
        }
      );
    }

    if (stage == null) // everything below assumes stage!=null
      return;

    if (msgtype === "alpha")
    {
      let bin = parseInt(val[0]);
      if (bin < shapebufs.length)
      {
        alphas[bin] = parseFloat(val[1]);
        shapebufs[bin].setParameters({ opacity: alphas[bin] });
        if (expansion_shapebufs.length)
          for (let g=0; g < nrots; g++ )
            expansion_shapebufs[bin][g].setParameters({opacity: alphas[bin]});
        RenderRequest();
      }
    }

    if (msgtype === "colour")
    {
      let bin = parseInt(val[0]);
      if (bin < shapebufs.length)
      {
        let si =  parseInt(val[1]);
        colours[bin][3*si] = parseFloat(val[2]);
        colours[bin][3*si+1] = parseFloat(val[3]);
        colours[bin][3*si+2] = parseFloat(val[4]);
        shapebufs[bin].setAttributes({ color: colours[bin] });

        if (expansion_shapebufs.length)
          for (let g=0; g < nrots; g++ )
          {
            expansion_colours[bin][3*si] = parseFloat(val[2]);
            expansion_colours[bin][3*si+1] = parseFloat(val[3]);
            expansion_colours[bin][3*si+2] = parseFloat(val[4]);
            expansion_shapebufs[bin][g].setAttributes({ color: expansion_colours[bin] });
          }
        RenderRequest();
      }
    }

    if (msgtype === "DisplayTooltips")
    {
      let displaytooltips = val[0];
      stage.signals.hovered.removeAll();
      if (displaytooltips == "hover")
        stage.signals.hovered.add( HoverPickingProxyfunc );
    }

    if (msgtype === "ShowThisTooltip")
    {
      if (ttipid != "")
      {
        current_ttip = eval(datval[1]).split("\n\n")[0];
        current_ttip_ids = eval(datval[1]).split("\n\n")[1];

	      tooltip.innerText = current_ttip;
        tooltip.style.bottom = canvaspos.y + 7 + "px";
        tooltip.style.left = canvaspos.x + 8 + "px";
        tooltip.style.fontSize = fontsize.toString() + "pt";
        tooltip.style.display = "block";
      }
    }

    if (msgtype === "TooltipOpacity")
    {
      div_annotation_opacity = val[0];
      Object.assign(tooltip.style, {
        backgroundColor: "rgba(255, 255, 255, " + val[0] + " )",
      });
    }

    if (msgtype === "Redraw")
    {
      RenderRequest("notify_cctbx").then(()=> {
          SendOrientationMsg("Redraw");
          WebsockSendMsg( 'Redrawn ');
        }
      );
    }

    if (msgtype === "ReOrient")
    {
      WebsockSendMsg( 'Reorienting ');
      let sm = new Float32Array(16);
      for (let j=0; j<16; j++)
      {
        sm[j] = parseFloat(val[j]);
        if (isNaN( sm[j] ))
          return; // do nothing just in case
      }

      let m = new NGL.Matrix4();
      m.fromArray(sm);
      stage.viewerControls.orient(m);
      RenderRequest().then(()=> {
          SendOrientationMsg("ReOrient");
        }
      );
    }

    if (msgtype.includes("Expand") && shapeComp != null)
    {
      if (msgtype == "Expand" && expstate == "")
        return;

      if (msgtype == "ExpandP1" && expstate == "isP1Expanded")
        return;

      if (msgtype == "ExpandFriedel" && expstate == "isFriedelExpanded")
        return;

      if (msgtype == "ExpandP1Friedel" && expstate == "isP1FriedelExpanded")
        return;

      WebsockSendMsg('Expanding data...');
      // delete the shapebufs[] that holds the positions[] arrays
      shapeComp.removeRepresentation(repr);
      // remove shapecomp from stage first
      stage.removeComponent(shapeComp);

      expansion_positions = [];
      expansion_colours = [];
      expansion_radii = [];
      expansion_ttips = [];
      expansion_shapebufs = [];
      let nexpandrefls = 0;

      //alert('rotations:\n' + val);
      // Rotation matrices for the spacegroup come as a string of floats
      // separated by line breaks between each roation matrix
      let rotationstrs = datval[1].split("\n");
      let Rotmats = [];
      let r = new NGL.Vector3();

      for (let rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
      {
        Rotmats.push( new NGL.Matrix3() );
        // convert string of rotation matrix elements into a Matrix3
        let elmstrs = rotationstrs[rotmxidx].split(",");
        for (let j=0; j<9; j++)
          Rotmats[rotmxidx].elements[j] = parseFloat(elmstrs[j]);
      }

      let Imx = new NGL.Matrix3();
      Imx.identity(); // for testing
      if ( !(msgtype.includes("P1")) && rotationstrs.length == 1 && Rotmats[0].equals(Imx) )
        throw "Only the identity matrix is provided. That means no P1 expansion of reflections!";

      for (let bin=0; bin<nbins; bin++)
      {
        let nsize = positions[bin].length/3; // number of reflections in each bin
        let csize = nsize*3;
        let nsize3 = nsize*3;
        let anoexp = false;

        if (msgtype.includes("Friedel") )
        {
          anoexp = true;
          csize = nsize*6;
        }
        expansion_positions.push( [] );
        expansion_shapebufs.push( [] );
        expansion_colours.push( [] );
        expansion_radii.push( [] );
        expansion_ttips.push( [] );

        expansion_colours[bin] = colours[bin];
        expansion_radii[bin] = radii[bin];
        if (anoexp)
        {
          let colarr = [];
          let cl = colours[bin].length;
          for (let i=0; i<cl; i++)
          {
            colarr[i] = colours[bin][i];
            colarr[i+cl] = colours[bin][i];
          }
          expansion_colours[bin] = new Float32Array(colarr);

          let radiiarr = [];
          let rl = radii[bin].length;
          for (let i=0; i<rl; i++)
          {
            radiiarr[i] = radii[bin][i];
            radiiarr[i+rl] = radii[bin][i];
          }
          expansion_radii[bin] = new Float32Array(radiiarr);
        }

        nrots = 0;
        nexpandrefls = 0;
        for (let rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
        {
          if (rotationstrs[rotmxidx].length < 1 )
            continue;
          nrots++;
          expansion_positions[bin].push( [] );
          expansion_shapebufs[bin].push( [] );
          expansion_ttips[bin].push([]);
          if (anoexp) {
            expansion_ttips[bin][rotmxidx].ids = new Array(2 * nsize + 2);
            expansion_ttips[bin][rotmxidx].cartpos = new Array(2 * nsize);
          }
          else {
            expansion_ttips[bin][rotmxidx].ids = new Array(nsize + 1);
            expansion_ttips[bin][rotmxidx].cartpos = new Array(nsize);
          }
          expansion_ttips[bin][rotmxidx].ids[0] = rotmxidx; // id number of rotation. Used by PickingProxyfunc
          expansion_positions[bin][rotmxidx] = new Float32Array(csize);
          nexpandrefls += csize;

          for (let i=0; i<nsize; i++)
          {
            let idx= i*3;
            r.x = positions[bin][idx];
            r.y = positions[bin][idx+1];
            r.z = positions[bin][idx+2];
            r.applyMatrix3(Rotmats[rotmxidx]);
            expansion_positions[bin][rotmxidx][idx] = r.x;
            expansion_positions[bin][rotmxidx][idx + 1] = r.y;
            expansion_positions[bin][rotmxidx][idx + 2] = r.z;
            expansion_ttips[bin][rotmxidx].cartpos[i] = [r.x,r.y,r.z];
            expansion_ttips[bin][rotmxidx].ids[i + 1] = ttips[bin].ids[i + 1];
            if (anoexp)
            {
              r.negate(); // inversion for anomalous pair
              expansion_positions[bin][rotmxidx][nsize3 + idx] = r.x;
              expansion_positions[bin][rotmxidx][nsize3 + idx + 1] = r.y;
              expansion_positions[bin][rotmxidx][nsize3 + idx + 2] = r.z;
              expansion_ttips[bin][rotmxidx].cartpos[nsize + i] = [r.x, r.y, r.z];
 // indicate this hkl is a friedel mate with negative id of original hkl
              expansion_ttips[bin][rotmxidx].ids[nsize + 1 + i] = -ttips[bin].ids[i + 1];
            }
          }

          expansion_shapebufs[bin][rotmxidx] = new NGL.SphereBuffer({
              position: expansion_positions[bin][rotmxidx],
              color: expansion_colours[bin],
              radius: expansion_radii[bin],
// rotmxidx works as the id of the rotation of applied symmetry operator when creating tooltip for an hkl
              picking: expansion_ttips[bin][rotmxidx],
              }, { disableImpostor: iswireframe } );
          shape.addBuffer(expansion_shapebufs[bin][rotmxidx]);
          WebsockSendMsg('Expanded rotation operator ' + rotmxidx.toString());
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
      shapeComp = stage.addComponentFromObject(shape);
      MakeHKL_Axis();
      repr = shapeComp.addRepresentation('buffer', { wireframe: iswireframe });

      for (let bin=0; bin<nbins; bin++)
        for (let rotmxidx=0; rotmxidx < nrots; rotmxidx++ )
          expansion_shapebufs[bin][rotmxidx].setParameters( { opacity: alphas[bin] } ); 

      RotateAxisComponents(componentaxis, componenttheta); // apply any component rotation if specified on the GUI
      requestedby = "ExpandedInBrowser";
      stage.viewer.requestRender();
      //RenderRequest("ExpandedInBrowser"); //.then(()=> {
      WebsockSendMsg('Done ExpandedInBrowser ' + msgtype ); // "ExpandedInBrowser checked for in python
     // });
     }
    
    if (msgtype === "DisableZoomDrag") {
      WebsockSendMsg('using camerazoom');

      stage.mouseControls.remove("drag-shift-right");
      stage.mouseControls.remove("drag-shift-left");
	    stage.mouseControls.remove("drag-middle");
	    stage.mouseControls.add("drag-middle", CameraZoom);
      stage.mouseControls.add("drag-shift-right", CameraZoom);
      stage.mouseControls.add("drag-shift-left", CameraZoom);
      stage.mouseControls.add("scroll", CameraZoom);
      rotationdisabled = true;
    }

    if (msgtype === "EnableZoomDrag") {
      WebsockSendMsg('using mouse zoomdrag ');

      stage.mouseControls.remove("drag-shift-right");
      stage.mouseControls.remove("drag-shift-left");
	    stage.mouseControls.remove("drag-middle");
	    stage.mouseControls.add("drag-middle", NGL.MouseActions.zoomDrag);
      stage.mouseControls.add("drag-shift-right", NGL.MouseActions.zoomDrag);
      stage.mouseControls.add("drag-shift-left", NGL.MouseActions.zoomDrag);
      rotationdisabled = false;
    }

    if (msgtype === "DisableMouseRotation") {
      WebsockSendMsg('Fix mouse rotation');

      stage.mouseControls.remove("drag-left");
      stage.mouseControls.remove("scroll");
      stage.mouseControls.remove("scroll-ctrl");
      stage.mouseControls.remove("scroll-shift");
      stage.mouseControls.remove("drag-shift-right");
      stage.mouseControls.remove("drag-shift-left");
      rotationdisabled = true;
    }

    if (msgtype === "EnableMouseRotation") {
      WebsockSendMsg('Can mouse rotate ');
      stage.mouseControls.add("drag-left", NGL.MouseActions.rotateDrag);
      stage.mouseControls.add("scroll-ctrl", NGL.MouseActions.scrollCtrl);
      stage.mouseControls.add("scroll-shift", NGL.MouseActions.scrollShift);
      stage.mouseControls.remove("drag-shift-right");
      stage.mouseControls.remove("drag-shift-left");
      rotationdisabled = false;
    }

    if (msgtype === "RotateStage")
    { // rotate stage and its components
      while (isAutoviewing)
      {
        sleep(100);
      }

      WebsockSendMsg('Rotating stage ');

      let sm = new Float32Array(9);
      let m4 = new NGL.Matrix4();

      for (let j = 0; j < 9; j++)
        sm[j] = parseFloat(val[j]);

      // GL matrices are the transpose of conventional rotation matrices
      m4.set(sm[0], sm[3], sm[6], 0.0,
        sm[1], sm[4], sm[7], 0.0,
        sm[2], sm[5], sm[8], 0.0,
        0.0, 0.0, 0.0, 1.0
      );
      stage.viewerControls.orient(m4);
      stage.viewer.camera.zoom = parseFloat(val[9]);
      if (val[10] == "verbose")
        postrotmxflag = true;
      RenderRequest().then(()=> {
          ReturnClipPlaneDistances("RotateStage");
          SendOrientationMsg("RotateStage");
          WebsockSendMsg('Done RotateStage');
        }
      );
    }

    if (msgtype === "RotateAxisStage")
    { // rotate stage and its components
      WebsockSendMsg('Rotating stage around axis');

      let sm = new Float32Array(9);
      let m4 = new NGL.Matrix4();
      let axis = new NGL.Vector3();
      let theta = parseFloat(val[3]);
      axis.x = parseFloat(val[0]);
      axis.y = parseFloat(val[1]);
      axis.z = parseFloat(val[2]);
      m4.makeRotationAxis(axis, theta);
      stage.viewerControls.applyMatrix(m4);
      if (val[4] == "verbose")
        postrotmxflag = true;
      RenderRequest().then(()=> {
          ReturnClipPlaneDistances("RotateAxisStage");
          SendOrientationMsg("RotateAxisStage");
        }
      );
    }

    if (msgtype === "RotateComponents" && shapeComp != null)
    {
      WebsockSendMsg('Rotating components ');

      let sm = new Float32Array(9);
      let m4 = new NGL.Matrix4();
      stm4 = stage.viewerControls.getOrientation().elements;

      for (let j = 0; j < 9; j++)
        sm[j] = parseFloat(val[j]);

      // GL matrices are the transpose of conventional rotation matrices
      m4.set(sm[0], sm[3], sm[6], stm4[3],
        sm[1], sm[4], sm[7], stm4[7],
        sm[2], sm[5], sm[8], stm4[11],
        stm4[12], stm4[13], stm4[14], stm4[15]
      );
      shapeComp.setTransform(m4);
      if (val[9] == "verbose")
        postrotmxflag = true;
      RenderRequest().then(()=> {
          SendComponentRotationMatrixMsg();
        }
      );      
    }

    if (msgtype === "RotateAxisComponents" && shapeComp != null)
    { // rotating components from "Rotate around Vector" GUI control. Stage remains still
      WebsockSendMsg('Rotating components around axis ');
      let sm = new Float32Array(9);
      let m4 = new NGL.Matrix4();
      componenttheta = parseFloat(val[3]);
      componentaxis.x = parseFloat(val[0]);
      componentaxis.y = parseFloat(val[1]);
      componentaxis.z = parseFloat(val[2]);
      RotateAxisComponents(componentaxis, componenttheta);

      if (val[4] == "verbose")
        postrotmxflag = true;
      RenderRequest().then(()=> {
          SendComponentRotationMatrixMsg();
        }
      );      
    }

    if (msgtype === "AnimateRotateAxisComponents" && shapeComp != null) {
      WebsockSendMsg('Animate rotating components around axis ');
      animationspeed = parseFloat(val[3])*0.05;
      animaaxis.x = parseFloat(val[0]);
      animaaxis.y = parseFloat(val[1]);
      animaaxis.z = parseFloat(val[2]);

      if (animationspeed == 0.0)
        StopAnimateBtn.style.display = "None";
      else
        StopAnimateBtn.style.display = "Block";

      AnimateRotation(animaaxis, animatheta);

      SendComponentRotationMatrixMsg();
    }

    if (msgtype === "TranslateHKLpoints" && shapeComp != null)
    {
      WebsockSendMsg( 'Translating HKLs ' );
      let strs = datval[1].split("\n");
      let sm = new Float32Array(3);
      let elmstrs = strs[0].split(",");
      for (let j=0; j<3; j++)
        sm[j] = parseFloat(elmstrs[j]);
      shapeComp.setPosition([ sm[0], sm[1], sm[2] ]);
      RenderRequest().then(()=> {
          SendOrientationMsg("TranslateHKLpoints");
        }
      );
    }

    if (msgtype === "DrawSphere") {
      let pos = new Float32Array(3);
      let rgb = new Float32Array(3);
      for (let j = 0; j < 3; j++) {
        pos[j] = parseFloat(val2[j]);
        rgb[j] = parseFloat(val2[j + 3]);
      }
      let radius = parseFloat(val2[6]);
      let iswireframe = parseInt(val2[8]);

      if (vectorshape == null)
      {
        if (iswireframe == 1)
          vectorshape = new NGL.Shape('vectorshape', { disableImpostor: true });
        else
          vectorshape = new NGL.Shape('vectorshape');
      }

      vectorshape.addSphere(pos, rgb, radius);
      // if reprname is supplied then make a representation named reprname
      // of this and all pending spheres stored in vectorshape and render them.
      // Otherwise just accummulate the new sphere
      let reprname = val2[7].trim();
      if (reprname != "") {
        DeletePrimitives(reprname); // delete any existing vectors with the same name
        vectorshapeComps.push(stage.addComponentFromObject(vectorshape));
        if (iswireframe == 1)
          vectorreprs.push(
            vectorshapeComps[vectorshapeComps.length - 1].addRepresentation('vecbuf',
              { name: reprname, wireframe: true })
          )
        else
          vectorreprs.push(
            vectorshapeComps[vectorshapeComps.length - 1].addRepresentation('vecbuf',
              { name: reprname })
          )
        vectorshapeComps[vectorshapeComps.length - 1].autoView(500) // half a second animation
        vectorshape = null;
        RenderRequest();
      }
    }

    if (msgtype === "DrawVector")
    {
      componenttheta = 0.0;
      componentaxis.x = 0.0;
      componentaxis.y = 0.0;
      componentaxis.z = 1.0;
      // unapply any component rotation if specified on the "rotate around vector" GUI control to avoid 
      // messing up the stage orientation with the component orientation when aligning stage with a vector
      RotateAxisComponents(componentaxis, componenttheta);

      let r1 = new Float32Array(3);
      let r2 = new Float32Array(3);
      let rgb = new Float32Array(3);
      for (let j=0; j<3; j++)
      {
        r1[j] = parseFloat(val2[j]);
        r2[j] = parseFloat(val2[j+3]);
        rgb[j]= parseFloat(val2[j+6]);
      }
      let labelpos = parseFloat(val2[12]);
      let label = val2[9].trim();
      let radius = parseFloat(val2[11])* vectorwidth;
      let reprname = val2[10].trim();
      let autozoom = val2[13];
      DrawVector(r1, r2, rgb, radius, label, labelpos, reprname, autozoom);
    }

    if (msgtype === "RemovePrimitives")
    {
      RemovePrimitives(datval[1]);
    }

    if (msgtype === "DefineHKL_Axes")
    {
      let strarrs = datval[1].split("\n\n");
      let hstart = eval(strarrs[0]);
      let hend = eval(strarrs[1]);
      let kstart = eval(strarrs[2]);
      let kend = eval(strarrs[3]);
      let lstart = eval(strarrs[4]);
      let lend = eval(strarrs[5]);
      let hlabelpos = eval(strarrs[6]);
      let klabelpos = eval(strarrs[7]);
      let llabelpos = eval(strarrs[8]);

      DefineHKL_Axes(hstart, hend, kstart, kend, 
                 lstart, lend, hlabelpos, klabelpos, llabelpos)
    }

    if (msgtype === "SetFontSize")
    {
      fontsize = parseFloat(val[0]);
      MakeColourChart();
      MakeButtons();
      MakePlusMinusButtons();
    }

    if (msgtype === "SetVectorWidth")
    {
      vectorwidth = parseFloat(val[0]);
      MakeHKL_Axis();
      RenderRequest();
    }

    if (msgtype === "GetReflectionsInFrustum")
    {
      RenderRequest();
      maxReflInFrustum = parseInt(val[0]); 
      GetReflectionsInFrustum();
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
      let near = parseFloat(val[0]);
      let far = parseFloat(val[1]);
      let origcameraZpos = parseFloat(val[2]);
      let zoom = parseFloat(val[3]);
      stage.viewer.parameters.clipMode = 'camera';
      // clipScale = 'absolute' means clip planes are using scene dimensions
      stage.viewer.parameters.clipScale = 'absolute';
      clipFixToCamPosZ = true;

      if (near >= far) { // default to no clipping if near >= far
        stage.viewer.parameters.clipMode = 'scene';
        // clipScale = 'relative' means clip planes are in percentage
        stage.viewer.parameters.clipScale = 'relative';
        clipFixToCamPosZ = false;
        near = 0;
        far = 100;
      }
      
      stage.viewer.parameters.clipNear = near;
      stage.viewer.parameters.clipFar = far;
      if (clipFixToCamPosZ === true) {
        stage.viewer.parameters.clipNear = near + (origcameraZpos - stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = far + (origcameraZpos - stage.viewer.camera.position.z);
      }

      if (Number.isNaN(zoom) == false)
        stage.viewer.camera.zoom = zoom;
// provide a string so async RenderRequest() to call GetReflectionsInFrustum() after rendering
      RenderRequest("getfrustum").then(()=> {
          WebsockSendMsg('ClipPlaneDistancesSet'); // releases clipplane_msg_sem semaphore in jsview_3d.py
          SendOrientationMsg("getfrustum");
          GetReflectionsInFrustum();
        }
      );
    }

    if (msgtype === "GetClipPlaneDistances")
      ReturnClipPlaneDistances("GetClipPlaneDistances");

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

    if (msgtype === "PrintInformation") {
      hklequationmsg = datval[1];
      let infofsize = fontsize + 1; // slightly bigger font of infobanner to make it more prominent
      let wp = getTextWidth(hklequationmsg, infofsize)+3;
      if (infobanner != null) {
        infobanner.remove(); // delete previous infobanner if any
        infobanner = null;
      }
      pmleft = 70 + wp;
      MakePlusMinusButtons();
      if (hklequationmsg == "")
        return;
      infobanner = addBottomDivBox(hklequationmsg, pmbottom+3, 65, wp + 2, 5 + infofsize, 
                                    "rgba(255, 255, 255, 1.0)", infofsize);
    }

    if (msgtype === "SetBrowserDebug") {
      isdebug = (val[0] === "true");
    }

    if (msgtype === "SimulateClick") {
      document.dispatchEvent(simulate_click_evt);
    }

    if (msgtype === "UseWireFrame") {
      iswireframe = (val[0] === "true");
    }

    if (msgtype === "BackgroundColour") {
      stage.setParameters( { backgroundColor: datval[1] } )
    }

    if (msgtype ==="RemoveStageObjects")
    {
      RemoveStageObjects();
    }

    if (msgtype === "AddHKLCoordinates") {
      coordarray = new Float32Array(e.data);
    }

    if (msgtype === "AddHKLColours") {
      colourarray = new Float32Array(e.data);
    }

    if (msgtype === "AddHKLRadii") {
      radiiarray = new Float32Array(e.data);
    }

    if (msgtype === "AddHKLTTipIds") {
      ttipids = Array.from(new Float32Array(e.data)); // convert to plain javascript array so we can concatenate additional elements
      // assuming the above arrays have been initialised create the shape buffer
      AddSpheresBin2ShapeBuffer(coordarray, colourarray, radiiarray, ttipids, iswireframe);
    }

    if (msgtype === "MakeColourChart")
    {
      let msg = datval[1].split("\n\n");
      millerlabel = msg[0];
      fomlabel = msg[1];
      colourgradvalarrays = eval(msg[2]);
      MakeColourChart();
      RenderRequest();
    }

    if (msgtype ==="RenderStageObjects")
    {
      if (shape != null) 
      {
        shapeComp = stage.addComponentFromObject(shape);
        MakeHKL_Axis();
        MakeXYZ_Axis();
        repr = shapeComp.addRepresentation('buffer', { wireframe: iswireframe } );
        WebsockSendMsg('RenderStageObjects'); // tell python it may release hkls_drawn_sem semaphore
      }
    }

    if (msgtype == "SetDefaultOrientation")
    {
      SetDefaultOrientation();
      WebsockSendMsg('DefaultOrientSet '); 
    }

    if (msgtype === "SetAutoView")
    {
      SetAutoview(shapeComp, 500);
      WebsockSendMsg('AutoViewSet ');
    }

    if (msgtype === "MakeImage") {
      filename = val[0];
      stage.viewer.makeImage({ // using NGL's builtin function for making an image blob. html div legends are stripped
        factor: 1,
        antialias: true,
        trim: false,
        transparent: false
      }).then(function (blob) {
        if (parseInt(val[1]) < 3) {
          // Using websocket_server in python2 which doesn't allow streaming large compressed data
          // So use NGL's download image function
          NGL.download(blob, filename);
        }
        else { // websockets in python3 which supports streaming large blobs
          WebsockSendMsg('Imageblob', false);
          WebsockSendMsg(blob);
        }

        WebsockSendMsg('ImageWritten ');
      });
    }

    if (msgtype === "MakeImage2") {
      filename = val[0];
      //CHROME ONLY
      // html2canvas retains div legends when creaing an image blob
      let oldbtnstyle = StopAnimateBtn.style.display;
      ResetViewBtn.style.display = "None"; // hide buttons and other GUL controls on this webpage
      StopAnimateBtn.style.display = "None";
      html2canvas(document.getElementById("viewport")).then(function (canvas) {
        //blob = canvas.toDataURL("image/jpeg", 0.9);
        if (canvas.toBlob) {
          canvas.toBlob(function (blob) {
            if (parseInt(val[1]) < 3) {
              // Using websocket_server in python2 which doesn't allow streaming large compressed data
              // So use NGL's download image function
              NGL.download(blob, filename);
            }
            else { // websockets in python3 which supports streaming large blobs
              WebsockSendMsg('Imageblob', false);
              WebsockSendMsg(blob);
            }
          }, 'image/png')
        }
      });
      ResetViewBtn.style.display = "Block";
      StopAnimateBtn.style.display = oldbtnstyle;
      WebsockSendMsg('ImageWritten ');
    }

    if (msgtype === "MakeBrowserDataColumnComboBox")
    {
      if (columnSelect != null)
        columnSelect.remove(); 

      let msg = datval[1].split("\n\n");
      let columnSelect = createElement("select", {
        onchange: function (e)
        {
          WebsockSendMsg('SelectedBrowserDataColumnComboBox: ' + e.target.value);
        },
      }, { top: "25px", right: "10px", width: "130px", position: "absolute" }, fsize = fontsize);

      for (i = 0; i < msg.length - 1; i++) // last element is index of currently selected item
      {
        labelval = msg[i].split("\n");
        columnSelect.add(createElement("option", { text: labelval[0], value: labelval[1] }, fsize = fontsize));
      }
      addElement(columnSelect);
      columnSelect.options[ parseInt(msg[msg.length - 1]) ].selected = "true"; // display the selected element

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
      */
    }

    if (isdebug)
      WebsockSendMsg('Browser.JavaScript Done: ' + msgtype); // tell server We are done
  }

  catch(err)
  {
    WebsockSendMsg('JavaScriptError: ' + err.stack );
  }

};



function timefunc() {
  let d = new Date();
  let now = d.getTime();
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
  let script=document.createElement('script');
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
  backgroundColor: "rgba(255, 255, 255, " + div_annotation_opacity + ")",
  color: "black",
  padding: "0.1em",
  fontFamily: "sans-serif"
});


function DrawVector(r1, r2, rgb, radius, label, labelpos, reprname, autozoom)
{
  if (vectorshape == null)
    vectorshape = new NGL.Shape('vectorshape');

  let rad = radius*0.05;
  vectorshape.addArrow(r1, r2, [rgb[0], rgb[1], rgb[2]], rad);
  if (label !== "")
  {
    let pos = new NGL.Vector3()
    let txtR = [
      r1[0] * (1.0 - labelpos) + r2[0] * labelpos,
      r1[1] * (1.0 - labelpos) + r2[1] * labelpos,
      r1[2] * (1.0 - labelpos) + r2[2] * labelpos
    ];
    pos.x = txtR[0];
    pos.y = txtR[1];
    pos.z = txtR[2];

    let r = parseInt(255*rgb[0]).toString();
    let g = parseInt(255*rgb[1]).toString();
    let b = parseInt(255*rgb[2]).toString();

    let elm = createDivElement(label, [r,g,b]);
    elm.style.color = "white";
    elm.style.fontSize = fontsize.toString() + "pt";
    elm.style.padding = "4px"

    annodivs.push([elm, pos]);  // store until we get a representation name
  }
  // if reprname is supplied with a vector then make a representation named reprname
  // of this and all pending vectors stored in vectorshape and render them.
  // Otherwise just accummulate the new vector
  if (reprname != "")
  {
    let wasremoved = DeletePrimitives(reprname); // delete any existing vectors with the same name
    let cmp = stage.addComponentFromObject(vectorshape);
    for (let i = 0; i < annodivs.length; i++)
    {
      let elm = annodivs[i][0];
      let pos = annodivs[i][1];
      cmp.addAnnotation(pos, elm);
    }
    vectorshapeComps.push(cmp);
    annodivs = [];
    vectorreprs.push(
      vectorshapeComps[vectorshapeComps.length-1].addRepresentation('vecbuf',
                                                                  { name: reprname} )
    );
    if (autozoom == "True" && !wasremoved) // don't animate if adding a vector that was just removed above
      vectorshapeComps[vectorshapeComps.length - 1].autoView(500) // half a second animation
      //SetAutoviewTimeout(vectorshapeComps[vectorshapeComps.length - 1], 500);
    vectorshape = null;
    RenderRequest();
  }

}


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

  Hlabelvec.x = hlabelpos[0];
  Hlabelvec.y = hlabelpos[1];
  Hlabelvec.z = hlabelpos[2];
  Klabelvec.x = klabelpos[0];
  Klabelvec.y = klabelpos[1];
  Klabelvec.z = klabelpos[2];
  Llabelvec.x = llabelpos[0];
  Llabelvec.y = llabelpos[1];
  Llabelvec.z = llabelpos[2];
};


function MakeHKL_Axis()
{
  if (Hstarstart == null || Hstarend == null || Kstarstart == null || Kstarend == null
    || Lstarstart == null || Lstarend == null)
    return;
  //blue-x
  DrawVector(Hstarstart, Hstarend , [ 0, 0, 1 ], vectorwidth, "h", 1.05, "Haxis", false);
  //green-y
  DrawVector(Kstarstart, Kstarend, [ 0, 1, 0 ], vectorwidth, "k", 1.05, "Kaxis", false);
  //red-z
  DrawVector( Lstarstart, Lstarend, [ 1, 0, 0 ], vectorwidth, "l", 1.05, "Laxis", false);
};


function getOrientMsg()
{
  let cvorientmx = null;
  cvorientmx = stage.viewerControls.getOrientation();
  if (cvorientmx == null || cvorientmx.determinant() == 0)
      return oldmsg; // don't return invalid matrix

  cvorient = cvorientmx.elements;
  for (let j=0; j<16; j++)
  {
    if (Number.isNaN( cvorient[j]) )
      return oldmsg; // don't return invalid matrix
  }
  let cameradist;
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
  let msg = String(cvorient);
  oldmsg = msg;
  return msg;
}

// Distinguish between click and hover mouse events.
function HoverPickingProxyfunc(pickingProxy) { PickingProxyfunc(pickingProxy, 'hover'); }
function ClickPickingProxyfunc(pickingProxy) { PickingProxyfunc(pickingProxy, 'click'); }

// listen to hover or click signal to show a tooltip at an hkl or to post hkl id for matching entry in 
// millerarraytable in GUI or for visualising symmetry mates of the hkl for a given rotation operator
function PickingProxyfunc(pickingProxy, eventstr) {
  // adapted from http://nglviewer.org/ngl/api/manual/interaction-controls.html#clicked
  if (pickingProxy
    && (Object.prototype.toString.call(pickingProxy.picker["ids"]) === '[object Array]')
    && displaytooltips) {
    canvaspos = pickingProxy.canvasPosition;
    let sym_id = -1;
    let hkl_id = -1;
    let ids = [];
    let is_friedel_mate = 0;
    if (pickingProxy.picker["ids"].length > 0) { // get stored id number of rotation applied to this hkl
      sym_id = pickingProxy.picker["ids"][0]; // id of rotation stored when expanding to P1
      ids = pickingProxy.picker["ids"].slice(1); // ids of reflection
      hkl_id = ids[pickingProxy.pid];
      if (hkl_id < 0)
        is_friedel_mate = 1;
    }
    // tell python the id of the hkl and id of the rotation operator
    rightnow = timefunc();
    if (rightnow - timenow > tdelay)
    { // only post every tdelay milli second as not to overwhelm python
      ttipid = String([hkl_id, sym_id]);
	    if (ttipid != "")
	    {
        // send this to python which will send back a tooltip text
        if (pickingProxy.mouse.buttons == 1 || eventstr == 'hover') // left click or hover for tooltips
          WebsockSendMsg(eventstr + '_tooltip_id: [' + ttipid + ']');
        if (pickingProxy.mouse.buttons == 2) // right click for matching hkls in table
          WebsockSendMsg('match_hkl_id: [' + ttipid + ']');
        if (isdebug)
          console.log("current_ttip_ids: " + String(current_ttip_ids) + ", ttipid: " + String(ttipid));
       }
	    
      timenow = timefunc();
    }
  }
  else {
    tooltip.style.display = "none";
    current_ttip = "";
    ttipid = "";
  }
};


function getTextWidth(text, fsize=8)
{
  // re-use canvas object for better performance
  let canvas = getTextWidth.canvas || (getTextWidth.canvas = document.createElement("canvas"));
  let context = canvas.getContext("2d");
  context.font = fsize.toString() + "pt sans-serif";
  let metrics = context.measureText(text);
  return metrics.width;
}


function MakeColourChart()
{
  /* colourgradvalarrays is a list of colour charts. If only one list then it's one colour chart.
  Otherwise it's usually a list of colour charts that constitute a gradient across colours,
  typically used for illustrating figure of merits attenuating phase values in map coefficients
  */
  if (millerlabel == null || colourgradvalarrays == null)
    return;

  let hfac = 60.0 / colourgradvalarrays[0].length;
  let ih = 3.0*hfac,
  topr = 25.0,
  topr2 = 0.0,
  lp = 2.0; // vertical margin between edge of white container and labels
  let ctop = 10;
  let cleft = 10;
  let maxnumberwidth = 0;
  for (let j = 0; j < colourgradvalarrays[0].length; j++)
  {
    let val = colourgradvalarrays[0][j][0];
    maxnumberwidth = Math.max( getTextWidth(val, fontsize), maxnumberwidth );
  }
  let wp = maxnumberwidth + 5,
  lp2 = lp + wp,
  gl = 2,
  wp2 = gl,
  fomlabelheight = 25;

  if (colourgradvalarrays.length === 1)
  {
    wp2 = 15;
    fomlabelheight = 0;
  }
  let wp3 = wp + colourgradvalarrays.length * wp2 + 2;
  let totalheight = ih * colourgradvalarrays[0].length + 35 + fomlabelheight;

  if (colourchart != null)
    colourchart.remove(); // delete previous colour chart if any
  colourchart = addDivBox(null, ctop, cleft, wp3, totalheight, "rgba(255, 255, 255, 1.0)");

  // make a white box on top of which boxes with transparent background are placed
  // containing the colour values at regular intervals as well as label legend of
  // the displayed miller array
  addDiv2Container(colourchart, null, topr2, lp, wp3, totalheight, 'rgba(255, 255, 255, 1.0)');

  // print label of the miller array used for colouring
  let lblwidth = getTextWidth(millerlabel, fontsize);
  addDiv2Container(colourchart, millerlabel, topr2, lp, lblwidth + 5, 20, 'rgba(255, 255, 255, 1.0)', fontsize);

  if (fomlabel != "" )
  {
    // print FOM label, 1, 0.5 and 0.0 values below colour chart
    let fomtop = topr2 + totalheight - 18;
    let fomlp = lp + wp;
    let fomwp = wp3;
    let fomtop2 = fomtop - 13;
    // print the 1 number
    addDiv2Container(colourchart, 1, fomtop2, fomlp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fontsize);
    // print the 0.5 number
    let leftp = fomlp + 0.48 * gl * colourgradvalarrays.length;
    addDiv2Container(colourchart, 0.5, fomtop2, leftp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fontsize);
    // print the FOM label
    addDiv2Container(colourchart, fomlabel, fomtop, fomlp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fontsize);
    // print the 0 number
    leftp = fomlp + 0.96 * gl * colourgradvalarrays.length;
    addDiv2Container(colourchart, 0, fomtop2, leftp, fomwp, 20, 'rgba(255, 255, 255, 0.0)', fontsize);
  }

  for (let j = 0; j < colourgradvalarrays[0].length; j++)
  {
    let val = colourgradvalarrays[0][j][0];
    let topv = j*ih + topr;
    let toptxt = topv - 5;
    // print value of miller array if present in colourgradvalarrays[0][j][0]
    addDiv2Container(colourchart,val, toptxt, lp, wp, ih, 'rgba(255, 255, 255, 0.0)', fontsize);
  }

  // if colourgradvalarrays is an array of arrays then draw each array next to the previous one
  for (let g = 0; g < colourgradvalarrays.length; g++)
  {
    let leftp = g*gl + lp + wp;
    for (let j = 0; j < colourgradvalarrays[g].length; j++)
    {
      let R = colourgradvalarrays[g][j][1];
      let G = colourgradvalarrays[g][j][2];
      let B = colourgradvalarrays[g][j][3];
      let rgbcol = 'rgba(' + R.toString() + ',' + G.toString() + ',' + B.toString() + ', 1.0)'
      let topv = j * ih + topr;
      // add an extra pixel to height to ensure no white lines accidentally emerge from roundoff errors
      addDiv2Container(colourchart, null, topv, leftp, wp2, ih + 1, rgbcol);
    }
  }

  colourchart.oncontextmenu = function (e) { // oncontextmenu captures right clicks
    e.preventDefault()
    //alert("in oncontextmenu")
    return false;
  };
  
  colourchart.onselectstart = function () { // don't select numbers or labels on chart when double clicking the coulour chart
    return false;
  }

  if (isHKLviewer == true)
    colourchart.style.cursor = "pointer";

  colourchart.onclick = function (e) {
    let sel = window.getSelection();
    sel.removeAllRanges(); // don't select numbers or labels on chart when double clicking the coulour chart
    if (isHKLviewer == true)
      WebsockSendMsg('onClick colour chart');
  };
}


function MakeXYZ_Axis() {
  // draw x and y arrows
  let linelength = 20;
  let linestart = 5;
  let vleft = 2;
  let hbottom = 2;
  let arrowhalfwidth = 2;
  let linewidth = 1;
  let arrowlength = 7;
  let labelfromarrow = linelength + linestart + arrowlength + 1;
  let labelwidth = getTextWidth("x", fontsize)

  if (XYZaxes != null)
    XYZaxes.remove(); // delete previous colour chart if any
  XYZaxes = addDivBox(null, 0, 10, 60, 60, "rgba(255, 255, 255, 0.0)");
  XYZaxes.style.top = ""; // offset from the bottom so leave top distance void
  //XYZaxes.style.bottom = "30px";
  XYZaxes.style.bottom = pmbottom.toString() + "px";

  let velm = createElement("div", { innerText: "" },
    {
      backgroundColor: "rgba(0, 255, 0, 1.0)", color: "rgba(0, 0, 0, 0.0)",
      bottom: linestart.toString() + "px",
      height: linelength.toString() + "px",
      left: (vleft + arrowhalfwidth - linewidth / 2.0).toString() + "px",
      width: linewidth.toString() + "px",
      position: "absolute"
    }, 10);
  addElement(velm);
  XYZaxes.append(velm);

  let uparrow = createElement("div", { innerText: "" }, {
    backgroundColor: "rgba(0, 0, 255, 0.0)", color: "rgba(0, 0, 0, 0.0)",
    bottom: (linelength + linestart).toString() + "px",
    left: vleft.toString() + "px",
    borderLeft: arrowhalfwidth.toString() + "px solid transparent",
    borderRight: arrowhalfwidth.toString() + "px solid transparent",
    borderBottom: arrowlength.toString() + "px solid rgba(0, 255, 0)",
    position: "absolute"
  }, 10);
  addElement(uparrow);
  XYZaxes.append(uparrow);

  let yelm = createElement("div", { innerText: "y" }, {
    backgroundColor: "rgba(0, 255, 0, " + div_annotation_opacity + ")", color: "rgba(255, 255, 255, 1.0)",
    left: (vleft - linewidth / 2.0).toString() + "px",
    bottom: labelfromarrow.toString() + "px",
    padding: "1px",
    position: "absolute"
  }, fontsize);
  addElement(yelm);
  XYZaxes.append(yelm);

  let helm = createElement("div", { innerText: "" }, {
    backgroundColor: "rgba(0, 0, 255, 1.0)", color: "rgba(0, 0, 0, 0.0)",
    left: linestart.toString() + "px",
    width: linelength.toString() + "px",
    bottom: (hbottom + arrowhalfwidth - linewidth / 2.0).toString() + "px",
    height: linewidth.toString() + "px",
    position: "absolute"
  }, 10);
  addElement(helm);
  XYZaxes.append(helm);

  let rightarrow = createElement("div", { innerText: "" }, {
    backgroundColor: "rgba(0, 0, 255, 0.0)", color: "rgba(0, 0, 0, 0.0)",
    left: (linelength + linestart).toString() + "px",
    bottom: hbottom.toString() + "px",
    borderTop: arrowhalfwidth.toString() + "px solid transparent",
    borderBottom: arrowhalfwidth.toString() + "px solid transparent",
    borderLeft: arrowlength.toString() + "px solid rgba(0, 0, 255)",
    position: "absolute"
  }, 10);
  addElement(rightarrow);
  XYZaxes.append(rightarrow);

  let xelm = createElement("div", { innerText: "x" }, {
    backgroundColor: "rgba(0, 0, 255, " + div_annotation_opacity + ")", color: "rgba(255, 255, 255, 1.0)",
    left: labelfromarrow.toString() + "px",
    bottom: (hbottom + arrowhalfwidth - linewidth / 2.0).toString() + "px",
    padding: "1px",
    position: "absolute"
  }, fontsize);
  addElement(xelm);
  XYZaxes.append(xelm);

  let arrowradius = 10;
  let zarrow = createElement("div", { innerText: "" }, {
    backgroundColor: "rgba(255,0 ,0, 1.0)", color: "rgba(0, 0, 0, 1.0)",
    left: (linelength / 2 + linestart).toString() + "px",
    bottom: (linelength / 2 + linestart).toString() + "px",
    height: arrowradius.toString() + "px",
    width: arrowradius.toString() + "px",
    borderRadius: "50%",
    position: "absolute"
  }, 10);
  addElement(zarrow);
  XYZaxes.append(zarrow);

  let zarrowtip = createElement("div", { innerText: "" }, {
    backgroundColor: "rgba(255,255,255, 1.0)", color: "rgba(0, 0, 0, 1.0)",
    left: (arrowradius / 2 + linelength / 2 + linestart - 1).toString() + "px",
    bottom: (arrowradius / 2 + linelength / 2 + linestart - 1).toString() + "px",
    height: "2px",
    width: "2px",
    borderRadius: "50%",
    position: "absolute"
  }, 10);
  addElement(zarrowtip);
  XYZaxes.append(zarrowtip);

  let zelm = createElement("div", { innerText: "z" }, {
    backgroundColor: "rgba(255, 0, 0, " + div_annotation_opacity + ")", color: "rgba(255, 255, 255, 1.0)",
    left: (arrowradius + linelength / 2 + linestart).toString() + "px",
    bottom: (arrowradius + linelength / 2 + linestart).toString() + "px",
    padding: "1px",
    position: "absolute"
  }, fontsize);
  addElement(zelm);
  XYZaxes.append(zelm);
}


function AddSpheresBin2ShapeBuffer(coordarray, colourarray, radiiarray, ttipids) 
{
  // Primary function for creating the reflections to be displayed. There will be 
  // as many SphereBuffer elements in shapebufs as there are reflection bins requested 
  // by the user (default is 1).
  // Tooltip ids is a list of numbers matching the array index of the radiiarray 
  let ttiplst = [-1].concat(ttipids);
  // Prepend this list with -1. This value will be reassigned with an id nummber of 
  // a rotation operator when expanding to P1. PickingProxyfunc() will send back to cctbx.python the 
  // id number of the rotation operator and number in ttiplst matching the reflection that was clicked.
  let posarray = new Array(radiiarray.length);
  for (let i = 0; i < posarray.length; i++)
    posarray[i] = coordarray.slice(3 * i, 3 * i + 3);
  // cartpos: posarray is used in our own frustum culling function as well as for debugging
  // since it is simpler to access than shapebufs[0].geometry.attributes.position.array[4 * 3 * i]
  ttips.push( { ids: ttiplst, cartpos: posarray,
       getPosition: function() { return { x:0, y:0 }; } // dummy function to avoid crash
  }  );
  positions.push( new Float32Array( coordarray ) );
  colours.push( new Float32Array( colourarray ) );
  radii.push( new Float32Array( radiiarray ) );
  let curridx = positions.length - 1;
  
  shapebufs.push( new NGL.SphereBuffer({
    position: positions[curridx], // 1dim array [x0, y0, z0, x1, y1, z1,...] for all reflections
    color: colours[curridx], // 1dim array [R0, G0, B0, R1, G1, B1,...] for all reflections
    radius: radii[curridx], // 1dim array [r0, r1, r3,...]for all reflections
    picking: ttips[curridx],
  }, { disableImpostor: iswireframe })
  );

  if (shape == null)
    shape = new NGL.Shape('shape');
  shape.addBuffer(shapebufs[curridx]);

  alphas.push(1.0);
  nbins = nbins + 1;
}


function webgl2CanvasPosition(x,y,z) 
{
  let r = new NGL.Vector3();
  r.x = x;
  r.y = y;
  r.z = z;
  r.add(stage.viewer.translationGroup.position);
  r.applyMatrix4(stage.viewer.rotationGroup.matrix);
  r.applyMatrix4(stage.viewer.camera.projectionMatrix);
  let canvasX = (1 - r.x)*0.5 * stage.viewer.width;
  let canvasY = (1 + r.y)*0.5 * stage.viewer.height;
  return [canvasX, canvasY];
}


function GetReflectionsInFrustumFromBuffer(buffer) {
  // For the simple case where clip planes are parallel with the screen as in clipFar, clipNear.
  // Use cartesian coordinates of reflections stored in shapebufs[0].picking.cartpos
  let hkls_infrustums = [];
  let rotid_infrustum = [];
  if (buffer.parameters.opacity < 0.3) // use the same threshold as when tooltips won't show
    return [hkls_infrustums, rotid_infrustum];

  let nrefl = 0;
  for (let i = 0; i < buffer.picking.cartpos.length; i++)
  {
    let radius = 0.05* Math.abs(stage.viewer.parameters.clipFar -stage.viewer.parameters.clipNear)
    // possibly needs a better default value in case of using wireframe
    if (iswireframe == false) // no radius array available if using wireframe
      radius = buffer.geometry.attributes.radius.array[i];
    //radius=0.0;
    let x = buffer.picking.cartpos[i][0];
    let y = buffer.picking.cartpos[i][1];
    let z = buffer.picking.cartpos[i][2];

    let hklpos = new NGL.Vector3(x, y, z);
    let m = stage.viewer.modelGroup.children[0].matrixWorld;
    let currenthklpos = hklpos.applyMatrix4(m);
    let childZ = currenthklpos.z - stage.viewer.camera.position.z;
    // do rough exclusion of reflections from frustum with clipplanes for the sake of speed
    if ((childZ - radius) <stage.viewer.parameters.clipFar && (childZ + radius) > stage.viewer.parameters.clipNear)
    {
      if (nrefl >= maxReflInFrustum) // limit to avoid calling stage.viewer.pick() excessively
        return [hkls_infrustums, rotid_infrustum];
      let cv = webgl2CanvasPosition(x,y,z); 
      // stage.viewer.pick() calling readRenderTargetPixels() evaluates points in frustum
      let ret = stage.viewer.pick(cv[0], cv[1]); 
      if (ret.pid !== 0)
      {
        nrefl++;
        let hklid = buffer.picking.ids[i + 1];
        let rotid = buffer.picking.ids[0];
        hkls_infrustums.push(hklid);
        rotid_infrustum.push(rotid);
      }
    }
  }
  return [hkls_infrustums, rotid_infrustum];
}


function GetReflectionsInFrustum() {
  if (stage.viewer.parameters.clipScale != 'absolute')
    return; // we rely on clipNear clipFar being absolute values and not percentage
  let hkls_infrustums = [];
  let rotid_infrustum = [];
  for (let bin = 0; bin < shapebufs.length; bin++) {
    let arr = GetReflectionsInFrustumFromBuffer(shapebufs[bin]);
    hkls_infrustums = hkls_infrustums.concat(arr[0]);
    rotid_infrustum = rotid_infrustum.concat(arr[1]);
    if (expansion_shapebufs.length > 0)
      for (let rotmxidx = 0; rotmxidx < expansion_shapebufs[bin].length; rotmxidx++) {
        arr = GetReflectionsInFrustumFromBuffer(expansion_shapebufs[bin][rotmxidx]);
        hkls_infrustums = hkls_infrustums.concat(arr[0]);
        rotid_infrustum = rotid_infrustum.concat(arr[1]);
      }
  }
  WebsockSendMsg('InFrustum:[' + hkls_infrustums + ']:[' + rotid_infrustum + ']');
}


function MakePlusMinusButtons() {
  if (MinusBtn != null || PlusBtn != null) {
    stage.viewer.container.removeChild(MinusBtn);
    stage.viewer.container.removeChild(PlusBtn);
    PlusBtn.remove();
    MinusBtn.remove();
    MinusBtn = null;
    PlusBtn = null;
  }

  if (hklequationmsg == "")
    return;

  let btnwidth = getTextWidth("+", fontsize);
  PlusBtn = createElement("input", {
    value: "+",
    type: "button",
    onclick: function () {
      WebsockSendMsg('MoveClipPlanesUp');
    },
  },
    {
    bottom: pmbottom.toString() + "px", left: pmleft.toString() + "px",
        width: btnwidth.toString + "px", position: "absolute", fontWeight: "700" // bold font
    }, fontsize);
  addElement(PlusBtn);

  let left2 = btnwidth + pmleft + 20;
  btnwidth = getTextWidth("-", fontsize);
  MinusBtn = createElement("input", {
    value: "-",
    type: "button",
    onclick: function () {
      WebsockSendMsg('MoveClipPlanesDown');
    },
  }, {
      bottom: pmbottom.toString() + "px", left: left2.toString() + "px",
      width: btnwidth.toString + "px", position: "absolute", fontWeight: "700" // bold font
  }, fontsize);
  addElement(MinusBtn);
};


function MakeButtons() {
  if (ResetViewBtn != null)
    ResetViewBtn.remove();

  let btnwidth = getTextWidth("Reset View", fontsize);
  pmbottom = 5 + 3*fontsize; // used by MakePlusMinusButtons() 
  ResetViewBtn = createElement("input", {
    value: "Reset view",
    type: "button",
    onclick: function () {
      SetDefaultOrientation(false); // don't release python semaphore when SetDefaultOrientation isn't' called from python
      RenderRequest().then(()=> {
          SendOrientationMsg("MakeButtons");
        }
      );
    },
  }, { bottom: "10px", left: "10px", width: btnwidth.toString + "px", position: "absolute" }, fontsize);
  addElement(ResetViewBtn);

  if (StopAnimateBtn != null)
    StopAnimateBtn.remove();
  let leftoffset = 30 + btnwidth;
  btnwidth = getTextWidth("Toggle Animation", fontsize);
  StopAnimateBtn = createElement("input", {
    value: "Toggle Animation",
    type: "button",
    onclick: function () {
      animationspeed = -animationspeed;
      if (animationspeed > 0.0 && animaaxis.length() > 0.0)
        AnimateRotation(animaaxis, animatheta);
      WebsockSendMsg('ToggleAnimation');
    },
  }, { bottom: "10px", left: leftoffset.toString() + "px", width: btnwidth.toString + "px", position: "absolute" }, fontsize);
  addElement(StopAnimateBtn);
  StopAnimateBtn.style.display = "None";
}


function HKLscene()
{
  stage = new NGL.Stage('viewport', {  backgroundColor: "rgb(128, 128, 128)",
                                    tooltip:false, // create our own tooltip from a div element
                                    fogNear: 100, fogFar: 100 });
  stage.setParameters( { cameraType: camtype } );
// create tooltip element and add to the viewer canvas
  stage.viewer.container.appendChild(tooltip);
  // Always listen to click event as to display any symmetry hkls
  stage.signals.clicked.add(ClickPickingProxyfunc);
  
  stage.mouseObserver.signals.dragged.add(
    function ( deltaX, deltaY)
    {
      //let msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        postrotmxflag = true;
        ReturnClipPlaneDistances("mouseObserver.signals.dragged");
        SendOrientationMsg("MouseDragged");
        timenow = timefunc();
      }
      tooltip.style.display = "none";
    }
  );

  stage.mouseObserver.signals.clicked.add(
    function (x, y)
    {
      SendOrientationMsg("MouseClicked");
    }
  );

  stage.mouseObserver.signals.scrolled.add(
    function (delta)
    {
      //let msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        postrotmxflag = true;
        SendOrientationMsg("MouseScrolled");
        ReturnClipPlaneDistances("mouseObserver.signals.scrolled");
        timenow = timefunc();
      }
      tooltip.style.display = "none";
    }
  );

  stage.viewer.signals.rendered.add(
    function ()
    {
      if (postrotmxflag === true) 
      {
        SendOrientationMsg("ViewerRendered");
        postrotmxflag = false;
      }
      if (requestedby != "")
      {
        WebsockSendMsg(requestedby + '_AfterRendering');
        requestedby = "";
      }
    }
    
  );

  stage.viewerControls.signals.changed.add(
    function()
    {
      rightnow = timefunc();
      let t = 500;
      if (rightnow - timenow > t)
      { // only post every 250 milli second as not to overwhelm python
        //ReturnClipPlaneDistances();
        sleep(t).then(()=> {
            if (isAutoviewing)
              return;
            //SendOrientationMsg("CurrentView");
            //ReturnClipPlaneDistances("viewerControls.signals.changed");
          }
        );
        timenow = timefunc();
      }
    }
  );
 

  if (isdebug)
    stage.viewer.container.appendChild(debugmessage);

  // avoid NGL zoomFocus messing up clipplanes positions. So reassign those signals to zoomDrag
  stage.mouseControls.remove("drag-shift-right");
  stage.mouseControls.add("drag-shift-right", NGL.MouseActions.zoomDrag);
  stage.mouseControls.remove("drag-middle");
  stage.mouseControls.add("drag-middle", NGL.MouseActions.zoomDrag);
  stage.mouseControls.remove('clickPick-left'); // avoid undefined move-pick when clicking on a sphere
 
  //stage.mouseControls.remove("scroll");
  //stage.mouseControls.remove("scroll-ctrl");
  //stage.mouseControls.remove("scroll-shift");


  stage.viewer.requestRender();
  if (isdebug)
    debugmessage.innerText = dbgmsg;

  MakeButtons();

}


function OnUpdateOrientation()
{
  SendOrientationMsg('MouseMoved');
}


function PageLoad()
{
  try
  {
    let ret = MyWebGL_Detect(); // defined in webgl_check.js
    let msg = String(ret[1]);
    if (ret[0] == false) {
      let errmsg = "Critical WebGL problem! " + msg; 
      WebsockSendMsg(errmsg);
      throw new Error(errmsg, {cause: msg});
    }
    WebsockSendMsg(msg);
    //alert('In PageLoad');
    document.addEventListener('DOMContentLoaded', function () { 
      HKLscene(); 
    }, false );
    document.addEventListener('mouseup', function () {
      OnUpdateOrientation();
      GetReflectionsInFrustum();
    }, false);
    document.addEventListener('wheel', function (e) { OnUpdateOrientation(); }, false );
    document.addEventListener('scroll', function (e) { OnUpdateOrientation(); }, false );
    // mitigate flickering on some PCs when resizing
    document.addEventListener('resize', function () { RenderRequest(); }, false);

    document.addEventListener('simulate_click', function () { stage.viewer.requestRender(); }, false);
  }
  catch(err)
  {
    WebsockSendMsg('JavaScriptError: ' + err.stack );
  }
}


PageLoad();
