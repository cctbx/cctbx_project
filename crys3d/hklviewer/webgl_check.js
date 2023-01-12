// JavaScript source code

function MyWebGL_Detect(return_context)
{
  if (!!window.WebGLRenderingContext) {
    var canvas = document.createElement("canvas"),
         names = ["webgl2", "webgl", "experimental-webgl", "moz-webgl", "webkit-3d"],
       context = false;
    for(var i=0;i< names.length;i++)
    {
      try
      {
        context = canvas.getContext(names[i]);
        if (context && typeof context.getParameter == "function") {
          console.log('WebGL is enabled')
          if (return_context) {
            console.log(names[i])
            // return WebGL object if the function's argument is present
            return [true, {name:names[i], gl:context}];
          }
          return [true, 'WebGL is enabled'];
        }
      } catch(e) {}
    }
    console.log('WebGL is supported, but disabled')
    return [false, 'WebGL is supported, but disabled'];
  }
  console.log('WebGL is not supported')
  return [false, 'WebGL is not supported'];
}

retval = MyWebGL_Detect();
retmsg = "WebGL triage: " + String(retval);

if (MyWebGL_Detect() == false )
  console.log('WebGL error');
else
  console.log('WebGL works');
