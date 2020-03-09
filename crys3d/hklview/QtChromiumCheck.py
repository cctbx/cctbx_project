from __future__ import absolute_import, division, print_function

from PySide2.QtCore import QTimer

from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PySide2.QtWidgets import QApplication
import sys, os


jsstr = """
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
          // WebGL is enabled
          if (return_context) {
            // return WebGL object if the function's argument is present
            return {name:names[i], gl:context};
          }
          // else, return just true
          return true;
        }
      } catch(e) {}
    }
    // WebGL is supported, but disabled
    return false;
  }
  // WebGL not supported
  return false;
}


retval = MyWebGL_Detect();
retmsg = "WebGL triage: " + String(retval);
//alert(retmsg);

if (MyWebGL_Detect() == false )
{
  //alert('WebGL error');
  console.log('WebGL error');
}
else
{
  //alert('WebGL works');
  console.log('WebGL works');
}

"""



htmlstr1 = """

<html lang="en">
<head>
  <title>Detect webgl</title>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
</head>

<body>
This is a test for webgl
<div id="viewport" style="width:100%; height:100%;"></div>


</body>
</html>

"""




os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--single-process" # necessary for detecting webgl abilities on Mac


class MyQWebEnginePage(QWebEnginePage):
  def __init__(self, *args, **kwargs):
    #print('in MyQWebEnginePage')
    QWebEnginePage.__init__(self, *args, **kwargs)
    self.webglsupport = None
  def javaScriptConsoleMessage(self,level, message, lineNumber, sourceID):
    print(message)
    self.webglsupport = True
    if message == 'WebGL error':
      self.webglsupport = False
    return super(MyQWebEnginePage,self).javaScriptConsoleMessage(level, message, lineNumber, sourceID)
  def javaScriptAlert(self,securityOrigin,msg):
    print(msg)
    return super(MyQWebEnginePage,self).javaScriptAlert(securityOrigin,msg)


if (__name__ == "__main__") :
  app1 = QApplication(sys.argv)
  # give the browser time to instatiate and then close down gracefully
  QTimer.singleShot(10000, app1.quit ) # in case pageloadFinished() is never executed
  browser = QWebEngineView()
  webpage = MyQWebEnginePage(browser)
  browser.setPage(webpage)

  def pageloadFinished( ok):
    #print('webpage load has finished: ' + str(ok))
    browser.page().runJavaScript(jsstr)
    # give the script time to run and then close us down gracefully
    QTimer.singleShot(1000, app1.quit )

  webpage.loadFinished.connect(pageloadFinished)
  webpage.setHtml(htmlstr1)
  browser.hide() # show() or hide() is necessary for loading the html page
  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  #print('WebGL=' + str(webpage.webglsupport),)

  app1.exec_()
