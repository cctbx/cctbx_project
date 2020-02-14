from __future__ import absolute_import, division, print_function

from PySide2.QtCore import QTimer

from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEngineScript, QWebEnginePage
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


htmlstr = """

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

webglsupport = True

class MyQWebEnginePage(QWebEnginePage):
  webglsupport = True
  def javaScriptConsoleMessage(self,level, message, lineNumber, sourceID):
    print(message)
    if message == 'WebGL error':
      webglsupport = False

def CheckWebGL():
  script = QWebEngineScript()
  script.setInjectionPoint(QWebEngineScript.DocumentCreation)
  script.setSourceCode(jsstr)
  script.setInjectionPoint(QWebEngineScript.Deferred)

  browser = QWebEngineView()
  webpage = MyQWebEnginePage()
  webpage.setHtml(htmlstr)
  browser.setPage(webpage)
  #browser.show()
  browser.page().scripts().insert(script)
  print('WebGL=' + str(webglsupport),)
  # avoid "Release of profile requested but WebEnginePage still not deleted. Expect troubles !"
  webpage.deleteLater() 

if (__name__ == "__main__") :
  app1 = QApplication(sys.argv)
  # give the browser time to instatiate and after 1 second close down gracefully
  QTimer.singleShot(1000, app1.quit ) 
  CheckWebGL()
  app1.exec_()
