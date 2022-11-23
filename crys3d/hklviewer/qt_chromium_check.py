from __future__ import absolute_import, division, print_function

from PySide2.QtCore import QTimer
from PySide2.QtWebEngineWidgets import QWebEngineView
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
          console.log('WebGL is enabled')
          if (return_context) {
            console.log(names[i])
            // return WebGL object if the function's argument is present
            return {name:names[i], gl:context};
          }
          return true;
        }
      } catch(e) {}
    }
    console.log('WebGL is supported, but disabled')
    return false;
  }
  console.log('WebGL not supported')
  return false;
}

retval = MyWebGL_Detect();
retmsg = "WebGL triage: " + String(retval);

if (MyWebGL_Detect() == false )
   console.log('WebGL error');
else
  console.log('WebGL works');

"""



htmlstr1 = """

<html lang="en">
<head>
  <title>Detect webgl</title>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
</head>
<body>This is a test for webgl <div id="viewport" style="width:100%; height:100%;"></div>
</body>
</html>
"""

# some documentation for chromium flags
# https://peter.sh/experiments/chromium-command-line-switches

os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] += " --single-process" # necessary for detecting webgl abilities on Mac

class MyQWebEnginePage(QWebEnginePage):
  def __init__(self, *args, **kwargs):
    QWebEnginePage.__init__(self, *args, **kwargs)
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
  from PySide2.QtWidgets import QApplication
  app1 = QApplication(sys.argv)
  #app1.aboutToQuit.connect(lambda :print("about to quit"))
  print("QTWEBENGINE_CHROMIUM_FLAGS = %s" %os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] )
  # give the browser time to instatiate and then close down gracefully
  QTimer.singleShot(10000, app1.quit ) # in case pageloadFinished() is never executed
  browser = QWebEngineView()
  webpage = MyQWebEnginePage(browser)
  browser.setPage(webpage)

  def pageloadFinished( ok):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    browser.page().runJavaScript(jsstr)
    # give the script time to run and then close us down gracefully
    QTimer.singleShot(1000, app1.quit )

  webpage.loadFinished.connect(pageloadFinished)
  webpage.setHtml(htmlstr1)
  browser.hide() # show() or hide() is necessary for loading the html page
  app1.exec_()

