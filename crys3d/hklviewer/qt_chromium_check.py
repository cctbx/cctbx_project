from __future__ import absolute_import, division, print_function

from PySide2.QtCore import QTimer
from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
import sys, os


WeblglChecklibpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "webgl_check.js")

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

flgs = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = flgs + " --single-process" # necessary for detecting webgl abilities on Mac

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
  webglcheckscript = ""
  with open(WeblglChecklibpath, "r") as f:
    webglcheckscript = f.read()
    assert webglcheckscript != ""

  #app1.aboutToQuit.connect(lambda :print("about to quit"))
  print("QTWEBENGINE_CHROMIUM_FLAGS = %s" %os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] )
  # give the browser time to instatiate and then close down gracefully
  QTimer.singleShot(10000, app1.quit ) # in case pageloadFinished() is never executed
  browser = QWebEngineView()
  webpage = MyQWebEnginePage(browser)
  browser.setPage(webpage)

  def pageloadFinished( ok):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    browser.page().runJavaScript(webglcheckscript)
    # give the script time to run and then close us down gracefully
    QTimer.singleShot(1000, app1.quit )

  webpage.loadFinished.connect(pageloadFinished)
  webpage.setHtml(htmlstr1)
  browser.hide() # show() or hide() is necessary for loading the html page
  app1.exec_()
