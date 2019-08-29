from __future__ import absolute_import, division, print_function
import wx
import wx.glcanvas
from gltbx.gl import *
from gltbx.glu import *

class MyCanvasBase(wx.glcanvas.GLCanvas):
    def __init__(self, parent):
        wx.glcanvas.GLCanvas.__init__(self, parent, -1)
        self.init = False
        # initial mouse position
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)


    def OnEraseBackground(self, event):
        pass # Do nothing, to avoid flashing on MSW.


    def OnSize(self, event):
        size = self.GetClientSize()
        if self.GetContext():
          self.SetCurrent()
          glViewport(0, 0, size.width, size.height)
        event.Skip()

    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent()
        if not self.init:
          self.InitGL()
          self.init = True
        self.OnDraw()

    def OnMouseDown(self, evt):
        self.CaptureMouse()

    def OnMouseUp(self, evt):
        self.ReleaseMouse()

    def OnMouseMotion(self, evt):
        if evt.Dragging() and evt.LeftIsDown():
            self.lastx, self.lasty = self.x, self.y
            self.x, self.y = evt.GetPosition()
            self.Refresh(False)


class CubeCanvas(MyCanvasBase):
    def InitGL(self):
        #print "cube InitCL"
        # set viewing projection
        glMatrixMode(GL_PROJECTION)
        glOrtho(-9.5, 9.5, -9.5, 9.5, 1.0, 3.0)

        # position viewer
        glMatrixMode(GL_MODELVIEW)
        glTranslatef(0.0, 0.0, -2.0)

        # position object
        glRotatef(self.y, 1.0, 0.0, 0.0)
        glRotatef(self.x, 0.0, 1.0, 0.0)

        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)

        glShadeModel(GL_FLAT)
        glFlush()

    def OnDraw(self):
        #print "cube OnDraw"
        # clear color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT)
        glClear(GL_DEPTH_BUFFER_BIT)

        # draw six faces of a cube
        if (1):
          glBegin(GL_QUADS)
          glNormal3f( 0.0, 0.0, 1.0)
          glVertex3f( 0.5, 0.5, 0.5)
          glVertex3f(-0.5, 0.5, 0.5)
          glVertex3f(-0.5,-0.5, 0.5)
          glVertex3f( 0.5,-0.5, 0.5)

          glNormal3f( 0.0, 0.0,-1.0)
          glVertex3f(-0.5,-0.5,-0.5)
          glVertex3f(-0.5, 0.5,-0.5)
          glVertex3f( 0.5, 0.5,-0.5)
          glVertex3f( 0.5,-0.5,-0.5)

          glNormal3f( 0.0, 1.0, 0.0)
          glVertex3f( 0.5, 0.5, 0.5)
          glVertex3f( 0.5, 0.5,-0.5)
          glVertex3f(-0.5, 0.5,-0.5)
          glVertex3f(-0.5, 0.5, 0.5)

          glNormal3f( 0.0,-1.0, 0.0)
          glVertex3f(-0.5,-0.5,-0.5)
          glVertex3f( 0.5,-0.5,-0.5)
          glVertex3f( 0.5,-0.5, 0.5)
          glVertex3f(-0.5,-0.5, 0.5)

          glNormal3f( 1.0, 0.0, 0.0)
          glVertex3f( 0.5, 0.5, 0.5)
          glVertex3f( 0.5,-0.5, 0.5)
          glVertex3f( 0.5,-0.5,-0.5)
          glVertex3f( 0.5, 0.5,-0.5)

          glNormal3f(-1.0, 0.0, 0.0)
          glVertex3f(-0.5,-0.5,-0.5)
          glVertex3f(-0.5,-0.5, 0.5)
          glVertex3f(-0.5, 0.5, 0.5)
          glVertex3f(-0.5, 0.5,-0.5)
          glEnd()

        glRotatef((self.x - self.lastx)/1., 1.0, 0.0, 0.0)
        glRotatef((self.y - self.lasty)/1., 0.0, 1.0, 0.0)
        glFlush()

        self.SwapBuffers()


class RunDemoApp(wx.App):
    def __init__(self, name="wx.glCanvas + gltbx"):
        self.name = name
        wx.App.__init__(self, redirect=False)

    def OnInit(self):
        wx.Log_SetActiveTarget(wx.LogStderr())
        self.SetAssertMode(wx.PYAPP_ASSERT_DIALOG)
        frame = wx.Frame(
          None, -1, self.name,
          pos=(50,50),
          size=(640,480),
          style=wx.DEFAULT_FRAME_STYLE)
        frame.CreateStatusBar()
        menuBar = wx.MenuBar()
        menu = wx.Menu()
        item = menu.Append(-1, "E&xit\tAlt-X", "Exit demo")
        self.Bind(wx.EVT_MENU, self.OnExitApp, item)
        menuBar.Append(menu, "&File")
        frame.SetMenuBar(menuBar)
        frame.Show(True)
        frame.SetSize((640, 480))
        win = CubeCanvas(frame)
        win.SetFocus()
        self.SetTopWindow(frame)
        self.frame = frame
        return True

    def OnExitApp(self, evt):
        self.frame.Close(True)

if __name__ == '__main__':
  RunDemoApp().MainLoop()
