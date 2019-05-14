from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
page_origin = (20.,220.)
boxedge = 500.

class PointTransform:
  '''provide the necessary transformation to go from image pixel coordinates
     to coordinates on the printed page of the .pdf report'''

  def __init__(self,detector_edge):
    self.boxedge = boxedge
    self.page_origin = page_origin
    self.size1 = detector_edge
    self.size2 = detector_edge
    self.subwindow_origin=[0.,0.]
    self.subwindow_fraction=1.


  def toPage(self, image_pixel_xy):

    image_fractional_coords = ((1.-image_pixel_xy[0]/self.size1),
                               image_pixel_xy[1]/self.size2)
    image_subwindow_coords = ((image_fractional_coords[1]-self.subwindow_origin[1])/
                              self.subwindow_fraction,
                              (image_fractional_coords[0]-self.subwindow_origin[0])/
                              self.subwindow_fraction)
    if 0.<image_subwindow_coords[0]<1. and 0.<image_subwindow_coords[1]<1.:
      page_coords = (image_subwindow_coords[0]*self.boxedge + self.page_origin[0],
                     (1.-image_subwindow_coords[1])*self.boxedge + self.page_origin[1]
      )
      return page_coords
    return None


from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import cm,mm
class Graph:

  def __init__(self,fileout):
    self.c = Canvas(fileout,pagesize=letter)

  def title(self,text):
    print(text)
    lines = text.split('\n')
    self.c.setFont('Helvetica',12)
    self.c.drawString(2*cm,26*cm,lines[0])
    if len(lines)>1:
      self.c.drawString(2*cm,25.5*cm,lines[1])

  def setTransform(self,detector_edge):
    #given the raw image fractional coordinates of the subwindow origin
    self.T = PointTransform(detector_edge)

  def __del__(self):
    self.c.save()

class PDF:
  def __init__(self,filename):
    self.R = Graph(filename)

  def make_image_plots_detail(self,ray_sim):

    normal = ray_sim.sim.tracing_impacts

    self.R.setTransform(ray_sim.detector.raw.focus()[0])
    self.R.title(
    "%.3f bandpass + %.3f degrees mosaicity (full widths); perfect optics"%(
      ray_sim.sim.bandpass,
      ray_sim.sim.mosaicity)+
"\nEnergy %4.1f KeV;  Detector distance %6.1f mm;   Limiting resolution %6.2f Angstrom"%(
      (12.398/(ray_sim.camera.lambda0*1E10)),
      ray_sim.camera.distance*1000.,
      ray_sim.structure.limiting_resolution))
    data_array = 255-ray_sim.image
    import numpy
    try:
      import PIL.Image as Image
    except ImportError:
      import Image
    imageout = Image.frombuffer("L",data_array.focus(),
      data_array.as_numpy_array().astype(numpy.uint8).tostring(),
      "raw","L",0,1
      )
    self.R.c.drawInlineImage(imageout,x=2*cm,y=9*cm, width=15*cm, height=15*cm)
    self.R.c.showPage()
    return self


if __name__=="__main__":
   data_array = flex.double(flex.grid((768,768)),1.0)
   print(data_array.focus())
   data_array = flex.double(flex.grid((7,7)),255)
   for x in range(7):
     data_array[(3,x)] = 0.
     data_array[(x,3)] = 0.
   try:
     import PIL.Image as Image
   except ImportError:
     import Image
   import numpy
   args = ("L",0,1)
   imageout = Image.frombuffer("L",data_array.focus(),
     data_array.as_float().as_numpy_array().astype(numpy.uint8).tostring(),
     "raw","L",0,1)

   imageout.save("newfile.png","PNG")
