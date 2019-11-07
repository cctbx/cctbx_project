from __future__ import absolute_import, division, print_function
from six.moves import range
from rstbx_ext import SpotClass
from scitbx.array_family import flex
from rstbx.indexing_api import outlier_detection

class Graph:

  def __init__(self,fileout):
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.lib.pagesizes import letter
    self.c = Canvas(fileout,pagesize=letter)

class OutlierPlotPDF:
  def __init__(self,filename):
    self.R = Graph(filename)

  def setTitle(self,title):
    self.title=title

  def labelit_screen_output(self,lines):
    # outlier detection output
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch

    margin = 0.25
    block_text = self.R.c.beginText(margin*inch,letter[1] - margin*inch)
    block_text.setFont('Courier',7)
    for line in lines.split("\n"):
      block_text.textLine(line.rstrip())
    self.R.c.drawText(block_text)
    self.R.c.showPage()

def main_go(index_engine,verbose=False,phil_set=None):
    # first round of minimization
    if phil_set.indexing.outlier_detection.pdf is not None:
      phil_set.__inject__("writer",OutlierPlotPDF(phil_set.indexing.outlier_detection.pdf))

    # do some 12G parameter refinement here

    if verbose: print("Before outlier rejection, triclinic rmsd %.3f on %d spots"%(
      index_engine.residual(), index_engine.count_GOOD()))

    # outlier detection
    od = outlier_detection.find_outliers(ai=index_engine,verbose=verbose,
                                         horizon_phil=phil_set)
    if phil_set.indexing.outlier_detection.pdf is not None:
      od.make_graphs(canvas=phil_set.writer.R.c,left_margin=0.5)

    # do some 12G parameter refinement here, using the od.cache_status flags

    if verbose: print("-After outlier rejection, triclinic rmsd %.3f on %d spots"%(
      index_engine.residual(), index_engine.count_GOOD()))

    # update outlier graphs
    od.update(phil_set,ai=index_engine,status_with_marked_outliers=od.cache_status)

    if phil_set.indexing.outlier_detection.pdf is not None:
      od.make_graphs(canvas=phil_set.writer.R.c,left_margin=4.5)
      phil_set.writer.R.c.showPage()
      phil_set.writer.R.c.save()

    # estimate unit cell error

    if phil_set.indexing.outlier_detection.switch==True:
      # placeholder implementation XXX return to this later
      raw_spot_input = flex.vec3_double()
      assert len(process_dictionary['indexing'])==len(index_engine.get_observed_spot_positions(False))
      aipos = index_engine.get_observed_spot_positions(False)
      for ij in range(len(process_dictionary['indexing'])):
        spot = process_dictionary["indexing"][ij]
        if index_engine.get_status(ij)==SpotClass.OUTLIER:
          raw_spot_input.append((spot[0],spot[1],spot[2]))
      process_dictionary["indexing"] = raw_spot_input
      if len(raw_spot_input) < phil_set.distl_minimum_number_spots_for_indexing:
        message = "The remaining %d spots are less than the minimum %d set for indexing."%(
          len(raw_spot_input) , phil_set.distl_minimum_number_spots_for_indexing)
        print("Raising exception",message)
        raise Exception(message)
      print("Reindexing on the %d outlying spots to hunt for a second lattice"%len(raw_spot_input))
      index_engine = AutoIndexOne(process_dictionary,opt_rawframes = opt_frames,
                                  horizon_phil=phil_set)
      # do some 3G parameter refinement here
      # do some 12G parameter refinement here
      print("Reindexed OK")
"""Migration process:
   3) rid of dependency on process_dictionary and opt_ choices (DONE)
   4) get rid of saga spot status.  Use return values exclusively. (DONE)
   5) then implement efficiency by pushing to C++ the raw_spor_positions_mm_to_recip_space_xyz
   save 12 out of 36 seconds.
"""
