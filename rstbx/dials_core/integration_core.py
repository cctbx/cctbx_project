from __future__ import absolute_import, division, print_function
from six.moves import range
from cctbx.array_family import flex
from scitbx import matrix

def show_observations(obs,out=None, n_bins=12):
  if out==None:
    import sys
    out = sys.stdout
  from libtbx.str_utils import format_value

  if n_bins is not None: # if None, assume binner already set up
    obs.setup_binner(n_bins = n_bins)
  result = []
  counts_given = obs.binner().counts_given()
  counts_complete = obs.binner().counts_complete()
  for i_bin in obs.binner().range_used():
    sel_w = obs.binner().selection(i_bin)
    sel_fo_all = obs.select(sel_w)
    d_max_,d_min_ = sel_fo_all.d_max_min()
    d_range = obs.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=True)
    sel_data = obs.select(sel_w).data()
    sel_sig = obs.select(sel_w).sigmas()

    if(sel_data.size() > 0):
      bin = resolution_bin(
        i_bin        = i_bin,
        d_range      = d_range,
        mean_I       = flex.mean(sel_data),
        n_work       = sel_data.size(),
        mean_I_sigI  = flex.mean(sel_data/sel_sig),
        d_max_min    = (d_max_, d_min_),
        completeness = (counts_given[i_bin], counts_complete[i_bin]))
      result.append(bin)
  print("\n Bin  Resolution Range  Compl.         <I>     <I/sig(I)>", file=out)
  for bin in result:
    fmt = " %s %s    %s  %s"
    print(fmt%(
      format_value("%3d",   bin.i_bin),
      format_value("%-17s", bin.d_range),
      format_value("%8.1f", bin.mean_I),
      format_value("%8.1f", bin.mean_I_sigI),
      ), file=out)
  return result

class resolution_bin(object):
  def __init__(self,
               i_bin         = None,
               d_range       = None,
               completeness  = None,
               alpha_work    = None,
               beta_work     = None,
               mean_I        = None,
               mean_I_sigI   = None,
               target_work   = None,
               target_free   = None,
               n_work        = None,
               n_free        = None,
               mean_f_obs    = None,
               fom_work      = None,
               scale_k1_work = None,
               pher_work     = None,
               pher_free     = None,
               sigmaa        = None,
               d_max_min     = None):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

from rstbx.apps import simple_integration
class integration_core(simple_integration):
  def __init__(self): simple_integration.__init__(self)
  """Integration concept.  Focus on each predicted spot position S.
     I(S): the integration mask for S is constructed by superimposing the
           bodypixels of the 10 nearest spotfinder spots; thus getting
           the maximum envelope.
     B(S): the background mask for S is obtained by finding an equal number
           of pixels as are in I(S).  Choose the nearest pixels to I(S)
           that are not within an exclusion distance (guard=3 pixels).
     P(S): A positional correction for I(S)--adjust the position of the mask
           away from the predicted position, before integrating.  This is constructed
           by considering the 10 closest tuples of (spotfinder spot,prediction).
           P(S) is the average positional offset for this set of 10.
     O(S): The set of spots close enough to S that they must be included
           in calculating where B(S) can be sampled.
  """

### NKS TO do list
# 1) account for spots that lie partially on inactive areas yet interfere with
#    the background of good spots
# 2) account for dead pixels (==-2) on the Pilatus

  def integration_masks_as_xy_tuples(self):
    values = []
    for imsk in range(len(self.BSmasks)):
      smask_keys = self.get_ISmask(imsk)
      for ks in range(0,len(smask_keys),2):
        values.append((smask_keys[ks],smask_keys[ks+1]))
    return values

  def background_masks_as_xy_tuples(self):
    values = []
    for imsk in range(len(self.BSmasks)):
      bmask = self.BSmasks[imsk]
      for key in bmask.keys():
        values.append((key[0],key[1]))
    return values

  def user_callback(self,dc,wxpanel,wx):
    # arguments are a wx Device Context, an Xray Frame, and the wx Module itself
    # BLUE: predictions
    for ix,pred in enumerate(self.predicted):
        if self.BSmasks[ix].keys()==[]:continue
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pred[1]/self.pixel_size,
          pred[0]/self.pixel_size)
        dc.SetPen(wx.Pen('blue'))
        dc.SetBrush(wx.BLUE_BRUSH)
        dc.DrawCircle(x,y,1)

    for imsk in range(len(self.BSmasks)):
      smask_keys = self.get_ISmask(imsk)
      bmask = self.BSmasks[imsk]
      if len(bmask.keys())==0: continue

      # CYAN: integration mask
      for ks in range(0,len(smask_keys),2):
        x,y = wxpanel._img.image_coords_as_screen_coords(smask_keys[ks+1],
                                                         smask_keys[ks])
        dc.SetPen(wx.Pen('cyan'))
        dc.SetBrush(wx.CYAN_BRUSH)
        dc.DrawCircle(x,y,2)

      # YELLOW: background mask
      for key in bmask.keys():
        x,y = wxpanel._img.image_coords_as_screen_coords(key[1],key[0])
        dc.SetPen(wx.Pen('yellow'))
        dc.SetBrush(wx.CYAN_BRUSH)
        dc.DrawCircle(x,y,2)

    for spot in self.spotfinder.images[self.frame_numbers[self.image_number]]["inlier_spots"]:
      # RED: spotfinder spot pixels
      for pxl in spot.bodypixels:
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pxl.y,
          pxl.x)
        dc.SetPen(wx.Pen('red'))
        dc.SetBrush(wx.RED_BRUSH)
        dc.DrawCircle(x,y,1)

      # GREEN: spotfinder centers of mass
      x,y = wxpanel._img.image_coords_as_screen_coords(
        spot.ctr_mass_y(),
        spot.ctr_mass_x())
      dc.SetPen(wx.Pen('green'))
      dc.SetBrush(wx.GREEN_BRUSH)
      dc.DrawCircle(x,y,1)

  def user_callback1(self,dc,wxpanel,wx):
    x,y = wxpanel._img.image_coords_as_screen_coords(100,100)
    dc.SetPen(wx.Pen('green'))
    dc.SetBrush(wx.GREEN_BRUSH)
    dc.DrawCircle(x,y,10)

  def initialize_increments(self,image_number=0):
    #initialize a data structure that contains possible vectors
    # background_pixel - spot_center
    # consider a large box 4x as large as the presumptive mask.
    from scitbx.array_family import flex
    Incr = []
    Distsq = flex.double()
    self.sorted = [] # a generic list of points close in distance to a central point
    if self.mask_focus[image_number] == None: return
    for i in range(-self.mask_focus[image_number][0],1+self.mask_focus[image_number][0]):
      for j in range(-self.mask_focus[image_number][1],1+self.mask_focus[image_number][1]):
        Incr.append(matrix.col((i,j)))
        Distsq.append(i*i+j*j)
    order = flex.sort_permutation(Distsq)
    for i in range(len(order)):
      #print i,order[i],Distsq[order[i]],Incr[order[i]]
      self.sorted.append(Incr[order[i]])

