import cspad_tbx
import numpy
import math
from scitbx.array_family import flex

# alternate implementation of hitfinder, use the idea of running spotfinder
#   on the data from the innermost four sensors.  Once this is done, a hit is
#   defined as an image where there are >=16 spots whose peak values
#   exceed the defined threshold.

class distl_hitfinder(object):

  def distl_filter(self,cspad_img,distance,timestamp,wavelength):
    self.hitfinder_d["DATA"] = cspad_img
    self.hitfinder_d["DISTANCE"] = distance
    self.hitfinder_d["TIMESTAMP"] = timestamp
    self.hitfinder_d["WAVELENGTH"] = wavelength

    import xfel_targets
    args = ["indexing.data=dummy",
            "distl.detector_format_version=CXI 5.1",
            "distl.bins.verbose=False",
            self.asic_filter,
            ] + targets[self.m_xtal_target]

    from spotfinder.applications.xfel import cxi_phil
    horizons_phil = cxi_phil.cxi_versioned_extract(args)
    horizons_phil.indexing.data = self.hitfinder_d

    from xfel.cxi import display_spots
    display_spots.parameters.horizons_phil = horizons_phil

    from rstbx.new_horizons.index import pre_indexing_validation,pack_names
    pre_indexing_validation(horizons_phil)
    imagefile_arguments = pack_names(horizons_phil)

    from spotfinder.applications import signal_strength
    info = signal_strength.run_signal_strength_core(horizons_phil,imagefile_arguments)

    imgdata = info.Files.images[0].linearintdata

    active_data = self.get_active_data(info.Files.images[0],horizons_phil)

    peak_heights = flex.int( [
      imgdata[ spot.max_pxl_x(), spot.max_pxl_y() ]
      for spot in info.S.images[info.frames[0]]["spots_total"]
    ])

    outscale = 256
    corrected = peak_heights.as_double() * self.correction
    outvalue = outscale *(1.0-corrected)
    outvalue.set_selected(outvalue<0.0,0.)
    outvalue.set_selected(outvalue>=outscale,int(outscale)-1)
    outvalue = flex.int(outvalue.as_numpy_array().astype(numpy.int32))
    # essentially, select a peak if the peak's ADU value is > 2.5 * the 90-percentile pixel value

    #work = display_spots.wrapper_of_callback(info)
    #work.display_with_callback(horizons_phil.indexing.data)
    return peak_heights,outvalue

  def get_active_data(self,imgobj,phil):
    active_areas = imgobj.get_tile_manager(phil).effective_tiling_as_flex_int()
    data = imgobj.linearintdata

    active_data = flex.double()
    for tile in xrange(len(active_areas)//4):
      block = data.matrix_copy_block(
          i_row=active_areas[4*tile+0],i_column=active_areas[4*tile+1],
          n_rows=active_areas[4*tile+2]-active_areas[4*tile+0],
          n_columns=active_areas[4*tile+3]-active_areas[4*tile+1]).as_1d().as_double()
      active_data = active_data.concatenate(block)

    #print "The mean is ",flex.mean(active_data),"on %d pixels"%len(active_data)
    order = flex.sort_permutation(active_data)
    #print "The 90-percentile pixel is ",active_data[order[int(0.9*len(active_data))]]
    #print "The 99-percentile pixel is ",active_data[order[int(0.99*len(active_data))]]

    adjlevel = 0.4
    brightness = 1.0
    percentile90 = active_data[order[int(0.9*len(active_data))]]
    if percentile90 > 0.:
      self.correction = brightness * adjlevel / percentile90
    else: self.correction = 1.0
    return active_data

  def set_up_hitfinder(self):
    self.hitfinder_d = cspad_tbx.dpack(
      active_areas    = self.active_areas,
      beam_center_x   = cspad_tbx.pixel_size * self.beam_center[0],
      beam_center_y   = cspad_tbx.pixel_size * self.beam_center[1],
      data            = flex.int(flex.grid(1765,1765),0),
      xtal_target     = self.m_xtal_target)

    #figure out which asics are on the central four sensors
    B = self.active_areas
    assert len(B)%4 == 0
    asics = [(B[i],B[i+1],B[i+2],B[i+3]) for i in xrange(0,len(B),4)]

    center = self.beam_center
    distances = flex.double()
    for iasic in xrange(len(asics)):
      cenasic = ((asics[iasic][2] + asics[iasic][0])/2. ,
                 (asics[iasic][3] + asics[iasic][1])/2. )

      distances.append(math.hypot(cenasic[0]-center[0], cenasic[1]-center[1]))
    orders = flex.sort_permutation(distances)

    flags = flex.int(len(asics),0)
    #Use the central 8 asics (central 4 sensors)
    for i in xrange(8): flags[orders[i]]=1
    self.asic_filter = "distl.tile_flags="+",".join(["%1d"%b for b in flags])
