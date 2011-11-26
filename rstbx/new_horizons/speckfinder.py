import math
from scitbx.array_family import flex
from libtbx import adopt_init_args

# alternate implementation of spotfinder; use the idea of identifying
# bright pixels on the inner 32 asics. (intensities > 2.0 * 90th percentile)
# Then cluster the brights into spots and use for autoindexing.
# Filter out any spots with area > 20 pixels...since this is a PAD with LCLS beam.

class speckfinder:

  def get_active_data(self):
    data = self.imgobj.linearintdata
    indexing = []
    for asic in self.corners:
      block = data.matrix_copy_block(
          i_row=asic[0],i_column=asic[1],
          n_rows=asic[2]-asic[0],
          n_columns=asic[3]-asic[1])
      active_data = block.as_1d().as_double()

      order = flex.sort_permutation(active_data)
      if self.verbose:
        print "The mean is ",flex.mean(active_data),"on %d pixels"%len(active_data)
        print "The 90-percentile pixel is ",active_data[order[int(0.9*len(active_data))]]
        print "The 99-percentile pixel is ",active_data[order[int(0.99*len(active_data))]]

      percentile90 = active_data[order[int(0.9*len(active_data))]]
      maximas = flex.vec2_double()
      for idx in xrange(len(active_data)-1, int(0.9*len(active_data)), -1):
        if active_data[order[idx]] > 2.0 * percentile90:
          if self.verbose: print "    ", idx, active_data[order[idx]]
          irow = order[idx] // (asic[3]-asic[1])
          icol = order[idx] % (asic[3]-asic[1])
          #self.green.append((asic[0]+irow, asic[1]+icol))
          maximas.append((irow, icol))
      CLUS = clustering(maximas)
      #coords = CLUS.as_spot_max_pixels(block,asic)
      coords = CLUS.as_spot_center_of_mass(block,asic,percentile90)
      intensities = CLUS.intensities
      for coord,height in zip(coords,intensities):
        self.green.append(coord)
        indexing.append( (
          coord[0] * float(self.inputpd["pixel_size"]),
          coord[1] * float(self.inputpd["pixel_size"]),
          0.0, # 0 -degree offset for still image
          height)
        )
    return indexing

  def __init__(self,imgobj,phil,inputpd,verbose=False):
    adopt_init_args(self,locals())
    self.active_areas = imgobj.get_tile_manager(phil).effective_tiling_as_flex_int()
    B = self.active_areas

    #figure out which asics are on the central four sensors
    assert len(self.active_areas)%4 == 0
    # apply an additional margin of 1 pixel, since we don't seem to be
    # registering the global margin.
    asics = [(B[i]+1,B[i+1]+1,B[i+2]-1,B[i+3]-1) for i in xrange(0,len(B),4)]

    from scitbx.matrix import col
    centre_mm = col((float(inputpd["xbeam"]),float(inputpd["ybeam"])))
    centre = centre_mm / float(inputpd["pixel_size"])
    distances = flex.double()
    cenasics = flex.vec2_double()
    self.corners = []
    for iasic in xrange(len(asics)):
      cenasic = ((asics[iasic][2] + asics[iasic][0])/2. ,
                 (asics[iasic][3] + asics[iasic][1])/2. )
      cenasics.append(cenasic)
      distances.append(math.hypot(cenasic[0]-centre[0], cenasic[1]-centre[1]))
    orders = flex.sort_permutation(distances)

    self.flags = flex.int(len(asics),0)
    #Use the central 8 asics (central 4 sensors)
    self.green = []
    for i in xrange(32):
      #self.green.append( cenasics[orders[i]] )
      self.corners.append (asics[orders[i]])
      #self.green.append((self.corners[-1][0],self.corners[-1][1]))
      self.flags[orders[i]]=1
    self.asic_filter = "distl.tile_flags="+",".join(["%1d"%b for b in self.flags])

class clustering:
  def __init__(self,maxima):
    self.verbose=False
    #for the next spot, define the universe of possible pixels (working targets)
    # and the map of which pixels have been visited already (pixel_visited)
    working_targets = list(xrange(len(maxima)))
    pixel_visited = flex.bool(len(working_targets),False)
    self.spots = []
    while len(working_targets) > 0:
      # for this spot, indices of pixels known to be in the spot (pixel_members)
      # and a stack of indices of pixels still in process of testing for connections
      #  (connection_stack).  Pop/push operates on the end of stack
      pixel_members = [working_targets[0]]
      connection_stack = [0]
      assert len(pixel_visited)==len(working_targets)
      pixel_visited[0]=True
      while len(connection_stack) > 0:
        idx_current = connection_stack[-1]
        for idx_target in xrange(len(working_targets)):
          if not pixel_visited[idx_target]:
            distance = math.hypot( maxima[working_targets[idx_current]][0]-
                                   maxima[working_targets[idx_target]][0],
                                   maxima[working_targets[idx_current]][1]-
                                   maxima[working_targets[idx_target]][1])
            if distance >= 2.0: continue
            pixel_visited[idx_target]=True
            pixel_members.append(working_targets[idx_target])
            connection_stack.append(idx_target)
        if connection_stack[-1] == idx_current: connection_stack.pop()
      if self.verbose: print "new spot with %d pixels"%len(pixel_members),pixel_members
      for idx in pixel_members:#[working_targets[i] for i in pixel_members]:
        working_targets.remove(idx)
      pixel_visited = flex.bool(len(working_targets),False)
      self.spots.append( [maxima[i] for i in pixel_members] )

  def as_spot_max_pixels(self,active_block,asic):
    maxima = flex.vec2_double()
    for spot in self.spots:
      if self.verbose:print [(int(row)+asic[0],int(col)+asic[1]) for row,col in spot]
      pixel_values = [active_block[(int(row),int(col))] for row,col in spot]
      if self.verbose:print "PIXEL_VALUES",pixel_values
      addr = pixel_values.index(max(pixel_values))
      maxima.append( ( int(spot[addr][0]) + asic[0], int(spot[addr][1]) + asic[1] ))
      if self.verbose:print
    return maxima

  def as_spot_center_of_mass(self,active_block,asic,percentile90):
    from scitbx.matrix import col
    maxima = flex.vec2_double()
    self.intensities = flex.double()
    for spot in self.spots:
      pixels = [col(((frow)+asic[0],(fcol)+asic[1])) for frow,fcol in spot]
      pixel_values = [active_block[(int(frow),int(fcol))] for frow,fcol in spot]
      numerator = col((0.0,0.0)); denominator = 0.0
      for ispot in xrange(len(spot)):
        numerator += (pixel_values[ispot]-percentile90)*pixels[ispot]
        denominator += pixel_values[ispot]-percentile90
      if len(spot) < 20:
        #any spot with more than 20 pixels is a clear outlier
        maxima.append( numerator/denominator )
        self.intensities.append(denominator)
    return maxima
