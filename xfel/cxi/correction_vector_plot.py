import os,math
from scitbx import matrix
from cctbx.array_family import flex
from xfel import correction_vector_store,get_correction_vector_xy

class lines(correction_vector_store):
  def __init__(self,params):
    correction_vector_store.__init__(self)
    self.params = params

  def literals(self):
    if self.params.run_numbers is None:
      path = self.params.outdir_template
      stream = open(path,"r")
      print path
      reading_setting = False
      for line in stream.readlines():
        if line.find("Cell in setting %d"%self.params.bravais_setting_id)>=0:
          reading_setting = True
        if line.find("-----END")>=0:
          reading_setting = False
        if reading_setting and line.find("CV OBSCENTER")==0:
          yield line.strip(),0
        if reading_setting and len(self.tiles)==0 and line.find("EFFEC")==0:
          self.tiles = flex.int([int(a) for a in line.strip().split()[2:]])
          assert len(self.tiles)==256
          print list(self.tiles)
          self.initialize_per_tile_sums()
      return
    for run in self.params.run_numbers:
        templ = self.params.outdir_template%run
        items = os.listdir(templ)
        for item in items:
          path = os.path.join(templ,item)
          stream = open(path,"r")
          print path
          reading_setting = False
          for line in stream.readlines():
            if line.find("Cell in setting %d"%self.params.bravais_setting_id)>=0:
              reading_setting = True
            if line.find("-----END")>=0:
              reading_setting = False
            if reading_setting and line.find("CV OBSCENTER")==0:
              yield line.strip(),run
            if reading_setting and self.tiles is None and line.find("EFFEC")==0:
              self.tiles = flex.int([int(a) for a in line.strip().split()[2:]])
              assert len(self.tiles)==256
              print list(self.tiles)
              self.initialize_per_tile_sums()

  def vectors(self):
    self.tile_rmsd = [0.]*64
    for line,run in self.literals():
     try:
      tokens = line.split()
      if len(tokens)!=13:
        continue # guard against subprocess contention or unflushed output
      self.register_line( float(tokens[2]),float(tokens[3]),
                       float(tokens[5]),float(tokens[6]),
                       float(tokens[8]),float(tokens[9]),
                       float(tokens[11]),float(tokens[12]) )

      yield "OK"
     except ValueError:
       print "Valueerror"
    for x in xrange(64):
      if self.tilecounts[x]==0: continue
      self.radii[x]/=self.tilecounts[x]
      sum_cv = matrix.col(self.mean_cv[x])
      self.mean_cv[x] = sum_cv/self.tilecounts[x]
      mean_cv = matrix.col(self.mean_cv[x])
      selection = (self.master_tiles == x)
      selected_cv = self.master_cv.select(selection)
      if len(selected_cv)>0:
        self.tile_rmsd[x] = math.sqrt(
      flex.mean(flex.double([ (matrix.col(cv) - mean_cv).length_sq() for cv in selected_cv ]))
      )
      else: self.tile_rmsd[x]=0.
    self.overall_N = flex.sum(flex.int( [int(t) for t in self.tilecounts] ))
    self.overall_cv = matrix.col(self.overall_cv)/self.overall_N
    self.overall_rmsd = math.sqrt( self.sum_sq_cv / self.overall_N )

def run_correction_vector_plot(working_phil):

  L = lines(working_phil)
  for line in L.vectors():
    pass # pull out the information, lines class does all the work

  close_x = flex.double()
  close_y = flex.double()
  far_x = flex.double()
  far_y = flex.double()
  master_coords = L.master_coords
  master_cv = L.master_cv
  master_tiles = L.master_tiles
  for idx in xrange(0,len(master_coords),10):
    if matrix.col(master_cv[idx]).length() < L.tile_rmsd[ master_tiles[idx] ]:
      pass
      #close_x.append(master_coords[idx][0])
      #close_y.append(master_coords[idx][1])
    else:
      far_x.append(master_coords[idx][0])
      far_y.append(master_coords[idx][1])
      close_x.append(master_coords[idx][0]+master_cv[idx][0])
      close_y.append(master_coords[idx][1]+master_cv[idx][1])
  if working_phil.show_plots is True:
    from matplotlib import pyplot as plt
    plt.plot(close_x,close_y,"r.")
    plt.plot(far_x,far_y,"g.")
    plt.axes().set_aspect("equal")
    plt.show()


  sort_radii = flex.sort_permutation(flex.double(L.radii))
  tile_rmsds = flex.double()
  for idx in xrange(64):
    x = sort_radii[idx]
    print "Tile %2d: radius %7.2f, %6d observations, delx %5.2f  dely %5.2f, rmsd = %5.2f"%(
      x, L.radii[x], L.tilecounts[x], L.mean_cv[x][0], L.mean_cv[x][1],
      L.tile_rmsd[x]
      )
  print "\nOverall                 %8d observations, delx %5.2f  dely %5.2f, rmsd = %5.2f"%(
      L.overall_N, L.overall_cv[0], L.overall_cv[1], L.overall_rmsd)
  print "Average tile rmsd %5.2f"%flex.mean(flex.double(L.tile_rmsd))
  print "Average tile displacement %5.2f"%(flex.mean(
    flex.double([matrix.col(cv).length() for cv in L.mean_cv])))

  if working_phil.show_plots is True:
    plt.plot([(L.tiles[4*x+0]+L.tiles[4*x+2])/2. for x in xrange(64)],[(L.tiles[4*x+1]+L.tiles[4*x+3])/2. for x in xrange(64)],"go")
    for x in xrange(64):
      plt.text(10+(L.tiles[4*x+0]+L.tiles[4*x+2])/2.,10+(L.tiles[4*x+1]+L.tiles[4*x+3])/2.,"%d"%x)
    plt.show()

    for idx in xrange(64):
      x = sort_radii[idx]
      print "Tile %2d: radius %7.2f, %6d observations, delx %5.2f  dely %5.2f, rmsd = %5.2f"%(
        x, L.radii[x], L.tilecounts[x], L.mean_cv[x][0], L.mean_cv[x][1],
        L.tile_rmsd[x]
        ),
      if L.tilecounts[x] < 3:
        print
      else:
        wtaveg = L.weighted_average_angle_deg_from_tile(x)
        print "Tile rotation %.2f degrees"%wtaveg

      from matplotlib import pyplot as plt
      xcv,ycv = get_correction_vector_xy(L,x)
      plt.plot(xcv,ycv,"r.")
      plt.plot([L.mean_cv[x][0]],[L.mean_cv[x][1]],"go")
      plt.axes().set_aspect("equal")
      plt.show()
