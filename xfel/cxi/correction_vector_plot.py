import os,math
from scitbx import matrix
from cctbx.array_family import flex

class lines:
  def __init__(self,params):
    self.params = params
    self.tiles = None
  def literals(self):
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
              self.tiles = [int(a) for a in line.strip().split()[2:]]
              assert len(self.tiles)==256
              print self.tiles

  def vectors(self):
    self.tilecounts = [0.0001]*64
    self.radii = [0.]*64
    self.mean_cv = [matrix.col((0.,0.))]*64
    self.all_cv = [[] for x in xrange(64)]
    self.overall_cv = matrix.col((0.,0.))
    self.sum_sq_cv = 0.
    self.tile_rmsd = [0.]*64
    for line,run in self.literals():
     try:
      tokens = line.split()
      if len(tokens)>13:
        continue # guard against subprocess contention
      observcen = matrix.col((float(tokens[2]),float(tokens[3])))
      refinecen = matrix.col((float(tokens[5]),float(tokens[6])))
      observspo = matrix.col((float(tokens[8]),float(tokens[9])))
      predicspo = matrix.col((float(tokens[11]),float(tokens[12])))
      prediction = predicspo-refinecen
      #observation = observspo-observcen
      cv = predicspo-observspo
      for x in xrange(len(self.tiles)//4):
        if self.tiles[4*x+0]<predicspo[0]<self.tiles[4*x+2] and \
           self.tiles[4*x+1]<predicspo[1]<self.tiles[4*x+3]:
           break
      if cv.length() > 10:
        print "problem",cv.length(), cv, x, line, run
      itile = x
      # $$$ Acount for shadow obstruction in thermolysin runs
      # now handled by xfel_targets implementation; not needed here
        #if (run in [72,73]) and (itile in [41,43,46,47,53,55,57]): continue
        #if (run in [21,22,23,24,25,26,27,72,73]) and (itile in [40,42,44,45,52,54,58,59]): continue
      # $$$ Done accounting for shadow
      self.tilecounts[itile]+=1
      self.radii[itile]+=prediction.length()
      self.mean_cv[itile] = self.mean_cv[itile] + cv
      self.all_cv[itile].append(cv)
      self.overall_cv += cv
      self.sum_sq_cv += cv.length_sq()

      yield prediction, cv, itile
     except ValueError:
       print "Valueerror"
    for x in xrange(64):
      self.radii[x]/=self.tilecounts[x]
      self.mean_cv[x]/=self.tilecounts[x]
      if len(self.all_cv[x])>0:
        self.tile_rmsd[x] = math.sqrt(
      flex.mean(flex.double([ (cv - self.mean_cv[x]).length_sq() for cv in self.all_cv[x] ]))
      )
      else: self.tile_rmsd[x]=0.
    self.overall_N = flex.sum(flex.int( [int(t) for t in self.tilecounts] ))
    self.overall_cv /= self.overall_N
    self.overall_rmsd = math.sqrt( self.sum_sq_cv / self.overall_N )

def run_correction_vector_plot(working_phil):

  L = lines(working_phil)
  master_coords = []
  master_cv     = []
  master_tiles  = flex.int()
  for line in L.vectors():
    pass # pull out the information, lines class does all the work
    #assemble a master plot
    master_coords.append(line[0])
    master_cv.append(line[1])
    master_tiles.append(line[2])
  close_x = flex.double()
  close_y = flex.double()
  far_x = flex.double()
  far_y = flex.double()
  for idx in xrange(0,len(master_coords),10):
    if master_cv[idx].length() < L.tile_rmsd[ master_tiles[idx] ]:
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
  print "Average tile displacement %5.2f"%(flex.mean(flex.double([cv.length() for cv in L.mean_cv])))

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
        )

      from matplotlib import pyplot as plt
      xcv = [point[0] for point in L.all_cv[x]]
      ycv = [point[1] for point in L.all_cv[x]]
      plt.plot(xcv,ycv,"r.")
      plt.plot([L.mean_cv[x][0]],[L.mean_cv[x][1]],"go")
      plt.show()
