from __future__ import division
from __future__ import print_function
from six.moves import range
import os,math
from scitbx import matrix
from cctbx.array_family import flex
from xfel import correction_vector_store,get_correction_vector_xy
from xfel import get_radial_tangential_vectors

class manage_sql:
  def __init__(self,params):
    self.params = params
    self.have_db = False
    if self.params.mysql.runtag is not None:
      self.have_db = True
      from xfel.merging.database.merging_database import manager
      self.manager = manager

  def get_cursor(self):
    CART = self.manager(self.params)
    db = CART.connection()
    cursor = db.cursor()
    return cursor

  def initialize_tables_and_insert_command(self):
    if self.have_db is not True: return
    CART = self.manager(self.params)
    db = CART.connection()
    cursor = db.cursor()
    new_tables = CART.positional_refinement_schema_tables(self.params.mysql.runtag)
    for table in new_tables:
      cursor.execute("DROP TABLE IF EXISTS %s;"%table[0])
      cursor.execute("CREATE TABLE %s "%table[0]+table[1].replace("\n"," ")+" ;")
    from six.moves import cStringIO as StringIO
    self.query = StringIO()
    self.query.write("INSERT INTO %s_spotfinder VALUES "%self.params.mysql.runtag)
    self.firstcomma = ""

  def get_frame_dictionary(self):
    if self.have_db is not True: return dict()
    CART = self.manager(self.params)
    db = CART.connection()
    cursor = db.cursor()
    cursor.execute("SELECT unique_file_name,frame_id_1_base FROM %s_frame"%(self.params.mysql.runtag))
    del CART
    frame_dict = {}
    for ftuple in cursor.fetchall():
      frame_dict[ftuple[0]] = ftuple[1]
    return frame_dict

  def insert(self,run,itile,tokens):
    self.query.write(self.firstcomma); self.firstcomma=","
    self.query.write("('%d','%d','%10.2f','%10.2f','%10.2f','%10.2f','%10.2f','%10.2f','%10.2f','%10.2f',"%(
          run,itile,float(tokens[2]),float(tokens[3]),float(tokens[5]),float(tokens[6]),
          float(tokens[8]),float(tokens[9]),float(tokens[11]),float(tokens[12]) ))
    self.query.write("'%4d','%4d','%4d','%6.3f','%6.3f')"%(
      int(tokens[14]),int(tokens[15]),int(tokens[16]),float(tokens[19]),float(tokens[21])))


  def send_insert_command(self):
    if self.have_db is not True: return
    cursor = self.get_cursor()
    cursor.execute( self.query.getvalue() )

class lines(correction_vector_store):
  def __init__(self,params):
    correction_vector_store.__init__(self)
    self.params = params
    self.database = manage_sql(self.params)

  def literals(self):
    frame_dict = self.database.get_frame_dictionary()

    if self.params.run_numbers is None:
      path = self.params.outdir_template
      stream = open(path,"r")
      print(path)
      for line in stream:
        if line.find("XFEL processing:") == 0:
           tokens = line.strip().split("/")
           picklefile = line.strip().split()[2]
           frame_id = frame_dict.get(picklefile, None)
           if frame_id is not None:
             print("FETCHED",frame_id)
        if line.find("CV OBSCENTER")==0 and line.find("Traceback")==-1:
          potential_tokens = line.strip().split()
          if len(potential_tokens)==22 and \
                 potential_tokens[17].isdigit() and \
                 int(potential_tokens[17])==self.params.bravais_setting_id:
            yield frame_id,potential_tokens
        if len(self.tiles)==0 and line.find("EFFEC")==0:
          self.tiles = flex.int([int(a) for a in line.strip().split()[2:]])
          assert len(self.tiles)==256
          print(list(self.tiles))
          self.initialize_per_tile_sums()
      return
    for run in self.params.run_numbers:
        templ = self.params.outdir_template%run
        items = os.listdir(templ)
        for item in items:
          path = os.path.join(templ,item)
          stream = open(path,"r")
          print(path)
          for line in stream:
            if line.find("CV OBSCENTER")==0:
              potential_tokens = line.strip().split()
              if len(potential_tokens)==22 and \
                 potential_tokens[17].isdigit() and \
                 int(potential_tokens[17])==self.params.bravais_setting_id:
                yield None,potential_tokens
            if len(self.tiles)==0 and line.find("EFFEC")==0:
              self.tiles = flex.int([int(a) for a in line.strip().split()[2:]])
              assert len(self.tiles)==256
              print(list(self.tiles))
              self.initialize_per_tile_sums()

  def vectors(self):
    self.database.initialize_tables_and_insert_command()

    self.tile_rmsd = [0.]*64

    for run,tokens in self.literals():
     try:
      itile = self.register_line( float(tokens[2]),float(tokens[3]),
                       float(tokens[5]),float(tokens[6]),
                       float(tokens[8]),float(tokens[9]),
                       float(tokens[11]),float(tokens[12]) )
      if run is not None:
        self.database.insert(run,itile,tokens)
      yield "OK"
     except ValueError:
       print("Valueerror")

    self.database.send_insert_command()
    for x in range(64):
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
  for idx in range(0,len(master_coords),10):
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
  radial_sigmas = flex.double(64)
  tangen_sigmas = flex.double(64)
  for idx in range(64):
    x = sort_radii[idx]
    print("Tile %2d: radius %7.2f, %6d observations, delx %5.2f  dely %5.2f, rmsd = %5.2f"%(
      x, L.radii[x], L.tilecounts[x], L.mean_cv[x][0], L.mean_cv[x][1],
      L.tile_rmsd[x]
        ), end=' ')
    if L.tilecounts[x] < 3:
      print()
      radial = (0,0)
      tangential = (0,0)
      rmean,tmean,rsigma,tsigma=(0,0,1,1)
    else:
      wtaveg = L.weighted_average_angle_deg_from_tile(x)
      print("Tile rotation %6.2f deg"%wtaveg, end=' ')
      radial,tangential,rmean,tmean,rsigma,tsigma = get_radial_tangential_vectors(L,x)
      print("%6.2f %6.2f"%(rsigma,tsigma))
    radial_sigmas[x]=rsigma
    tangen_sigmas[x]=tsigma
  rstats = flex.mean_and_variance(radial_sigmas,L.tilecounts.as_double())
  tstats = flex.mean_and_variance(tangen_sigmas,L.tilecounts.as_double())

  print("\nOverall                 %8d observations, delx %5.2f  dely %5.2f, rmsd = %5.2f"%(
      L.overall_N, L.overall_cv[0], L.overall_cv[1], L.overall_rmsd))
  print("Average tile rmsd %5.2f"%flex.mean(flex.double(L.tile_rmsd)))
  print("Average tile displacement %5.2f"%(flex.mean(
    flex.double([matrix.col(cv).length() for cv in L.mean_cv]))))
  print("Weighted average radial sigma %6.2f"%rstats.mean())
  print("Weighted average tangential sigma %6.2f"%tstats.mean())

  if working_phil.show_plots is True:
    plt.plot([(L.tiles[4*x+0]+L.tiles[4*x+2])/2. for x in range(64)],[(L.tiles[4*x+1]+L.tiles[4*x+3])/2. for x in range(64)],"go")
    for x in range(64):
      plt.text(10+(L.tiles[4*x+0]+L.tiles[4*x+2])/2.,10+(L.tiles[4*x+1]+L.tiles[4*x+3])/2.,"%d"%x)
    plt.show()

    for idx in range(64):
      x = sort_radii[idx]
      print("Tile %2d: radius %7.2f, %6d observations, delx %5.2f  dely %5.2f, rmsd = %5.2f"%(
        x, L.radii[x], L.tilecounts[x], L.mean_cv[x][0], L.mean_cv[x][1],
        L.tile_rmsd[x]
        ), end=' ')
      if L.tilecounts[x] < 3:
        print()
        radial = (0,0)
        tangential = (0,0)
        rmean,tmean,rsigma,tsigma=(0,0,1,1)
      else:
        wtaveg = L.weighted_average_angle_deg_from_tile(x)
        print("Tile rotation %6.2f deg"%wtaveg, end=' ')
        radial,tangential,rmean,tmean,rsigma,tsigma = get_radial_tangential_vectors(L,x)
        print("%6.2f %6.2f"%(rsigma,tsigma))

      if working_phil.colormap:
        from pylab import imshow, axes, colorbar, show
        import numpy

        xcv,ycv = get_correction_vector_xy(L,x)
        _min = min(min(xcv),min(ycv))
        _max = max(max(xcv),max(ycv))

        hist,xedges,yedges = numpy.histogram2d(xcv,ycv,bins=40,range=[[_min,_max],[_min,_max]])
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]

        imshow(hist.T,extent=extent,interpolation='nearest',origin='lower')
        from matplotlib.patches import Ellipse
        ell = Ellipse(xy=(L.mean_cv[x][0],L.mean_cv[x][1]),
                      width=2.*rsigma, height=2.*tsigma,
                      angle=math.atan2(-(radial[1]),-(radial[0]))*180./math.pi,
                      edgecolor="y", linewidth=2, fill=False, zorder=100)
        axes().add_artist(ell)
        colorbar()
        show()

      else:
        from matplotlib import pyplot as plt
        xcv,ycv = get_correction_vector_xy(L,x)
        if len(xcv)==0 or len(ycv)==0: continue
        plt.plot(xcv,ycv,"r.")
        plt.plot([L.mean_cv[x][0]],[L.mean_cv[x][1]],"go")
        plt.plot([L.mean_cv[x][0]+radial[0]],[L.mean_cv[x][1]+radial[1]],"yo")
        plt.plot([L.mean_cv[x][0]+tangential[0]],[L.mean_cv[x][1]+tangential[1]],"bo")
        from matplotlib.patches import Ellipse
        ell = Ellipse(xy=(L.mean_cv[x][0],L.mean_cv[x][1]),
                      width=2.*rsigma, height=2.*tsigma,
                      angle=math.atan2(-(radial[1]),-(radial[0]))*180./math.pi,
                      edgecolor="y", linewidth=2, fill=False, zorder=100)
        plt.axes().add_artist(ell)
        plt.axes().set_aspect("equal")
        _min = min(min(xcv),min(ycv))
        _max = max(max(xcv),max(ycv))
        plt.axes().set_xlim(_min,_max)
        plt.axes().set_ylim(_min,_max)
        plt.show()
