from __future__ import absolute_import, division, print_function
from six.moves import range
import math
import scitbx.math
from scitbx.array_family import flex
from rstbx.outlier_spots.fit_distribution import fit_cdf
from rstbx.indexing_api import rayleigh_cpp
from rstbx_ext import * # gets us SpotClass
from six.moves import zip

def format_data(x_data=None,y_data=None):
  """
  =============================================================================
  Function converts a list of separate x and y coordinates to a format suitable
  for plotting in ReportLab

  Arguments:
    x_data - a list of x coordinates (or any object that can be indexed)
    y_data - a list of y coordinates (or any object that can be indexed)

  Returns:
    A tuple of tuples containing the paired data

  Notes:
    The output from this function should be in a list for plotting in ReportLab.
    Multiple items in the list represent multiple data sets
  -----------------------------------------------------------------------------
  """
  assert(x_data is not None)
  assert(y_data is not None)
  assert(len(x_data) == len(y_data))

  # store data
  data = []
  for x,y in zip(x_data,y_data):
    data.append((x,y))
  data = tuple(data)

  return data

def make_histogram_data(data=None,n_bins=25):
  """
  =============================================================================
  Function generates a histogram given the data

  Arguments:
    data - data to be counted
    n_bins - number of bins to use in the histogram

  Returns:
    Two tuples are returned: upper limit of each bin, count in each bin

  Notes:
    data should already be sorted (minimum to maximum)
    the returned tuples will be of length n_bins
  -----------------------------------------------------------------------------
  """
  # determine binning
  min_value = float(math.floor(data[0]))
  max_value = 1.0
  d_bins = (max_value - min_value)/float(n_bins)
  bins = [0.0 for i in range(n_bins)]
  bins[0] = min_value + d_bins
  for i in range(1,n_bins):
    bins[i] = bins[i-1] + d_bins
  bins.append(float(math.ceil(data[-1])))
  n_bins = n_bins + 1
  # sort data into bins
  current_bin = 0
  histogram_data = [0 for i in range(n_bins)]
  for i in range(len(data)):
    if (data[i] < bins[current_bin]):
      histogram_data[current_bin] = histogram_data[current_bin] + 1
    else:
      while(data[i] > bins[current_bin]):
        current_bin = current_bin + 1

      histogram_data[current_bin] = histogram_data[current_bin] + 1
  if ((data[-1] < bins[-1]) and (data[-1] > bins[-2])):
    histogram_data[-1] = histogram_data[-1] + 1

  # fractionalize output
  for i in range(len(histogram_data)):
    histogram_data[i] = float(histogram_data[i])/float(len(data))

  return bins,histogram_data

class find_outliers:

  def __init__(self,ai=None,fraction=0.4,verbose=True,horizon_phil=None):

    # initialize variables
    self.fraction = fraction
    self.observed_spots = None
    self.predicted_spots = None
    self.fraction_spot_indices = None
    self.good = None

    self.sorted_observed_spots = None
    self.dr = None
    self.x = None
    self.not_good_dr = None
    self.dx = None
    self.dy = None

    self.plot_dxdy_data = None
    self.plot_cdf_data = None
    self.plot_pdf_data = None
    self.verbose = verbose
    self.cache_status = self.update(horizon_phil,ai,status_with_marked_outliers=None)

  def update(self,horizon_phil,ai,status_with_marked_outliers,verbose=False):
    # first time through: status with marked outliers is None; no input information
    # second time through after re-refinement: status is a list of SpotClass instances
    first_time_through = status_with_marked_outliers is None

    # update spot positions
    self.observed_spots = ai.raw_spot_input
    self.predicted_spots,current_status = ai.get_predicted_spot_positions_and_status(
                                          old_status = status_with_marked_outliers )

    # store Miller indices
    self.hkl = [0 for i in range(len(self.observed_spots))]
    for i in range(len(self.observed_spots)):
      self.hkl[i] = ai.get_hkl(i)
    return self.update_detail(horizon_phil,current_status,first_time_through,verbose)

  def update_detail(self,horizon_phil,current_status,first_time_through,verbose):
    assert(len(self.observed_spots) == len(self.predicted_spots))

    if horizon_phil.indexing.outlier_detection.verbose:
      classes=[str(current_status[i]) for i in range(len(self.observed_spots))]
      class_types = set(classes)
      class_counts = dict([[item,classes.count(item)] for item in class_types])
      flex_counts = flex.int(class_counts.values())
      assert flex.sum(flex_counts) == len(self.observed_spots)
      #for pair in class_counts.items():
      #  print "%10s %6d"%pair
      #print "%10s %6d"%("TOTAL",len(self.observed_spots))
      if status_with_marked_outliers == None:
        # status_with_marked_outliers==None is shorthand for identifying the first run through
        print("""After indexing on a subset of %d spots (from all images), %d were reclassified as
      either lying on the spindle, or potential overlapped spots or ice rings."""%(
      len(self.observed_spots),len(self.observed_spots)-class_counts["GOOD"]))
      else:
        print("""Rerefinement on just the well-fit spots followed by spot reclassification
      leaves %d good spots on which to calculate a triclinic rmsd."""%(class_counts["GOOD"]))

    # check good spots
    if (self.good is not None):
      match = 0
      for i in range(len(self.observed_spots)):
        if ((current_status[i] == SpotClass.GOOD) and self.good[i]):
          match = match + 1
      if self.verbose:print("Number of GOOD spots matched with previous model =",match)

    # calculate differences for all spots
    self.sorted_observed_spots = {}
    self.dr = flex.double()
    self.not_good_dr = flex.double()
    self.dx = [0.0 for i in range(len(self.observed_spots))]
    self.dy = [0.0 for i in range(len(self.observed_spots))]
    for i in range(len(self.observed_spots)):
      o = self.observed_spots[i]
      p = self.predicted_spots[i]
      self.dx[i] = o[0] - p[0]
      self.dy[i] = o[1] - p[1]
      self.sorted_observed_spots[
        math.sqrt(self.dx[i]*self.dx[i] + self.dy[i]*self.dy[i])] = i

    # separate GOOD spots
    spotclasses = {SpotClass.GOOD:0,SpotClass.SPINDLE:0,SpotClass.OVERLAP:0,SpotClass.ICE:0,SpotClass.OUTLIER:0,SpotClass.NONE:0}
    for key in sorted(self.sorted_observed_spots.keys()):
      spotclass = current_status[self.sorted_observed_spots[key]]
      spotclasses[spotclass]+=1
      if (current_status[self.sorted_observed_spots[key]] == SpotClass.GOOD):
        self.dr.append(key)
      else:
        self.not_good_dr.append(key)
    if verbose: print(", ".join(["=".join([str(i[0]),"%d"%i[1]]) for i in spotclasses.items()]), end=' ')
    totalsp = sum(spotclasses.values())
    if verbose: print("Total=%d"%(totalsp),"# observed spots",len(self.observed_spots))
    assert totalsp == len(self.observed_spots), "Some spot pairs have the same predicted-observed distances. Do you have duplicated images?"

    self.x = flex.double(len(self.dr))
    for i in range(len(self.x)):
      self.x[i] = float(i)/float(len(self.x))

    limit = int(self.fraction*len(self.dr))
    if limit < 4: return # Basic sanity check, need at least a few good spots to fit the distribution
    fitted_rayleigh = fit_cdf(x_data=self.dr[0:limit],
                              y_data=self.x[0:limit],distribution=rayleigh_cpp)
    if False:
        y_data=self.x[0:limit]
        inv_cdf = [fitted_rayleigh.distribution.inv_cdf(cdf) for cdf in y_data]
        from matplotlib import pyplot as plt
        plt.plot(self.dr[0:limit],self.x[0:limit],"r+")
        plt.plot(inv_cdf,y_data,"b.")
        plt.show()

    # store indices for spots used for fitting
    self.fraction_spot_indices = []
    for dr in self.dr[0:limit]:
      self.fraction_spot_indices.append(self.sorted_observed_spots[dr])

    # generate points for fitted distributions
    rayleigh_cdf_x = flex.double(range(500))
    rayleigh_cdf_x /= float(len(rayleigh_cdf_x))
    rayleigh_cdf = fitted_rayleigh.distribution.cdf(x=rayleigh_cdf_x)

    # generate points for pdf
    dr_bins,dr_histogram = make_histogram_data(data=self.dr,n_bins=100)
    rayleigh_pdf = flex.double(len(dr_bins))
    for i in range(len(dr_bins)):
      rayleigh_pdf[i] = fitted_rayleigh.distribution.pdf(x=dr_bins[i])
    rayleigh_pdf = rayleigh_pdf/flex.sum(rayleigh_pdf)
    dr_bins = flex.double(dr_bins)
    dr_histogram = flex.double(dr_histogram)

    # standard deviation for cdf
    sd = math.sqrt((4.0-math.pi)/(2.0)*
                   fitted_rayleigh.x[0]*fitted_rayleigh.x[0])
    if self.verbose:print('Standard deviation of Rayleigh fit = %4.3f'%sd)
    sd_data = None
    radius_outlier_index = None
    limit_outlier = None
    # --- Quoted code superceeded by extension module call to find_green_bar
    """
    for i in range(len(rayleigh_cdf_x)):
      mx = rayleigh_cdf_x[i]
      my = rayleigh_cdf[i]
      for j in range(1,len(self.dr)):
        upper_x = self.dr[j]
        upper_y = self.x[j]
        lower_x = self.dr[j-1]
        lower_y = self.x[j-1]
        if ((my >= lower_y) and (my < upper_y)):
          if ((sd <= (upper_x - mx)) and ((lower_x - mx) > 0.0)):
            sd_data = ((mx,my),(lower_x,lower_y))
            radius_outlier_index = j-1
            limit_outlier = lower_x
            if self.verbose:print "Width of green bar = %4.3f"%(lower_x - mx)
            break
        if (sd_data is not None):
          break
    """
    from rstbx.indexing_api import find_green_bar
    green = find_green_bar(rayleigh_cdf_x = rayleigh_cdf_x,
                           rayleigh_cdf = rayleigh_cdf,
                           dr = self.dr, x = self.x, sd = sd)
    if green.is_set:
      #assert radius_outlier_index == green.radius_outlier_index
      #assert limit_outlier == green.limit_outlier
      #assert sd_data[0][0] == green.sd_mx
      #assert sd_data[0][1] == green.sd_my
      #assert sd_data[1][0] == green.sd_lower_x
      #assert sd_data[1][1] == green.sd_lower_y
      if self.verbose:print("Width of green bar = %4.3f"%(green.sd_lower_x - green.sd_mx))
      radius_outlier_index = green.radius_outlier_index
      limit_outlier = green.limit_outlier
      sd_data = ((green.sd_mx, green.sd_my), (green.sd_lower_x, green.sd_lower_y))

    if (radius_outlier_index is None):
      radius_outlier_index = len(self.dr)
    if (limit_outlier is None):
      limit_outlier = self.dr[-1]
    radius_95 = None
    for i in range(len(rayleigh_cdf)):
      if (rayleigh_cdf[i] >= 0.95):
        radius_95 = rayleigh_cdf_x[i]
        break
    if (radius_95 is None):
      radius_95 = rayleigh_cdf_x[-1]
    upper_circle = []
    lower_circle = []
    d_radius = 2.0*radius_95/100.0
    x = -radius_95
    r2 = radius_95*radius_95
    for i in range(100):
      y = math.sqrt(r2 - x*x)
      upper_circle.append((x,y))
      lower_circle.append((x,-y))
      x = x + d_radius
    y = 0.0
    upper_circle.append((x,y))
    lower_circle.append((x,-y))
    self.sqrtr2 = math.sqrt(r2)

    # color code dx dy
    dxdy_fraction = []
    dxdy_inliers = []
    dxdy_outliers = []

    limit = self.dr[int(self.fraction*len(self.dr))]

    trifold = dict(fraction=0,inlier=0,outlier=0,total=0)
    for key in self.dr:
      trifold["total"]+=1
      i = self.sorted_observed_spots[key]
      if (key < limit):
        trifold["fraction"]+=1
        if (not ((self.dx[i] > 1.0) or (self.dx[i] < -1.0) or
                 (self.dy[i] > 1.0) or (self.dy[i] < -1.0))):
          dxdy_fraction.append((self.dx[i],self.dy[i]))
      elif (key < limit_outlier):
        trifold["inlier"]+=1
        if (not ((self.dx[i] > 1.0) or (self.dx[i] < -1.0) or
                 (self.dy[i] > 1.0) or (self.dy[i] < -1.0))):
          dxdy_inliers.append((self.dx[i],self.dy[i]))
      else:
        trifold["outlier"]+=1
        if (not ((self.dx[i] > 1.0) or (self.dx[i] < -1.0) or
                 (self.dy[i] > 1.0) or (self.dy[i] < -1.0))):
          dxdy_outliers.append((self.dx[i],self.dy[i]))
    if verbose: print(", ".join(["=".join([str(i[0]),"%d"%i[1]]) for i in trifold.items()]))

    # color code observed fractions
    o_fraction = []
    o_inliers = []
    o_outliers = []
    mr = format_data(x_data=rayleigh_cdf_x,y_data=rayleigh_cdf)
    limit = int(self.fraction*len(self.dr))
    for i in range(len(self.dr)):
      if (self.dr[i] <= 1.0):
        if (i < limit):
          o_fraction.append((self.dr[i],self.x[i]))
        elif (i < radius_outlier_index):
          o_inliers.append((self.dr[i],self.x[i]))
        else:
          o_outliers.append((self.dr[i],self.x[i]))
    if horizon_phil.indexing.outlier_detection.verbose:
      o_outliers_for_severity = []
      for i in range(radius_outlier_index, len(self.dr)):
        o_outliers_for_severity.append((self.dr[i],self.x[i]))

    # limit data range
    for i in range(len(dr_bins)):
      if (dr_bins[i] > 1.0):
        dr_bins.resize(i)
        dr_histogram.resize(i)
        rayleigh_pdf.resize(i)
        break
    ho = format_data(x_data=dr_bins,y_data=dr_histogram)
    hr = format_data(x_data=dr_bins,y_data=rayleigh_pdf)

    # format data for graphing
    self.plot_dxdy_data = [dxdy_fraction,dxdy_inliers,dxdy_outliers,
                           [(0.0,0.0)],[],[],[],[]]

    self.framework = {4:dict(status=SpotClass.SPINDLE),
                      5:dict(status=SpotClass.OVERLAP),
                      6:dict(status=SpotClass.OUTLIER),
                      7:dict(status=SpotClass.ICE),
    }
    for key in self.not_good_dr:
      i = self.sorted_observed_spots[key]
      status = current_status[i]
      if (not ((self.dx[i] > 1.0) or (self.dx[i] < -1.0) or
               (self.dy[i] > 1.0) or (self.dy[i] < -1.0))):
        statuskey = [k for k in self.framework.keys() if self.framework[k]["status"]==status][0]
        self.plot_dxdy_data[statuskey].append((self.dx[i],self.dy[i]))

    self.plot_cdf_data = [mr,o_fraction,o_inliers,o_outliers]
    if (sd_data is not None):
      self.plot_cdf_data.append(sd_data)
    self.plot_pdf_data = [ho,hr]

    # mark outliers
    if (first_time_through): #i.e., first time through the update() method
      if (radius_outlier_index < len(self.dr)):
        for i in range(radius_outlier_index,len(self.dr)):
          current_status[self.sorted_observed_spots[self.dr[i]]] = SpotClass.OUTLIER

    # reset good spots
    self.good = [False for i in range(len(self.observed_spots))]
    for i in range(len(self.observed_spots)):
      if (current_status[i] == SpotClass.GOOD):
        self.good[i] = True

    count_outlier = 0
    count_good = 0
    for i in range(len(self.observed_spots)):
      if (current_status[i] == SpotClass.OUTLIER):
        count_outlier = count_outlier + 1
      elif (current_status[i] == SpotClass.GOOD):
        count_good = count_good + 1
    if self.verbose:print('Old GOOD =', len(self.dr),\
          'OUTLIER =', count_outlier,\
          'New GOOD =', count_good)
    if horizon_phil.indexing.outlier_detection.verbose and status_with_marked_outliers is None:
      print("\nOf the remaining %d spots, %.1f%% were lattice outliers, leaving %d well-fit spots"%(
       len(self.dr),100.*count_outlier/len(self.dr), count_good ))
      if count_outlier==0:return
      #width of green bar is sd
      delta_spread = o_outliers_for_severity[1][1]-o_outliers_for_severity[0][1]
      severity = 0.
      for item in o_outliers_for_severity:
        delta_r = item[0] # obs - predicted deviation in mm
        spread = item[1] # order of observed deviation on a scale from 0 to 1
        # now invert the cdf to find expected delta r:
        expected_delta_r = fitted_rayleigh.distribution.sigma * math.sqrt(
          -2.* math.log(1.-spread) )
        #print item, expected_delta_r, (delta_r - expected_delta_r) / sd
        severity += ((delta_r - expected_delta_r) / sd)
      severity *= delta_spread
      print("The outlier severity is %.2f sigma [defined in J Appl Cryst (2010) 43, p.611 sec. 4].\n"%severity)
    return current_status


  # make plots
  def make_graphs(self,canvas=None,left_margin=None):#text=None):
    from reportlab.graphics import renderPDF
    from reportlab.lib.pagesizes import letter
    from reportlab.graphics.shapes import Drawing,String
    from reportlab.graphics.charts.legends import Legend
    from reportlab.graphics.charts.lineplots import LinePlot
    from reportlab.graphics.widgets.markers import makeMarker
    from reportlab.lib import colors
    from reportlab.lib.units import inch
    #help(colors)

    self.framework = {4:dict(status=SpotClass.SPINDLE,color=colors.black),
                      5:dict(status=SpotClass.OVERLAP,color=colors.limegreen),
                      6:dict(status=SpotClass.OUTLIER,color=colors.greenyellow),
                      7:dict(status=SpotClass.ICE,color=colors.skyblue),
    }


    # set size and position
    width,height = letter
    #letter_landscape = (width,height)
    plot_dim = 3.0*inch

    # construct scatter plot
    plot_dxdy_pos = (left_margin*inch,height - plot_dim - 0.5*inch)
    plot_dxdy = LinePlot()
    plot_dxdy.data = self.plot_dxdy_data

    std_colors = {0:colors.darkred, 1:colors.red, 2:colors.salmon}
    for key in std_colors.keys():
      plot_dxdy.lines[key].strokeColor = None
      plot_dxdy.lines[key].symbol = makeMarker('Circle')
      plot_dxdy.lines[key].symbol.strokeColor = None
      plot_dxdy.lines[key].symbol.fillColor = std_colors[key]
      plot_dxdy.lines[key].symbol.size = 1.2

    for key in self.framework.keys():
      plot_dxdy.lines[key].strokeColor = None
      plot_dxdy.lines[key].symbol = makeMarker('Circle')
      plot_dxdy.lines[key].symbol.strokeColor = None
      plot_dxdy.lines[key].symbol.fillColor = self.framework[key]["color"]
      plot_dxdy.lines[key].symbol.size = 1.2

    plot_dxdy.lines[3].strokeColor = None
    plot_dxdy.lines[3].symbol = makeMarker('Circle')
    plot_dxdy.lines[3].symbol.strokeColor = colors.blue
    plot_dxdy.lines[3].symbol.fillColor = None
    plot_dxdy.lines[3].symbol.strokeWidth = 0.6
    plot_dxdy.lines[3].symbol.size = plot_dim*(self.sqrtr2)
    #print plot_dxdy.lines[3].symbol.getProperties()
    plot_dxdy.width = plot_dim
    plot_dxdy.height = plot_dim
    plot_dxdy.xValueAxis.valueMax = 1.0
    plot_dxdy.xValueAxis.valueMin = -1.0
    plot_dxdy.xValueAxis.joinAxis = plot_dxdy.yValueAxis
    plot_dxdy.xValueAxis.joinAxisMode = 'value'
    plot_dxdy.xValueAxis.joinAxisPos = -1.0
    plot_dxdy.yValueAxis.valueMax = 1.0
    plot_dxdy.yValueAxis.valueMin = -1.0
    d_dxdy = Drawing(plot_dim,plot_dim)
    d_dxdy.add(plot_dxdy)

    # construct cdf plot
    plot_cdf_pos = (left_margin*inch, height - 2.0*(plot_dim + 0.5*inch))
    plot_cdf = LinePlot()
    plot_cdf.data = self.plot_cdf_data
    plot_cdf.lines[0].strokeColor = colors.blue

    for key in std_colors.keys():
      plot_cdf.lines[key+1].strokeColor = None
      plot_cdf.lines[key+1].symbol = makeMarker('Circle')
      plot_cdf.lines[key+1].symbol.strokeColor = None
      plot_cdf.lines[key+1].symbol.fillColor = std_colors[key]
      plot_cdf.lines[key+1].symbol.size = 1.2

    if (len(self.plot_cdf_data) == 5):
      plot_cdf.lines[4].strokeColor = colors.green
    plot_cdf.width = plot_dim
    plot_cdf.height = plot_dim
    plot_cdf.xValueAxis.valueMax = 1.0
    plot_cdf.xValueAxis.valueMin = 0.0
    plot_cdf.yValueAxis.valueMax = 1.0
    plot_cdf.yValueAxis.valueMin = 0.0
    d_cdf = Drawing(plot_dim,plot_dim)
    d_cdf.add(plot_cdf)

    # construct pdf plot
    plot_pdf_pos = (left_margin*inch, height - 3.0*(plot_dim + 0.5*inch))
    plot_pdf = LinePlot()
    plot_pdf.data = self.plot_pdf_data
    plot_pdf.lines[1].strokeColor = colors.blue
    plot_pdf.lines[0].strokeColor = None
    plot_pdf.lines[0].symbol = makeMarker('Circle')
    plot_pdf.lines[0].symbol.strokeColor = colors.red
    plot_pdf.lines[0].symbol.size = 1
    plot_pdf.width = plot_dim
    plot_pdf.height = plot_dim
    plot_pdf.xValueAxis.valueMax = 1.0
    plot_pdf.xValueAxis.valueMin = 0.0
    d_pdf = Drawing(2*plot_dim,plot_dim)
    d_pdf.add(plot_pdf)

    # add legend
    legend = Legend()
    legend.alignment = 'right'
    legend.colorNamePairs = [(std_colors[0],'Inliers (%d'%int(self.fraction*100.0) + '% used for fit)'),
                             (std_colors[1],'Other inliers'),
                             (std_colors[2],'Outliers, reject next round'),]
    for key in self.framework.keys():
      legend.colorNamePairs.append(  (self.framework[key]["color"], "%s"%self.framework[key]["status"]  )  )

    legend.x = plot_dim - 1.0*inch
    legend.y = plot_dim
    legend.columnMaximum = 8
    d_pdf.add(legend)

    # add titles
    title_pos = (plot_dim/2.0,plot_dim + 0.25*inch)
    title_dxdy = String(title_pos[0],title_pos[1],'dx vs. dy (all)')
    title_dxdy.fontSize = 15
    title_dxdy.textAnchor = 'middle'
    d_dxdy.add(title_dxdy)
    title_cdf = String(title_pos[0],title_pos[1],'cdf (good)')
    title_cdf.fontSize = 15
    title_cdf.textAnchor = 'middle'
    d_cdf.add(title_cdf)
    title_pdf = String(title_pos[0],title_pos[1],'pdf (good)')
    title_pdf.fontSize = 15
    title_pdf.textAnchor = 'middle'
    d_pdf.add(title_pdf)

    # draw everything
    renderPDF.draw(d_dxdy,canvas,plot_dxdy_pos[0],plot_dxdy_pos[1])
    renderPDF.draw(d_cdf,canvas,plot_cdf_pos[0],plot_cdf_pos[1])
    renderPDF.draw(d_pdf,canvas,plot_pdf_pos[0],plot_pdf_pos[1])

class find_outliers_from_matches(find_outliers):
  """Specifically rewritten ('RRR') to take dials refinery reflection manager matches.
     May have to be refactored again depending on what input list is finally chosen."""

  def get_new_status(self, old_status = None ):
    if old_status == None:
      return ( [SpotClass.GOOD] * len(self.predicted_spots) )
    else:
      return old_status

  def get_cache_status(self):
    value = flex.bool()
    if self.cache_status is None:
      raise Exception('No reflections left after removing outliers!')
    else:
      for status in self.cache_status:
        value.append( status==SpotClass.GOOD )
      return value

  def update(self,horizon_phil,matches,status_with_marked_outliers,verbose=False):
    # first time through: status with marked outliers is None; no input information
    # second time through after re-refinement: status is a list of SpotClass instances
    first_time_through = status_with_marked_outliers is None

    # update spot positions
    self.observed_spots = flex.vec2_double([(m.x_obs,m.y_obs) for m in matches])
    self.predicted_spots = flex.vec2_double([(m.x_calc,m.y_calc) for m in matches])

    # store Miller indices
    self.hkl = [m.miller_index for m in matches]

    current_status = self.get_new_status( old_status = status_with_marked_outliers )
    return self.update_detail(horizon_phil,current_status,first_time_through,verbose)
