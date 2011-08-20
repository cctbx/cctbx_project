from scitbx.array_family import flex
from scitbx import matrix
"""Standalone program gives tile translations for the CXI pad detector;
     based on lysozyme test powder arcs.
   Fixed coordinates for the detector center (850,850) are simply based on CXI output
     file size.
   distl.detector_tiling limits are based on a separate run of the program
     spotfinder.find_active_area
   coordinates of pixels on the powder arcs of all four quadrants were obtained by
     manual inspection.  200 XFEL shots with the highest count of spotfinder spots
     were summed and inspected with the phenix.image_viewer.
"""

init_ctr_slow = 850
init_ctr_fast = 850
init_ctr = matrix.col([init_ctr_slow,init_ctr_fast])

class quadrant:
  def __init__(self,arcs):
    self.arcs = arcs

  def refine_center_from_arcs(self):
    #initial guess for center, based on given size of Numpy arra
    init_ctr_slow = 850
    init_ctr_fast = 850
    init_ctr = matrix.col([init_ctr_slow,init_ctr_fast])
    self.x = flex.double([init_ctr_slow,init_ctr_fast])
    #add a radius parameter for each arc
    for arc in self.arcs:
      sample_point = matrix.col(arc.points[0])
      guess_radius = (sample_point - init_ctr).length()
      self.x.append(guess_radius)

    self.initial_functional = None
    import scitbx.lbfgs
    scitbx.lbfgs.run(target_evaluator=self)
    self.final_functional = self.compute_functional_and_gradients(
      functional_only=True)

  def compute_functional_and_gradients(self, functional_only=False):
    def get_f():
      center = matrix.col([self.x[0],self.x[1]])
      residual = 0.0
      for arc,radius in zip(self.arcs,self.x[2:]):
        for point in arc.points:
          ppoint = matrix.col(point)
          rad_vector = ppoint - center
          length = rad_vector.length()
          residual += (length-radius)*(length-radius)
      return residual

    f = get_f()
    if (self.initial_functional is None):
      self.initial_functional = f
    if (functional_only):
      return f
    g = flex.double()
    g.reserve(len(self.x))
    eps = 1e-6
    for i in xrange(len(self.x)):
      xi = self.x[i]
      self.x[i] = xi+eps
      f_eps = get_f()
      self.x[i] = xi
      g.append((f_eps-f)/eps)
    return f, g

  def show_summary(self):
    print "The center is %7.2f %7.2f; "%(self.x[0],self.x[1]),
    print "radii are",[round(xx,2) for xx in self.x[2:]]

  def get_tile_translation(self):
    slow_translation = int(float(init_ctr_slow) - self.x[0])
    fast_translation = int(float(init_ctr_fast) - self.x[1])
    return (slow_translation,fast_translation)

def parse(string):
  words = string.split()
  values = []
  for x in xrange(len(words)//2):
    values.append((int(words[2*x]),int(words[2*x+1])))
  return values

class powder_arc:
  def __init__(self,points):
    self.points = flex.vec2_double(parse(points))
    print list(self.points)

def lysozyme_calibration():
  quad_1_UL = quadrant([
    powder_arc("""
    846 728
    814 731
    790 740
    760 761"""),
    powder_arc("""
    836 755
    822 757
    802 763
    788 771"""),
    powder_arc("""
    787 794
    799 784
    810 776
    825 772"""),
    powder_arc("""
    846 793
    832 794
    822 797
    812 802"""),
  ])
  quad_2_UR = quadrant([
    powder_arc("""
    800 853
    800 875
    803 885
    812 898
    """),
    powder_arc("""
    775 864
    778 884
    787 902
    799 915
    """),
    powder_arc("""
    759 863
    763 889
    773 908
    783 922
    """),
    powder_arc("""
    734 873
    740 900
    748 918
    809 968
    """),
    powder_arc("""
    697 884
    718 943
    768 992
    802 1005
    """),
  ])
  quad_3_LR = quadrant([
    powder_arc("""
    854 925
    867 925
    886 919
    896 914"""),
    powder_arc("""
    867 937
    880 933
    896 927
    909 919
    """),
    powder_arc("""
    923 941
    904 954
    885 959
    860 963
    """),
    powder_arc("""
    973 943
    959 959
    939 977
    883 998
    """),
  ])
  quad_4_LL = quadrant([
    powder_arc("""
    914 700
    948 719
    998 790
    1007 831
    """),
    powder_arc("""
    929 753
    952 780
    964 807
    968 831
    """),
    powder_arc("""
    942 841
    941 815
    931 795
    897 762
    """),
    powder_arc("""
    897 776
    913 789
    928 815
    931 838
    """),
  ])
  quadrants = [quad_1_UL,quad_2_UR,quad_3_LR,quad_4_LL]
  for quad in quadrants:
    quad.refine_center_from_arcs()
    quad.show_summary()
  derive_tile_translations(quadrants)

def derive_tile_translations(quads):
  params = get_initial_cxi_scope()
  tile_list=[]
  tile_translations=[]
  tile_flags=[]
  for itile in xrange(len(params.distl.detector_tiling)//4):
    corner_UL = matrix.col([params.distl.detector_tiling[itile*4],
                            params.distl.detector_tiling[itile*4+1]])
    corner_LR = matrix.col([params.distl.detector_tiling[itile*4+2],
                            params.distl.detector_tiling[itile*4+3]])
    middle = (corner_UL + corner_LR)/2.
    if middle[0]<init_ctr_slow and middle[1]<init_ctr_fast:
      tile_list.append(1)      #quadrant 1, UL

    elif middle[0]<init_ctr_slow and middle[1]>init_ctr_fast:
      tile_list.append(2)      #quadrant 2, UR

    elif middle[0]>init_ctr_slow and middle[1]>init_ctr_fast:
      tile_list.append(3)      #quadrant 3, LR

    else:
      tile_list.append(4)      #quadrant 4, LL
    translation = quads[tile_list[-1]-1].get_tile_translation()
    tile_translations.append(translation[0])
    tile_translations.append(translation[1])

    #finished with tile translations; now output flags for the actual tiles
    # where powder arcs have been observed (the inner four tiles)
    this_quad = quads[tile_list[-1]-1]
    first_point = this_quad.arcs[0].points[0]
    if corner_UL[0]<first_point[0] and first_point[0]<corner_LR[0] and\
       corner_UL[1]<first_point[1] and first_point[1]<corner_LR[1]:
         tile_flags.append(1)
    else: tile_flags.append(0)
  print "distl.tile_translations=%s"%(",".join([str(t) for t in tile_translations]))
  print "distl.tile_flags=%s"%(",".join([str(t) for t in tile_flags]))
  print "OK"

def get_initial_cxi_scope():
  limits="""(1479, 1515) (1672, 1699)
(1281, 1515) (1474, 1699)
(1092, 1506) (1249, 1699)
(1065, 1506) (1090, 1699)
(853, 1505) (1037, 1698)
(1650, 1081) (1672, 1487)
(1479, 1303) (1604, 1487)
(1281, 1303) (1474, 1487)
(1466, 1082) (1650, 1274)
(1065, 1308) (1249, 1501)
(1253, 1080) (1437, 1273)
(853, 1307) (1037, 1500)
(622, 1465) (815, 1649)
(424, 1465) (617, 1649)
(213, 1473) (397, 1666)
(1048, 1098) (1241, 1282)
(850, 1098) (1043, 1282)
(623, 1252) (816, 1436)
(425, 1252) (618, 1436)
(213, 1275) (397, 1468)
(1466, 883) (1650, 1076)
(1506, 664) (1699, 848)
(1253, 882) (1437, 1075)
(1048, 885) (1241, 1069)
(1308, 664) (1501, 848)
(1099, 656) (1283, 849)
(850, 885) (1043, 1069)
(634, 1048) (818, 1241)
(421, 1048) (605, 1241)
(198, 1064) (391, 1248)
(1506, 451) (1699, 635)
(1308, 451) (1501, 635)
(1099, 458) (1283, 651)
(1514, 231) (1698, 424)
(1085, 262) (1278, 446)
(1514, 33) (1698, 226)
(885, 656) (1069, 849)
(885, 458) (1069, 651)
(887, 262) (1080, 446)
(634, 850) (818, 1043)
(421, 850) (605, 1043)
(656, 626) (849, 810)
(656, 414) (849, 598)
(664, 198) (848, 391)
(664, 0) (848, 193)
(458, 626) (651, 810)
(198, 851) (391, 1035)
(458, 414) (651, 598)
(451, 198) (635, 391)
(451, 0) (635, 193)
(1, 1473) (185, 1666)
(1, 1275) (185, 1468)
(0, 1064) (193, 1248)
(0, 851) (193, 1035)
(263, 617) (447, 810)
(50, 616) (234, 809)
(263, 419) (447, 612)
(50, 418) (234, 611)
(231, 213) (424, 397)
(33, 213) (226, 397)"""
  ilimits = []
  for line in limits.split("\n"):
    for ituple in line.split(") ("):
      for dint in ituple.split(" "):
        ilimits.append(dint.replace(")","").replace("(","").replace(",",""))

  args = [
          "distl.bins.verbose=True",
          "distl.minimum_spot_area=3",
          "distl.detector_tiling=%s"%str(",".join(ilimits)),
          "distl.peripheral_margin=1",
          "force_method2_resolution_limit=2.1",
          "distl_highres_limit=2.1",
          "distl.res.outer=2.1",
          ]

  from labelit import preferences
  from rstbx.command_line.index import special_defaults_for_new_horizons
  new_horizons_phil = preferences.RunTimePreferences()
  special_defaults_for_new_horizons( new_horizons_phil )
  new_horizons_phil.merge_command_line(args)
  horizons_phil=new_horizons_phil.command_extractor
  return horizons_phil

if __name__=="__main__":
  lysozyme_calibration()
