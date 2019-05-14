from __future__ import absolute_import, division, print_function
from six.moves import range
import math,copy
import scitbx.math
from scitbx.array_family import flex
from xfel.metrology.mark0 import get_phil,correction_vectors
from scitbx import matrix
from xfel import get_radial_tangential_vectors
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in

class fit_translation(correction_vectors,lbfgs_with_curvatures_mix_in):
  def __init__(self,params):
    correction_vectors.__init__(self)
    self.read_data(params)
    self.x = flex.double([0.]*128) # x & y displacements for each of 64 tiles
    lbfgs_with_curvatures_mix_in.__init__(self,
      min_iterations=0,
      max_iterations=1000,
      use_curvatures=True)

  def curvatures(self):
    curvs = flex.double([0.]*128)
    for x in range(64):
      curvs[2*x] = 2. * self.selection_counts[x]
      curvs[2*x+1]=2. * self.selection_counts[x]
    return curvs

  def compute_functional_and_gradients(self):
    # HATTNE does never enter this function
    print("HATTNE entering mark1.compute_functional_and_gradients")

    self.model_calcx = self.spotcx.deep_copy()
    self.model_calcy = self.spotcy.deep_copy()

    for x in range(64):
      selection = self.selections [x]
      self.model_calcx.set_selected(selection, self.model_calcx + self.x[2*x])
      self.model_calcy.set_selected(selection, self.model_calcy + self.x[2*x+1])

    squares = self.delrsq_functional(self.model_calcx, self.model_calcy)

    f = flex.sum(squares)
    calc_obs_diffx = self.model_calcx - self.spotfx
    calc_obs_diffy = self.model_calcy - self.spotfy

    gradients = flex.double([0.]*128)
    for x in range(64):
      selection = self.selections [x]
      gradients[2*x] = 2. * flex.sum( calc_obs_diffx.select(selection) )
      gradients[2*x+1]=2. * flex.sum( calc_obs_diffy.select(selection) )
    print("Functional ",math.sqrt(flex.mean(squares)))

    return f,gradients

  def post_min_recalc(self):

    print("ENTERING post_min_recalc cx, cy", \
      flex.mean(self.spotcx), \
      flex.mean(self.spotcy), \
      len(self.spotcx), \
      len(self.spotcy))
    print("ENTERING post_min_recalc fx, fy", \
      flex.mean(self.spotfx), \
      flex.mean(self.spotfy), \
      len(self.spotfx), \
      len(self.spotfy))
    print("HATTNE check input 0", \
      flex.mean(self.model_calcx), \
      flex.mean(self.model_calcy), \
      len(self.model_calcx), \
      len(self.model_calcy))

    self.delrsq = self.delrsq_functional(calcx = self.model_calcx, calcy = self.model_calcy)
    self.tile_rmsd = [0.]*(len(self.tiles) // 4)
    self.asymmetric_tile_rmsd = [0.]*(len(self.tiles) // 4)

    self.correction_vector_x = self.model_calcx -self.spotfx
    self.correction_vector_y = self.model_calcy -self.spotfy

    self.post_mean_cv = []

    print("HATTNE post_min_recalc input CV x,y", \
      flex.mean(flex.pow(self.correction_vector_x, 2)), \
      flex.mean(flex.pow(self.correction_vector_y, 2)))
    print("HATTNE post_min_recalc input model ", \
      flex.mean(self.model_calcx), \
      flex.mean(self.model_calcy))

    for x in range(len(self.tiles) // 4):
        #if self.tilecounts[x]==0: continue
        selection = self.selections[x]
        selected_cv = self.master_cv.select(selection)

        if selection.count(True) == 0:
          self.post_mean_cv.append(matrix.col((0, 0)))
          self.asymmetric_tile_rmsd[x] = 0
          self.tile_rmsd[x] = 0
        else:
          self.post_mean_cv.append(
            matrix.col([flex.mean(self.correction_vector_x.select(selection)),
                        flex.mean(self.correction_vector_y.select(selection))]))
          self.asymmetric_tile_rmsd[x] = math.sqrt(flex.mean(self.delrsq.select(selection)))

          sel_delx = self.correction_vector_x.select(selection)
          sel_dely = self.correction_vector_y.select(selection)
          symmetric_offset_x = sel_delx - self.post_mean_cv[x][0]
          symmetric_offset_y = sel_dely - self.post_mean_cv[x][1]
          symmetricrsq = symmetric_offset_x*symmetric_offset_x + symmetric_offset_y*symmetric_offset_y

          """
          if x == 14:
            print "STATS FOR TILE 14"

            aa = list(self.tiles[x * 4:(x + 1)*4])
            print "EFFECTIVE tiling", aa
            print "EFFECTIVE center", \
              ((aa[0] + aa[2]) / 2, (aa[1] + aa[3]) / 2), \
              (aa[2] - aa[0], aa[3] - aa[1])

            print "sel_delx", list(sel_delx)
            #print "sel_dely", list(sel_dely)
            print "model_calcx", list(self.model_calcx.select(selection))
            #print "model_calcy", list(self.model_calcy.select(selection))
            print "spotfx", list(self.spotfx.select(selection))
            #print "spotfy", list(self.spotfy.select(selection))
            print "spotcx", list(self.spotcx.select(selection))
            #print "spotcy", list(self.spotcy.select(selection))

            print "  sel_delx          ", flex.min(sel_delx), \
              flex.mean(sel_delx), flex.max(sel_delx)
            print "  sel_dely          ", flex.min(sel_dely), \
              flex.mean(sel_dely), flex.max(sel_dely)

            print "  symmetric_offset_x", flex.min(symmetric_offset_x), \
              flex.mean(symmetric_offset_x), flex.max(symmetric_offset_x)
            print "  symmetric_offset_y", flex.min(symmetric_offset_y), \
              flex.mean(symmetric_offset_y), flex.max(symmetric_offset_y)
            print "  symmetric rsq     ", flex.min(symmetricrsq), \
              flex.mean(symmetricrsq), flex.max(symmetricrsq)
            print "  rmsd              ", math.sqrt(flex.mean(symmetricrsq))

            #import sys
            #sys.exit(0)
          """

          self.tile_rmsd[x] =math.sqrt(flex.mean(symmetricrsq))

    self.overall_N = flex.sum(flex.int( [int(t) for t in self.tilecounts] ))
    self.overall_cv = matrix.col([flex.mean ( self.correction_vector_x),
                                  flex.mean ( self.correction_vector_y) ])
    self.overall_rmsd = math.sqrt(flex.mean(self.delrsq_functional(self.model_calcx, self.model_calcy)))

    print("HATTNE post_min_recalc post_min_recalc 1", list(self.overall_cv))
    print("HATTNE post_min_recalc post_min_recalc 2", self.overall_rmsd)


  def print_table(self):
    from libtbx import table_utils
    from libtbx.str_utils import format_value
    table_header = ["Tile","Dist","Nobs","aRmsd","Rmsd","delx","dely","disp","rotdeg",
                    "Rsigma","Tsigma","Transx","Transy","DelRot"]
    table_data = []
    table_data.append(table_header)
    sort_radii = flex.sort_permutation(flex.double(self.radii))
    tile_rmsds = flex.double()
    radial_sigmas = flex.double(64)
    tangen_sigmas = flex.double(64)

    wtaveg = [0.]*64
    for x in range(64):
      if self.tilecounts[x] >= 3:
        wtaveg[x] = self.weighted_average_angle_deg_from_tile(x, self.post_mean_cv[x], self.correction_vector_x,
          self.correction_vector_y)

    for idx in range(64):
      x = sort_radii[idx]
      if self.tilecounts[x] < 3:
        radial = (0,0)
        tangential = (0,0)
        rmean,tmean,rsigma,tsigma=(0,0,1,1)
      else:
        radial,tangential,rmean,tmean,rsigma,tsigma = get_radial_tangential_vectors(self,x)

      # paired rotations of two ASICS on the same sensor
      if x%2==0:
        delrot = "%5.2f"%(wtaveg[x]-wtaveg[x+1])
      else:
        delrot = ""

      radial_sigmas[x]=rsigma
      tangen_sigmas[x]=tsigma
      table_data.append(  [
        format_value("%3d",   x),
        format_value("%7.2f", self.radii[x]),
        format_value("%6d",  self.tilecounts[x]),
        format_value("%5.2f", self.asymmetric_tile_rmsd[x]),
        format_value("%5.2f", self.tile_rmsd[x]),
        format_value("%5.2f", self.post_mean_cv[x][0]),
        format_value("%5.2f", self.post_mean_cv[x][1]),
        format_value("%5.2f", matrix.col(self.post_mean_cv[x]).length()),
        format_value("%6.2f", wtaveg[x]),
        format_value("%6.2f", rsigma),
        format_value("%6.2f", tsigma),
        format_value("%5.2f", self.x[2*x]),
        format_value("%5.2f", self.x[2*x+1]),
        copy.copy(delrot)
      ])
    table_data.append([""]*len(table_header))
    rstats = flex.mean_and_variance(radial_sigmas,self.tilecounts.as_double())
    tstats = flex.mean_and_variance(tangen_sigmas,self.tilecounts.as_double())
    table_data.append(  [
        format_value("%3s",   "ALL"),
        format_value("%s", ""),
        format_value("%6d",  self.overall_N),
        format_value("%5.2f", math.sqrt(flex.mean(self.delrsq))),
        format_value("%5.2f", self.overall_rmsd),
        format_value("%5.2f", self.overall_cv[0]),
        format_value("%5.2f", self.overall_cv[1]),
        format_value("%5.2f", flex.mean(flex.double([cv.length() for cv in self.post_mean_cv]))),
        format_value("%s", ""),
        format_value("%6.2f", rstats.mean()),
        format_value("%6.2f", tstats.mean()),
        format_value("%s", ""),
        format_value("%s", ""),
        format_value("%s", ""),
      ])

    print()
    print(table_utils.format(table_data,has_header=1,justify='center',delim=" "))


def run(args):

  work_params = get_phil(args)
  C = fit_translation(work_params)
  C.post_min_recalc()
  C.print_table()

  return None

if (__name__ == "__main__"):

  result = run(args=["mysql.runtag=for_may060corner","mysql.passwd=sql789",
                     "mysql.user=nick","mysql.database=xfelnks",
                     "effective_tile_boundaries=713, 437, 907, 622, 516, 438, 710, 623, 713, 650, 907, 835, 516, 650, 710, 835, 509, 18, 694, 212, 507, 215, 692, 409, 721, 19, 906, 213, 720, 215, 905, 409, 86, 231, 280, 416, 283, 231, 477, 416, 85, 19, 279, 204, 283, 19, 477, 204, 106, 444, 291, 638, 106, 640, 291, 834, 318, 443, 503, 637, 318, 640, 503, 834, 434, 849, 619, 1043, 436, 1046, 621, 1240, 647, 848, 832, 1042, 648, 1045, 833, 1239, 18, 1066, 212, 1251, 214, 1065, 408, 1250, 17, 853, 211, 1038, 213, 853, 407, 1038, 229, 1474, 414, 1668, 229, 1277, 414, 1471, 15, 1474, 200, 1668, 16, 1278, 201, 1472, 442, 1464, 636, 1649, 638, 1464, 832, 1649, 441, 1252, 635, 1437, 638, 1252, 832, 1437, 846, 1134, 1040, 1319, 1042, 1133, 1236, 1318, 845, 922, 1039, 1107, 1042, 922, 1236, 1107, 1060, 1542, 1245, 1736, 1060, 1346, 1245, 1540, 848, 1543, 1033, 1737, 847, 1348, 1032, 1542, 1469, 1336, 1663, 1521, 1272, 1337, 1466, 1522, 1472, 1550, 1666, 1735, 1274, 1549, 1468, 1734, 1460, 1117, 1645, 1311, 1460, 921, 1645, 1115, 1247, 1117, 1432, 1311, 1248, 921, 1433, 1115, 1130, 718, 1315, 912, 1130, 522, 1315, 716, 918, 719, 1103, 913, 917, 523, 1102, 717, 1543, 514, 1737, 699, 1346, 513, 1540, 698, 1543, 725, 1737, 910, 1346, 725, 1540, 910, 1338, 94, 1523, 288, 1339, 290, 1524, 484, 1552, 93, 1737, 287, 1551, 289, 1736, 483, 1115, 114, 1309, 299, 918, 113, 1112, 298, 1115, 326, 1309, 511, 918, 326, 1112, 511".replace(" ",""),
  ])
  """mark0: no minimization; just evaluate tiles like cxi.plot_cv
     mark1: lsq fitting of translation of spot predictions, to bring them on top of
             the fictitious tile montage containing spotfinder spots.
     mark2: mark1 plus tile rotations following the translation step; spinning around tile center.
  """
