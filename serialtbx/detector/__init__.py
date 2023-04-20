from __future__ import division

def center(coords):
  """ Returns the average of a list of vectors
  @param coords List of vectors to return the center of
  """
  for c in coords:
    if 'avg' not in locals():
      avg = c
    else:
      avg += c
  return avg / len(coords)

class basis(object):
  """ Bucket for detector element information """
  def __init__(self, orientation = None, translation = None, panelgroup = None, homogenous_transformation = None, name = None):
    """
    Provide only orientation + translation or a panelgroup or a homogenous_transformation.

    @param orientation rotation in the form of a quarternion
    @param translation vector translation in relation to the parent frame
    @param panelgroup dxtbx panelgroup object whose local d matrix will represent the
    basis shift
    @param homogenous_transformation 4x4 matrix.sqr object representing a translation
    and a rotation. Must not also contain a scale as this won't be decomposed properly.
    @param name optional name for this basis shift
    """
    self.include_translation = True
    self.name = name

    if panelgroup is not None:
      d_mat = panelgroup.get_local_d_matrix()
      fast = matrix.col((d_mat[0],d_mat[3],d_mat[6])).normalize()
      slow = matrix.col((d_mat[1],d_mat[4],d_mat[7])).normalize()
      orig = matrix.col((d_mat[2],d_mat[5],d_mat[8]))

      v3 = fast.cross(slow).normalize()

      r3 = matrix.sqr((fast[0],slow[0],v3[0],
                       fast[1],slow[1],v3[1],
                       fast[2],slow[2],v3[2]))

      self.orientation = r3.r3_rotation_matrix_as_unit_quaternion()
      self.translation = orig

      if not self.name:
        self.name = panelgroup.get_name()

    elif orientation is not None or translation is not None:
      assert orientation is not None and translation is not None
      self.orientation = orientation
      self.translation = translation

    else:
      # Decompose the homegenous transformation assuming no scale factors were used
      h = homogenous_transformation
      self.orientation = matrix.sqr((h[0],h[1],h[2],
                                     h[4],h[5],h[6],
                                     h[8],h[9],h[10])).r3_rotation_matrix_as_unit_quaternion()
      self.translation = matrix.col((h[3],
                                     h[7],
                                     h[11]))
      assert h[12] == h[13] == h[14] == 0 and h[15] == 1

def iterate_detector_at_level(item, depth = 0, level = 0):
  """
  Iterate through all panel groups or panels of a detector object at a given
  hierarchy level
  @param item panel group or panel. Use detector.hierarchy().
  @param depth current dept for recursion. Should be 0 for initial call.
  @param level iterate groups at this level
  @return next panel or panel group object
  """
  if level == depth:
    yield item
  else:
    for child in item:
      for subitem in iterate_detector_at_level(child, depth+1, level):
        yield subitem

def iterate_panels(panelgroup):
  """
  Find and iterate all panels in the given panel group, regardless of the hierarchly level
  of this panelgroup
  @param panelgroup the panel group of interest
  @return the next panel
  """
  if panelgroup.is_group():
    for child in panelgroup:
      for subitem in iterate_panels(child):
        yield subitem
  else:
    yield panelgroup

def id_from_name(detector, name):
  """ Jiffy function to get the id of a panel using its name
  @param detector detector object
  @param name panel name
  @return index of panel in detector
  """
  return [p.get_name() for p in detector].index(name)

def get_center(pg):
  """ Find the center of a panel group pg, projected on its fast/slow plane """
  if pg.is_group():
    # find the average center of all this group's children
    children_center = col((0,0,0))
    count = 0
    for p in iterate_panels(pg):
      children_center += get_center(p)
      count += 1
    children_center /= count

    # project the children center onto the plane of the panel group
    pgf = col(pg.get_fast_axis())
    pgs = col(pg.get_slow_axis())
    pgn = col(pg.get_normal())
    pgo = col(pg.get_origin())

    return (pgf.dot(children_center) * pgf) + (pgs.dot(children_center) * pgs) + (pgn.dot(pgo) * pgn)
  else:
    s = pg.get_image_size()
    return col(pg.get_pixel_lab_coord((s[0]/2, s[1]/2)))
