from __future__ import division

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
