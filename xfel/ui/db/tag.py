from __future__ import division

class Tag(object):
  def __init__(self, tag_id = None, **kwargs):
    if tag_id is None:
      tag_id = 4 # create a new tag

    self.tag_id = tag_id

    # dummy values:
    self.name = ""
    self.comment = ""
