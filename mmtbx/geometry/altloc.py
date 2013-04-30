from __future__ import division

def Empty(obj, processor):
  """
  An object without an altloc identifier
  """

  processor.process_regular( obj )


class Alternate(object):
  """
  An object assigned with an altloc identifier
  """

  def __init__(self, identifier):

    self.identifier = identifier


  def __call__(self, obj, processor):

    processor.process_altloc( obj, identifier = self.identifier )


def from_identifier(entity):

  if entity:
    return Alternate( identifier = entity )

  else:
    return Empty


def from_atom(entity):

  parent = entity.parent()

  if parent:
    return from_identifier( parent.altloc )

  else:
    return Empty

