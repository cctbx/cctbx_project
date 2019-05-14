from __future__ import absolute_import, division, print_function

def Empty(data, coordinates, processor):
  """
  An object without an altloc identifier
  """

  processor.process_regular( data = data, coordinates = coordinates )


class Alternate(object):
  """
  An object assigned with an altloc identifier
  """

  def __init__(self, identifier):

    self.identifier = identifier


  def __call__(self, data, coordinates, processor):

    processor.process_altloc(
      data = data,
      coordinates = coordinates,
      identifier = self.identifier,
      )


def altid_for(atom):

  parent = atom.parent()

  if parent:
    return parent.altloc

  else:
    return Empty

# Indexing with altloc support
class Description(object):
  """
  An internal format to allow processing with the visitor pattern
  """

  def __init__(self, data, coordinates, altid):

    self.data = data
    self.coordinates = coordinates

    if altid:
      self.strategy = Alternate( identifier = altid )

    else:
      self.strategy = Empty


  def accept(self, processor):

    self.strategy(
      data = self.data,
      coordinates = self.coordinates,
      processor = processor,
      )


class Indexer(object):
  """
  Indexer that takes into account altloc
  """

  def __init__(self, factory):

    self.factory = factory
    self.regular = self.factory()

    self.altlocs = {}


  def add(self, altloc):

    self.altlocs[ altloc ] = self.factory()


# Visitor
class Inserter(object):
  """
  Fills up an indexer with data
  """

  def __init__(self, indexer):

    self.indexer = indexer


  def process_regular(self, data, coordinates):

    self.indexer.regular.add( object = data, position = coordinates )


  def process_altloc(self, data, coordinates, identifier):

    if identifier not in self.indexer.altlocs:
      self.indexer.add( altloc = identifier )

    self.indexer.altlocs[ identifier ].add( object = data, position = coordinates )


class Aggregator(object):
  """
  Queries the indexer and returns altloc-correct neighbours
  """

  def __init__(self, indexer):

    self.indexer = indexer
    self.ranges = []


  @property
  def entities(self):

    from itertools import chain
    return chain.from_iterable( self.ranges )


  def process_regular(self, data, coordinates):

    self.ranges.append( self.indexer.regular.close_to( centre = coordinates ) )

    for indexer in self.indexer.altlocs.values():
      self.ranges.append( indexer.close_to( centre = coordinates ) )


  def process_altloc(self, data, coordinates, identifier):

    self.ranges.append( self.indexer.regular.close_to( centre = coordinates ) )

    if identifier in self.indexer.altlocs:
      self.ranges.append(
        self.indexer.altlocs[ identifier ].close_to( centre = coordinates ),
        )
