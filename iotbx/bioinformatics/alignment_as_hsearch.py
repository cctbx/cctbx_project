"""Generate information for homology search from a sequence file"""

from __future__ import absolute_import, division, print_function

from iotbx import bioinformatics
from six.moves import zip

def get_extension(filename):

  import os.path
  return os.path.splitext( filename )[1]


def any_alignment(source):

  if isinstance( source, str ):
    infile = open( source )
    data = infile.read()
    infile.close()
    extension = get_extension( filename = source )

  else:
    data = source.read()

    try:
      name = source.name

    except AttributeError:
      extension = None

    else:
      extension = get_extension( filename = name )

  return bioinformatics.any_alignment_string( data = data, extension = extension )


class Result(object):
  """
  Homology search from alignment file
  """

  def __init__(self, alignment):

    self.alignment = alignment
    self.max_count = None


  def restrict(self, max_count):

    if max_count is None:
      self.max_count = None

    else:
      self.max_count = max_count + 1


  def hits(self):

    names = self.alignment.names[ 1 : self.max_count ]
    alignments = self.alignment.alignments[ 1 : self.max_count ]

    try:
      descriptions = self.alignment.descriptions[ 1 : self.max_count ]

    except AttributeError:
      descriptions = alignments

    base = self.alignment.alignments[0]

    for ( n, seq, desc ) in zip( names, alignments, descriptions ):
      pieces = n.split( "_" )
      assert 0 < len( pieces )

      if len( pieces ) == 1:
        pdb = pieces[0]
        chain = ""

      else:
        ( pdb, chain ) = pieces[:2]

      alignment = bioinformatics.clustal_alignment(
        names = [ "target", "%s_%s" % ( pdb, chain ) ],
        alignments = [ base, seq ],
        program = "<unknown>"
        )
      yield bioinformatics.homology_search_hit(
        identifier = pdb,
        chain = chain,
        annotation = desc,
        alignment = alignment
        )


  def __len__(self):

    if self.max_count is None:
      return self.alignment.multiplicity() - 1

    else:
      return self.max_count - 1


def parse(source):

  alignment = any_alignment( source = source )
  return Result( alignment = alignment )

