import re
import operator
import string
import os.path

# Wrap lines that are longer than 'width'
def wrap(text, width):

  return re.findall( "[^\n]{1,%d}" % width, text )


# Sequence formats
class sequence(object):
  """
  Sequence
  """

  def __init__(self, sequence, name):

    self.name = name
    self.sequence = "".join(
      [ char for char in sequence if not char.isspace() ]
      )


  def __eq__(self, other):

    return isinstance( other, sequence ) and self.sequence == other.sequence


  def __ne__(self, other): # for Python 2.2 compatibility

    return not self.__eq__(other)


  def __len__(self):

    return len( self.sequence )


class fasta_sequence(sequence):
  """
  Fasta sequence
  """

  def __init__(self, sequence, name, description):

    super( fasta_sequence, self ).__init__( sequence, name )
    self.description = description


  def format(self, width):

    return "\n".join(
      [ ( ">%s %s" % ( self.name, self.description ) )[:width] ]
      + wrap( self.sequence, width )
      )


  def __str__(self):

    return self.format( 70 )


class pir_sequence(sequence):
  """
  Pir sequence
  """

  def __init__(self, sequence, name, type, description):

    super( pir_sequence, self ).__init__( sequence, name )
    self.type = type
    self.description = description


  def format(self, width):

    return (
      ( ">%s;%s" % ( self.type, self.name ) )[:width] + "\n"
        + self.description[:width] + "\n"
        + "\n".join(
          [ "  " + line for line in wrap( self.sequence, width - 2 ) ]
          )
          + "*"
      )


  def __str__(self):

    return self.format( 70 )


# Alignment middle/match line calculation
class midline(object):
  """
  Alignment midline
  """

  def __init__(self, identical = "*", conserved = ":", semi = ".", differ = " "):

    self.identical = identical
    self.conserved = conserved
    self.semi = semi
    self.differ = differ


  def compare(self, alignments, gap = "-"):

    result = []
    if (len(alignments) != 0):
      for equi in zip( *alignments ):
        if gap not in equi:
          result.append(self.conservation_code( equi ))
        else:
          result.append(self.differ)
    return "".join(result)


  def conservation_code(self, residues):

    for r in residues:
      if (residues[0] != r): return self.differ
    return self.identical


# Alignment formats
class alignment(object):
  """
  Alignment
  """

  def __init__(self, alignments, names, gap = "-"):

    # The number of names should match the number of alignments
    if len( alignments ) != len( names ):
      raise ValueError, (
        "'alignments' and 'names' do not have the same length"
      )

    # All alignments should have the same length
    for o in alignments[1:]:
      if (len( alignments[0] ) != len( o )):
        raise ValueError, "'alignments' do not have the same length"

    self.names = names
    self.alignments = alignments
    self.gap = gap


  def identity_count(self):

    result = 0
    if (len(self.alignments) != 0):
      for equi in zip( *self.alignments ):
        if (self.gap in equi): continue
        for p in equi[1:]:
          if (p != equi[0]): break
        else:
          result += 1
    return result


  def identity_fraction(self):

    if not self.alignments:
      return 1.0

    alignment_length = len( self.alignments[0] )
    shortest_sequence_length = min(
      [ alignment_length - seq.count( self.gap ) for seq in self.alignments ]
      )

    if shortest_sequence_length == 0:
      return 1.0

    return float( self.identity_count() ) / shortest_sequence_length


  def midline(self):

    return midline().compare( self.alignments, self.gap )


  def sequence_strings(self):

    return [ "".join( [ c for c in seq if c != self.gap ] )
        for seq in self.alignments ]


  def length(self):

    if self.alignments:
      return len( self.alignments[0] )
    return 0


  def multiplicity(self):

    return len( self.alignments )


class fasta_alignment(alignment):
  """
  Fasta alignment
  """

  def __init__(self, alignments, names, descriptions, gap = "-"):

    # The number of names, types and description should be equal
    if len( names ) != len( descriptions ):
      raise ValueError, (
        "Inconsistent 'alignments' and 'descriptions' attributes"
        )

    super( fasta_alignment, self ).__init__( alignments, names, gap )
    self.descriptions = descriptions


  def format(self, width):

    return "\n\n".join( [
      fasta_sequence( sequence, name, description ).format( width )
      for ( sequence, name, description )
      in zip( self.alignments, self.names, self.descriptions )
      ] )


  def __str__(self):

    return self.format( 70 )


class pir_alignment(alignment):
  """
  Pir alignment
  """

  def __init__(self, alignments, names, types, descriptions, gap = "-"):

    # The number of names, types and description should be equal
    if ( len( names ) != len( types )
      or len( names ) != len( descriptions ) ):
      raise ValueError, (
        "Inconsistent 'alignments', 'types' and 'descriptions' attributes"
        )

    super( pir_alignment, self ).__init__( alignments, names, gap )
    self.types = types
    self.descriptions = descriptions


  def format(self, width):

    return "\n\n".join( [
      pir_sequence( sequence, name, type, description ).format( width )
      for ( sequence, name, type, description )
      in zip( self.alignments, self.names, self.types, self.descriptions )
      ] )


  def __str__(self):

    return self.format( 70 )


class clustal_alignment(alignment):
  """
  Clustal alignment
  """

  def __init__(self, alignments, names, program = "CLUSTAL 2.0.9", gap = "-"):

    super( clustal_alignment, self ).__init__( alignments, names, gap )
    self.program = program

  def make_aln_info(self, caption, alignment, aln_width):

    running_total = 0
    aln_info = []

    for line in wrap( alignment, aln_width ):
      running_total += len( line ) - line.count( self.gap )
      aln_info.append( ( caption, line, running_total ) )

    return aln_info


  def format(self, aln_width, caption_width):

    # All alignments
    aln_infos = [
      self.make_aln_info(
        caption = name.ljust( caption_width )[:caption_width ],
        alignment = alignment,
        aln_width = aln_width
        )
      for ( name, alignment ) in zip( self.names, self.alignments )
      ]

    # Midline
    aln_infos.append(
      [ ( " " * caption_width, line, "" )
        for line in wrap( self.midline(), aln_width ) ]
      )

    def fmt_num():
      if num: return " %s" % num
      return ""
    return (
      "%s multiple sequence alignment\n\n" % self.program
      + "\n\n".join(
        [
          "\n".join(
            [ "%s %s%s" % ( cap, ali, fmt_num() )
              for ( cap, ali, num ) in zipped_infos ]
            )
          for zipped_infos in zip( *aln_infos )
          ]
        )
      )


  def __str__(self):

    return self.format( 60, 15 )


# Sequence file parser
class generic_sequence_parser(object):
  """
  General-purpose sequence parser
  """

  def __init__(self, regex, type):

    self.regex = regex
    self.type = type


  def parse(self, text):

    objects = []
    non_compliant = []

    while True:
      match = self.regex.search( text )

      if match:
        i = text.find(match.group(0))
        assert i >= 0
        unknown, text = text[:i], text[i+len(match.group(0)):]

        if unknown and not unknown.isspace():
          non_compliant.append( unknown )

        objects.append( self.type( **match.groupdict() ) )

      else:
        break

    if text:
      non_compliant.append( text )

    return ( objects, non_compliant )


  def __call__(self, text):

    return self.parse( text )

# Sequence parser instances that can be used as functions
seq_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    ^ > \s*
    (?P<name> [^\n]* ) \n
    (?P<sequence> [^>^*]* )
    \*?\s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = sequence
  )

fasta_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    ^ >
    (?P<name> \S+ ) \s+
    (?P<description> [^\n]* ) \n
    (?P<sequence> [^>]* )
    \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = fasta_sequence
  )

pir_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    ^ >
    (?P<type> [PFDRN][13LC] ) ;
    (?P<name> \S+ ) \n
    (?P<description> [^\n]* ) \n
    (?P<sequence> [^>]* )
    \* \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = pir_sequence
  )


_implemented_sequence_parsers = {
  ".fasta": fasta_sequence_parse,
  ".pir": pir_sequence_parse,
  ".seq": seq_sequence_parse,
  }

def sequence_parser_for(file_name):

  ( name, extension ) = os.path.splitext( file_name )

  return _implemented_sequence_parsers.get( extension )


def known_sequence_formats():

  return _implemented_sequence_parsers.keys()


# Alignment file parsers
class generic_alignment_parser(object):
  """
  General purpose alignment parser
  """

  def fail(self, text):

    return ( None, text )


  def valid_alignments(self, alignments):

    for line2 in alignments[1:]:
      if (len( alignments[0] ) != len( line2 )):
        return False
    return True


  def extract(self, text):

    return [ match.groupdict() for match in self.regex.finditer( text ) ]


  def assess_parsing_results(self, data_dict, text):

    assert "alignments" in data_dict and "names" in data_dict

    # Remove whitespace from alignments
    alignments = [
      "".join( [ char for char in ali_str if not char.isspace() ] )
      for ali_str in data_dict[ "alignments" ]
      ]

    if not self.valid_alignments( alignments ):
      return ( None, text )

    data_dict[ "alignments" ] = alignments

    # Create alignment object
    return ( self.type( **data_dict ), "" )


  def __call__(self, text):

    return self.parse( text )


class sequential_alignment_parser(generic_alignment_parser):
  """
  Specific for sequential format alignments
  """

  def __init__(self, regex, type):

    self.regex = regex
    self.type = type


  def parse(self, text):

    data = self.extract( text )

    if text and not data:
      return self.fail( text )

    preprocessed_data = dict( [ ( name, [ info[ name ] for info in data ] )
      for name in self.regex.groupindex.keys() ] )

    return self.assess_parsing_results(
      data_dict = preprocessed_data,
      text = text
      )


class clustal_alignment_parser(generic_alignment_parser):
  """
  Specific for Clustal alignments
  """

  regex = re.compile(
    r"""
    ^
    (?P<name> [\S]+ ) \s+
    (?P<alignment> [A-Z\-]* )
    (?P<number> \s+ \d+ )? \s* \n
    """,
    re.VERBOSE | re.MULTILINE
    )
  HEADER = re.compile( r"\A(.*) multiple sequence alignment$", re.MULTILINE )

  def parse(self, text):

    match = self.HEADER.search( text )

    if not match:
      return self.fail( text )

    program = match.group( 1 )

    # Get names and data
    data = self.extract( text )

    # Merge data on names
    data_for = dict( [ ( info[ "name" ], [] ) for info in data ] )

    for info in data:
      data_for[ info[ "name" ] ].append( info[ "alignment" ] )

    # Check that alignments match
    alignments = [
      "".join( [ char for char in "".join( segments ) if not char.isspace() ] )
      for segments in data_for.values()
      ]

    if not self.valid_alignments( alignments ):
      return self.fail( text )

    # Original order need to be maintained
    alignment_for = dict( zip( data_for.keys(), alignments ) )
    unique_names = []

    for info in data:
      if info[ "name" ] not in unique_names:
        unique_names.append( info[ "name" ] )

    # Create alignment object
    return (
      clustal_alignment(
        names = unique_names,
        alignments = [ alignment_for[ name ] for name in unique_names ],
        program = program
        ),
      ""
      )


# Alignment parser instances that can be used as functions
clustal_alignment_parse = clustal_alignment_parser()
pir_alignment_parse = sequential_alignment_parser(
  regex = re.compile(
    r"""
    ^ >
    (?P<types> [PFDRN][13LC] ) ;
    (?P<names> \S+ ) \n
    (?P<descriptions> [^\n]* ) \n
    (?P<alignments> [^>]* )
    \* \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = pir_alignment
  )
fasta_alignment_parse = sequential_alignment_parser(
  regex = re.compile(
    r"""
    ^ >
    (?P<names> \S+ ) \s+
    (?P<descriptions> [^\n]* ) \n
    (?P<alignments> [^>]* )
    \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = fasta_alignment
  )
ali_alignment_parse = sequential_alignment_parser(
  regex = re.compile(
    r"""
    ^ > \s*
    (?P<names> [^\n]* ) \n
    (?P<alignments> [^>^*]* )
    \*?\s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = alignment
  )

_implemented_alignment_parsers = {
  ".pir": pir_alignment_parse,
  ".aln": clustal_alignment_parse,
  ".clustal": clustal_alignment_parse,
  ".fasta": fasta_alignment_parse,
  ".ali": ali_alignment_parse,
  ".fa": ali_alignment_parse,
  }

def alignment_parser_for(file_name):

  ( name, extension ) = os.path.splitext( file_name )

  return _implemented_alignment_parsers.get( extension )


def known_alignment_formats():

  return _implemented_alignment_parsers.keys()
