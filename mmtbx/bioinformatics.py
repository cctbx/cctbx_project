import re
import operator

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
    if any( [ len( alignments[0] ) != len( o ) for o in alignments[1:] ] ):
      raise ValueError, "'alignments' do not have the same length"

    self.names = names
    self.alignments = alignments
    self.gap = gap


  def identity_count(self):

    identical = [ equi for equi in zip( *self.alignments )
      if all( [ p == equi[0] for p in equi[1:] ] ) ]

    return len( [ id for id in identical if self.gap not in id ] )


  def midline(self):

    return midline().compare( self.alignments, self.gap )


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

  def __init__(self, alignments, names, program, version, gap = "-"):

    super( clustal_alignment, self ).__init__( alignments, names, gap )
    if program: self.program = program
    else:       self.program = ""
    self.version = version


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

    if (self.program): program = self.program + " "
    else:              program = ""
    def fmt_num(num):
      if num: return " %s" % num
      return ""
    return (
      "CLUSTAL %(program)s%(version)s multiple sequence alignment\n\n" % {
        "program": program,
        "version": self.version,
        }
      + "\n\n".join(
        [
          "\n".join(
            [ "%s %s%s" % ( cap, ali, fmt_num(num) )
              for ( cap, ali, num ) in zipped_infos ] )
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
        ( unknown, consumed, text ) = text.partition( match.group( 0 ) )

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


# Alignment file parsers
class generic_alignment_parser(object):
  """
  General purpose alignment parser
  """

  def valid_alignments(self, alignments):

    return all(
      [ len( line1 ) == len( line2 )
        for line1 in alignments
        for line2 in alignments[1:] ]
      )


  def fail(self, text):

    return ( None, text )


  def extract(self, text):

    return [ match.groupdict() for match in self.data.finditer( text ) ]


  def __call__(self, text):

    return self.parse( text )


class clustal_alignment_parser(generic_alignment_parser):
  """
  Specific for Clustal alignments
  """

  header = re.compile(
    r"""
    ^ CLUSTAL
    (?: \s (?P<program> [A-Z] ) )? \s
    ( \( )? (?P<version> [\d.]* ) (?(2) \) ) \s
    multiple \s sequence \s alignment \s* \n
    """,
    re.VERBOSE
    )

  data = re.compile(
    r"""
    ^
    (?P<name> \w+ ) \s+
    (?P<alignment> [A-Z\-]* )
    (?P<number> \s+ \d+ )? \s* \n
    """,
    re.VERBOSE | re.MULTILINE
    )

  def parse(self, text):

    # Header must match
    header_match = self.header.search( text )

    if not header_match:
      return self.fail( text )

    # Get names and data
    data = self.extract( text )

    # Sort data on names
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

    # Create alignment object
    header = header_match.groupdict()

    if header[ "program" ]: program = header[ "program" ]
    else:                   program = ""
    return (
      clustal_alignment(
        names = data_for.keys(),
        alignments = alignments,
        program = program,
        version = header[ "version" ]
        ),
      ""
      )

clustal_alignment_parse = clustal_alignment_parser()


class pir_alignment_parser(generic_alignment_parser):
  """
  Specific for pir-format alignments
  """

  data = re.compile(
    r"""
    ^ >
    (?P<type> [PFDRN][13LC] ) ;
    (?P<name> \S+ ) \n
    (?P<description> [^\n]* ) \n
    (?P<alignment> [^>]* )
    \* \s*
    """,
    re.MULTILINE | re.VERBOSE
    )

  def parse(self, text):

    data = self.extract( text )

    if text and not data:
      return self.fail( text )

    # Check that alignments match
    alignments = [
      "".join( [ char for char in info[ "alignment" ] if not char.isspace() ] )
      for info in data
      ]

    if not self.valid_alignments( alignments ):
      return self.fail( text )

    # Create alignment object
    return (
      pir_alignment(
        names = [ info[ "name" ] for info in data ],
        alignments = alignments,
        types = [ info[ "type" ] for info in data ],
        descriptions = [ info[ "description" ] for info in data ]
        ),
      ""
      )

pir_alignment_parse = pir_alignment_parser()
