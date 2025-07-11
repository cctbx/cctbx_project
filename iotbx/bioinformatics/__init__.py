"""
Tools for bioinformatics (reading/writing/comparison of sequence files.)
"""
from __future__ import absolute_import, division, print_function
from operator import itemgetter
from libtbx.utils import Abort
from libtbx import group_args
from six.moves import cStringIO as StringIO
import re
import operator
import os.path
import sys
from functools import reduce
from six.moves import range
from six.moves import zip
from libtbx.utils import Sorry

# Wrap lines that are longer than 'width'
def wrap(text, width):
  """Wrap lines that are longer than 'width'"""

  return re.findall( "[^\n]{1,%d}" % width, text )


# Sequence headers
class ebi_pdb_header(object):
  """
  EBI-style parsed header
  """

  regex = re.compile(
    r"""
    ^ > \s* PDB \s* :
    ( [^_]+ ) _
    ( \w* ) \s*
    """,
    re.VERBOSE
    )

  format = ">PDB:XXXX_X"

  def __init__(self, identifier, chain):

    self.identifier = identifier
    self.chain = chain


  def __str__(self):

    return ">PDB:%s_%s" % ( self.identifier, self.chain )


  def parse(cls, data):

    match = cls.regex.search( data )

    if not match:
      raise ValueError( cls.format, str( data ))

    ( identifier, chain ) = match.groups()
    return cls( identifier = identifier.upper(), chain = chain.upper() )

  parse = classmethod( parse )


class generic_sequence(object):
  """
  Sequence
  """

  def __init__(self, header, body):

    self.header = header
    self.body = body


  def format(self, width):

    return "\n".join(
      [ str( self.header )[:width] ] + wrap( self.body, width )
      )


  def reinterpret_header(self, header):

    self.header = header.parse( data = str( self.header ) )


  def __len__(self):

    return len( self.body )


  def __str__(self):

    return self.format( 70 )


# Sequence formats
class sequence(object):
  """
  Sequence
  """

  def __init__(self, sequence, name = ""):

    self.name = name
    self.sequence = "".join(
      [ char for char in sequence if not char.isspace() ]
      )


  def format(self, width):

    return "\n".join(
      [ ( ">%s" % self.name )[:width] ]
      + wrap( self.sequence, width )
      )


  def __str__(self):

    return self.format( 70 )

  def __hash__(self):
    '''
    Return a UID (hash) for the instanciated object. This method is required
    to be implemented if __eq__ is implemented and if the object is immutable.
    Immutable objects can be then used as keys for dictionaries.
    '''

    #NOTE: To my knowledge objects of this class are used as keyword in a dictionary only
    #by Sculptor.

    #NOTE: id(self) is already unique but hash(id(self)) distributes better the hashed values
    #this means that two consecutive id will have very different hash(id(self)) making
    #hashtable searches more robust when using them as keys of dictionaries. Difference in
    #performance are irrelevant.

    return hash(id(self))

  def __eq__(self, other):
    '''
    This method is used any time the equality among two instances of this class must be tested.
    Two sequences are equal if their sequence attribute is equal. Although self and other are
    two different objects. In particular two sequences might have a different name but they are equal
    if they have the same sequence.
    '''

    return isinstance( other, sequence ) and self.sequence == other.sequence


  def __ne__(self, other): # for Python 2.2 compatibility

    return not self.__eq__(other)


  def __len__(self):

    return len( self.sequence )


  def __getitem__(self, index):

    return self.sequence[ index ]


class fasta_sequence(sequence):
  """
  Fasta sequence
  """

  def __init__(self, sequence, name = "", description = ""):

    super( fasta_sequence, self ).__init__( sequence, name )
    self.description = description


  def format(self, width):

    return "\n".join(
      [ ( ">%s %s" % ( self.name, self.description ) )[:width] ]
      + wrap( self.sequence, width )
      )


class pir_sequence(sequence):
  """
  Pir sequence
  """

  def __init__(self, sequence, name = "", type = "P1", description = ""):

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


class pdb_sequence(sequence):
  """
  PDB sequence
  """

  def __init__(self, sequence, name = "", chain = ""):

    super( pdb_sequence, self ).__init__( sequence, name )
    self.chain = chain


  def format(self, width):

    return (
      ( ">PDB:%s_%s" % ( self.name, self.chain ) )[:width] + "\n"
        + "\n".join(
          [ "  " + line for line in wrap( self.sequence, width - 2 ) ]
          )
      )


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
      raise ValueError(
        "'alignments' and 'names' do not have the same length"
      )

    self._set_alignments( alignments = alignments )

    self.names = names
    self.alignments = alignments
    self.gap = gap


  def aligned_positions_count(self):

    result = 0

    for equi in zip( *self.alignments ):
      if self.gap not in equi:
        result += 1

    return result


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



  def shortest_seqlen(self): # returns the number of aligned residues

    if not self.alignments:
      return 0

    alignment_length = len( self.alignments[0] )
    shortest_sequence_length = min(
      [ alignment_length - seq.count( self.gap ) for seq in self.alignments ]
      )
    return shortest_sequence_length


  def identity_fraction(self):

    if not self.alignments:
      return 1.0

    shortest_sequence_length = self.shortest_seqlen()

    if shortest_sequence_length == 0:
      return 1.0

    return float( self.identity_count() ) / shortest_sequence_length


  def sequence_strings(self):

    return [ "".join( [ c for c in seq if c != self.gap ] )
        for seq in self.alignments ]


  def length(self):

    if self.alignments:
      return len( self.alignments[0] )
    return 0


  def multiplicity(self):

    return len( self.alignments )


  def format(self, width):

    return "\n\n".join( [ sequence( seq, name ).format( width )
      for ( seq, name ) in zip( self.alignments, self.names ) ] )


  def simple_format(self, width):

    return self.format( width = width )


  def extend(self, sequences):

    if len( sequences ) != self.multiplicity():
      raise ValueError("Inconsistent extension sequence set")

    pre_alis = []
    pre_gaps = []
    post_alis = []
    post_gaps = []

    for ( a, s ) in zip( self.alignments, sequences ):
      partial = "".join( [ c for c in a if c != self.gap ] )
      pos = s.find( partial )

      if pos == -1:
        raise ValueError("Alignment does not match sequence")

      assert 0 <= pos
      pre_alis.append( s[:pos] )
      pre_gaps.append( self.gap * pos )

      end = pos + len( partial )
      post_alis.append( s[ end : ] )
      post_gaps.append( self.gap * ( len( s ) - end ) )

    alis = []

    for ( index, a ) in enumerate( self.alignments ):
      alis.append(
        reduce(
          operator.add,
          pre_gaps[ : index ] + [ pre_alis[ index ] ] + pre_gaps[ index + 1 : ]
          )
        + a
        + reduce(
          operator.add,
          post_gaps[ : index ] + [ post_alis[ index ] ] + post_gaps[ index + 1 : ]
          )
        )

    self._set_alignments( alignments = alis )


  def copy(self, **kwd):

      keywords = {
          "alignments": self.alignments,
          "names": self.names,
          "gap": self.gap,
          }
      keywords.update( self.extra() )
      keywords.update( kwd )
      return self.__class__( **keywords )


  def extra(self):

      return {}


  def assign_as_target(self, index):

    if index < 0 or self.multiplicity() <= index:
      raise IndexError("Sequence not found in alignment")

    self.names = ( [ self.names[ index ] ] + self.names[ : index ]
      + self.names[ index + 1 : ] )
    self.alignments = ( [ self.alignments[ index ] ]
      + self.alignments[ : index ] + self.alignments[ index + 1 : ] )


  def _set_alignments(self, alignments):

    # All alignments should have the same length
    for o in alignments[1:]:
      if (len( alignments[0] ) != len( o )):
        raise ValueError("'alignments' do not have the same length")

    self.alignments = alignments


  def __str__(self):

    return self.format( 70 )


class fasta_alignment(alignment):
  """
  Fasta alignment
  """

  def __init__(self, alignments, names, descriptions, gap = "-"):

    # The number of names, types and description should be equal
    if len( names ) != len( descriptions ):
      raise ValueError(
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


  def simple_format(self, width):

    return self.format( width = width )


  def extra(self):

    return {
      "descriptions": self.descriptions,
      }


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
      raise ValueError(
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


  def simple_format(self, width):

    return self.format( width = width )


  def extra(self):

    return {
      "descriptions": self.descriptions,
      "types": self.types,
      }


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


  def format(self, aln_width, caption_width, number_width = None, middle_line=None):

    if not middle_line:
      middle_line = midline().compare(
        alignments = self.alignments,
        gap = self.gap
        )

    elif len( middle_line ) != self.length():
      raise ValueError("Incorrect midline length")

    if number_width is None:
      number_width = self.length_digits()

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
        for line in wrap( middle_line, aln_width ) ]
      )

    def fmt_num(num):
      if num: return " " + str( num ).rjust( number_width )
      return ""
    return (
      "%s multiple sequence alignment\n\n" % self.program
      + "\n\n".join(
        [
          "\n".join(
            [ "%s %s%s" % ( cap, ali, fmt_num(num) )
              for ( cap, ali, num ) in zipped_infos ]
            )
          for zipped_infos in zip( *aln_infos )
          ]
        )
      )


  def length_digits(self):

    import math
    return int( math.log10( self.length() ) ) + 1


  def simple_format(self, width):

    aln_width = int( 0.8 * width )
    num_width = min( self.length_digits(), width - aln_width )
    caption_width = max( width - ( aln_width + 1 ) - ( num_width + 1 ), 1 )
    return self.format(
      aln_width = aln_width,
      caption_width = caption_width,
      number_width = num_width,
      )


  def extra(self):

    return {
      "program": self.program,
      }


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


  def parse(self, text, **kwargs):

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

        objects.append(
          self.type( **dict( list(kwargs.items()) + list(match.groupdict().items()) ) )
          )

      else:
        break

    if text:
      non_compliant.append( text )

    return ( objects, non_compliant )


  def __call__(self, text, **kwargs):

    return self.parse( text, **kwargs )

# New-style sequence parsing
class SequenceFormat(object):
  """
  Bundle of independent parsing data
  """

  def __init__(self, regex, expected, create):

    self.regex = regex
    self.expected = expected
    self.create = create


st_separated_fasta = SequenceFormat(
  regex = re.compile(
    r"""
    ^
    ( > [^\n]* ) \n
    ( [^>^*]* )
    \*?\s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  expected = ">Header\nSEQUENCE",
  create = lambda headers, body: generic_sequence(
      header = headers[0],
      body  = body
      )
  )


def partition_string(data, regex):
  """
  Partition sequence files
  """

  position = 0

  for match in regex.finditer( data ):
    current = match.start()
    yield ( match.groups(), ( position, current, data[ position : current ] ) )
    position = match.end()

  remaining = data[ position : ].strip()

  if remaining:
    yield ( (), ( position, len( data ), remaining ) )

  raise StopIteration


def parse_sequence_str(data, format):
  """
  Generic purpose sequence header parser
  """

  partition = partition_string( data = data, regex = format.regex )

  for ( groups, ( n_start, n_end, n_data ) ) in partition:
    if n_data.strip():
      raise ValueError(
        "Uninterpretable block from %s to %s:\nExpected: %s\nFound: %s" % (
          n_start,
          n_end,
          format.expected,
          n_data,
          ))

    body = "".join( [ c for c in groups[-1] if not c.isspace() ] )
    yield format.create( headers = groups[:-1], body = body )


# Sequence parser instances that can be used as functions
seq_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    ^ > [ \t\v\f]*
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
    ^ > [ \t]*
    (?P<name> [^\s^:]+ ) [ :]? [ \t]*
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
    (?P<name> \S* ) \n
    (?P<description> [^\n]* ) \n
    (?P<sequence> [^>^\*]* )
    \* \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = pir_sequence
  )

tolerant_pir_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    ^ >
    (?P<name> [^\n]* ) \n
    (?P<description> [^\n]* ) \n
    (?P<sequence> [^>^\*]* )
    \*? \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = pir_sequence
  )

lineseparated_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    (?P<sequence> [^>^*]*? )
    \*?\s*\n(?:\n|\Z)
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = sequence
  )

dbfetch_sequence_parse = generic_sequence_parser(
  regex = re.compile(
    r"""
    ^ > \s* PDB \s* :
    (?P<name> [^_]+ ) _
    (?P<chain> \w* ) \s* \n
    (?P<sequence> [^>]* )
    \s*
    """,
    re.MULTILINE | re.VERBOSE
    ),
  type = pdb_sequence
  )

_tf_split_regex = re.compile(
  r"""
  \A \s*
  > ( [^\n]* ) \n
  ( .* )
  \Z
  """,
  re.MULTILINE | re.VERBOSE | re.DOTALL
  )

def tf_sequence_parse(text):
  """
  plain text sequence parser used for SOLVE/RESOLVE files
  """

  match = _tf_split_regex.search( text )

  if not match:
    return ( [], [ text ] )

  ( name, data ) = match.groups()

  return lineseparated_sequence_parse.parse(
    text = data,
    name = name.strip()
    )


_implemented_sequence_parsers = {
  ".fasta": fasta_sequence_parse,
  ".fa": fasta_sequence_parse,
  ".faa" : fasta_sequence_parse,
  ".pir": pir_sequence_parse,
  ".seq": seq_sequence_parse,
  ".dat": seq_sequence_parse,
  }

def sequence_parser_for(file_name):
  """Identify sequence parser suitable for file_name"""
  ( name, extension ) = os.path.splitext( file_name )

  return _implemented_sequence_parsers.get( extension )


def sequence_parser_for_extension(extension):
  """Identify sequence parser suitable for extension"""

  return _implemented_sequence_parsers.get( extension )


def known_sequence_formats():
  """List known sequence formats"""
  return list(_implemented_sequence_parsers.keys())


def parse_sequence(data):
  """Parse a sequence"""
  parsers = [
    fasta_sequence_parse,
    dbfetch_sequence_parse,
    pir_sequence_parse,
    tolerant_pir_sequence_parse,
    seq_sequence_parse,
    lineseparated_sequence_parse,
    tf_sequence_parse,
    ]

  results = []

  for p in parsers:
    ( seqs, junk ) = p( text = data )

    if not junk:
      return ( seqs, junk )

    else:
      results.append( ( seqs, junk ) )

  return min( results, key = lambda p: len( p[1] ) )

# XXX test needed
def any_sequence_format(file_name, assign_name_if_not_defined=False,
    data=None):
  """Parse any sequence from any format"""
  format_parser = sequence_parser_for(file_name)
  if data is None:
    with open(file_name, "r") as f:
      data = f.read()
  seq_object = None
  if (format_parser is not None):
    try :
      objects, non_compliant = format_parser.parse(data)
      assert (len(objects) > 0)
      assert (objects[0].sequence != "")
    except Exception as e :
      pass
    else :
      if (assign_name_if_not_defined):
        base_name = os.path.splitext(os.path.basename(file_name))[0]
        for k, seq in enumerate(objects):
          if (seq.name == ""):
            seq.name = "%s_%d" % (base_name, k+1)
            print(seq.name)
      return objects, non_compliant
  for other_parser in [fasta_sequence_parse.parse, pir_sequence_parse.parse,
                       seq_sequence_parse.parse, tf_sequence_parse] :
    if (other_parser is not format_parser):
      try :
        objects, non_compliant = other_parser(data)
        assert (len(objects) > 0)
        assert (objects[0].sequence != "")
      except Exception as e :
        pass
      else :
        if (assign_name_if_not_defined):
          base_name = os.path.splitext(os.path.basename(file_name))[0]
          for k, seq in enumerate(objects):
            if (seq.name == ""):
              seq.name = "%s_%d" % (base_name, k+1)
              print(seq.name)
        return objects, non_compliant
  # fallback: unformatted
  data = re.sub(r"\s", "", data)
  if (re.search(r"[^a-zA-Z\*]", data) is None):
    seq = sequence(data)
    if (assign_name_if_not_defined ):
      seq.name = os.path.splitext(os.path.basename(file_name))[0]
    return [seq], []
  return None, None

def merge_sequence_files(file_names, output_file, sequences=(),
    include_non_compliant=False, call_back_on_error=None,
    force_new_file=False):
  """Merge a set of sequence files into a single file"""
  assert (len(file_names) > 0) or (len(sequences) > 0)
  if (len(file_names)==1) and (len(sequences)==0) and (not force_new_file):
    return file_names[0]
  seq_out = StringIO()
  for seq_file in file_names :
    objects, non_compliant = any_sequence_format(seq_file)
    if (objects is None):
      msg = ("The file '%s' could not be parsed as one of the "+
          "standard formats.  This could either be a problem with the file, "+
          "a non-standard format, or a bug in our code.  For further help, "+
          "please email the developers at help@phenix-online.org.") % seq_file
      if (hasattr(call_back_on_error, "__call__")):
        msg += "  (This is not a fatal error, but the combined sequence "+\
                "file will be incomplete.)"
        if (not call_back_on_error(msg)):
          raise Abort()
        continue
      else :
        raise RuntimeError(msg)
    for seq_record in objects :
      name = str(seq_record.name)
      if (name == ""):
        name = "none"
      seq_out.write("> %s\n" % name)
      seq_out.write("%s\n" % seq_record.sequence)
  if (len(sequences) > 0):
    for k, seq in enumerate(sequences):
      seq_out.write("> seq%d\n" % (k+1))
      seq_out.write("%s\n" % seq)
  f = open(output_file, "w")
  f.write(seq_out.getvalue())
  f.close()
  return output_file

# Alignment file parsers
class generic_alignment_parser(object):
  """
  General purpose alignment parser
  """

  def fail(self, text):

    return ( None, text )


  def extract(self, text):

    return [ match.groupdict() for match in self.regex.finditer( text ) ]


  def assess_parsing_results(self, data_dict, text):

    assert "alignments" in data_dict and "names" in data_dict

    # Remove whitespace from alignments
    alignments = [
      "".join( [ char for char in ali_str if not char.isspace() ] )
      for ali_str in data_dict[ "alignments" ]
      ]

    if not check_alignments_are_valid( alignments ):
      return ( None, text )

    data_dict[ "alignments" ] = alignments

    # Create alignment object
    return ( self.type( **data_dict ), "" )


  def __call__(self, text):

    return self.parse( text )


def check_alignments_are_valid(alignments):
  """Check for valid alignments"""

  if not alignments:
    return True

  first = len( alignments[0] )

  for line in alignments[1:]:
    if first != len( line ):
      return False

  return True


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


CLUSTAL_BODY = re.compile(
  r"""
  ^
  (?P<name> .+? ) \s+
  (?P<alignment> [A-Z\-]+ )
  (?P<number> \s+ \d+ )? \s* $
  """,
  re.VERBOSE
  )
CLUSTAL_MIDLINE = re.compile( r"^ \s* [:.* ]+ \s* $", re.VERBOSE )
CLUSTAL_HEADER = re.compile( r"\A(.*) multiple sequence alignment$" )

def clustal_alignment_parse(text):
  """
  Specific for Clustal alignments
  """

  lines = list( reversed( text.splitlines() ) ) # create a stack

  if not lines:
    return ( None, text )

  assert 0 < len( lines )

  match = CLUSTAL_HEADER.search( lines.pop() )

  if not match:
    return ( None, text )

  program = match.group( 1 )

  # Read first block
  discard_clustal_empty_lines( lines )

  try:
    ( names, sequences ) = read_clustal_block( lines )

  except ValueError:
    return ( None, text )

  assert len( names ) == len( sequences )
  parts = [ sequences ]
  discard_clustal_empty_lines( lines )

  while lines:
    try:
      ( ns, sequences ) = read_clustal_block( lines )

    except ValueError:
      return ( None, text )

    assert len( ns ) == len( sequences )

    if names != ns:
      return ( None, text )

    assert len( names ) == len( sequences )

    parts.append( sequences )
    discard_clustal_empty_lines( lines )

  return (
    clustal_alignment(
      names = names,
      alignments = [
        "".join( seqs[ i ] for seqs in parts ) for i in range( len( names ) )
        ],
      program = program
      ),
    ""
    )


def read_clustal_block(lines):
  """Read a block from a clustal-format file"""

  names = []
  sequences = []

  while lines:
    if not lines[-1].strip():
      break

    # Get names and data
    match = CLUSTAL_BODY.search( lines[-1] )

    if match:
      info = match.groupdict()
      names.append( info[ "name" ] )
      sequences.append( info[ "alignment" ] )

    else:
      if not CLUSTAL_MIDLINE.search( lines[-1] ):
        raise ValueError("Line '%s' does not match expected body or midline formats" % lines[-1])

    lines.pop()

  return ( names, sequences )


def discard_clustal_empty_lines(lines):
  """Discard empty lines from clustal file"""
  while lines:
    if lines[ -1 ].strip():
      break

    lines.pop()


def hhalign_alignment_parse(text):
  """Parse an HHalign formatted file"""
  try:
    hhalign = hhalign_parser( output = text )

  except ValueError:
    return ( None, text )

  return ( hhalign.alignment(), "" )


# Alignment parser instances that can be used as functions
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
  ".hhr": hhalign_alignment_parse,
  }

def alignment_parser_for(file_name):
  """Identify alignment parser for file_name"""

  ( name, extension ) = os.path.splitext( file_name )

  return _implemented_alignment_parsers.get( extension )


def alignment_parser_for_extension(extension):
  """Identify alignment parser for extension"""

  return _implemented_alignment_parsers.get( extension )


def known_alignment_formats():
  """List known alignment formats"""
  return list(_implemented_alignment_parsers.keys())

def any_alignment_file(file_name):
  """Parse any alignment file"""
  base, ext = os.path.splitext(file_name)
  with open(file_name) as f:
    data = f.read()
  parser1 = None
  if (ext != ".hhr") and (ext in _implemented_alignment_parsers):
    parser1 = _implemented_alignment_parsers.get(ext)
    try :
      aln = parser1(data)[0]
      assert (len(aln.names) != 0)
    except KeyboardInterrupt : raise
    except Exception as e : pass
    else :
      return aln
  for parser in [pir_alignment_parse, clustal_alignment_parse,
                 fasta_alignment_parse, ali_alignment_parse] :
    if (parser is parser1):
      continue
    try :
      aln = parser(data)[0]
      assert (len(aln.names) != 0)
    except KeyboardInterrupt : raise
    except Exception as e : pass
    else :
      return aln
  raise RuntimeError("Not a recognizeable sequence alignment format!")


def any_alignment_string(data, extension = None):
  """Parse any alignment string"""

  tried = None

  if extension == ".hhr":
    return hhalign_parser( output = data )

  elif extension is not None:
    parser = alignment_parser_for_extension( extension = extension )
    aln = parser( text = data )[0]
    tried = parser

    if aln and aln.multiplicity() != 0:
      return aln

  for parser in [pir_alignment_parse, clustal_alignment_parse,
                 fasta_alignment_parse, ali_alignment_parse] :
    if parser is tried:
      continue

    aln = parser( text = data)[0]

    if aln and aln.multiplicity() != 0:
      return aln

  raise RuntimeError("Not a recognizeable sequence alignment format!")


class homology_search_hit(object):
  """
  A hit from a homology search
  """

  def __init__(self, identifier, chain, annotation, alignment):

    self.identifier = identifier
    self.chain = chain
    self.annotation = annotation
    self.alignment = alignment


  def target_alignment_index(self):

    return 0


  def model_alignment_index(self):

    return 1


  def target_alignment_sequence(self):

    return self.alignment.alignments[ self.target_alignment_index() ]


  def model_alignment_sequence(self):

    return self.alignment.alignments[ self.model_alignment_index() ]


class hhpred_parser(object):
  """
  Parses .hhr files from HHPred
  """

  SPLIT = re.compile(
    r"""
    ( Query .*? ) ( \n | \r\n | \r ) \2
    ( \s+ No \s+ Hit .*? ) \2 \2
    ( .*? )
    (?= (?:Done!) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL
    )

  HEADER = re.compile(
    r"""
    \A
    Query \s+ ( \S .* ) (?: \n | \r\n | \r )
    Match_columns \s+ ( \d+ ) (?: \n | \r\n | \r )
    No_of_seqs \s+ ( \d + ) \s out \s of \s ( \d+ ) (?: \n | \r\n | \r )
    Neff \s+ (\S[^\n^\r]*) (?: \n | \r\n | \r )
    Searched_HMMs \s+ ( \d+ ) (?: \n | \r\n | \r )
    Date \s+ (\S[^\n^\r]*) (?: \n | \r\n | \r )
    Command \s+ (\S.*)
    \Z
    """,
    re.VERBOSE
    )


  def __init__(self, output):

    res = self.SPLIT.search( output )

    if not res:
      raise ValueError("Incorrect file format")

    ( header, linefeed, summary, hits ) = res.groups()
    self.process_header( header = header )
    self.process_hits( hits = hits )


  def process_header(self, header):

    res = self.HEADER.search( header )

    if not res:
      raise ValueError("Incorrect header format")

    self.query = res.group( 1 )
    self.match_columns = int( res.group( 2 ) )
    self.no_of_sequences = ( int( res.group( 3 ) ), int( res.group( 4 ) ) )
    self.neff = res.group( 5 )
    self.searched_hmms = int( res.group( 6 ) )
    self.date = res.group( 7 )
    self.command = res.group( 8 )


  def process_hits(self, hits):

    self.setup_data_arrays()
    start = 0

    for m in self.HITS.finditer( hits ):
      if m.start() != start:
        assert start < m.start()

        unknown = hits[ start : m.start() ].strip()

        if unknown:
          raise ValueError("Uninterpretable block: %s" % repr( unknown ))

      start = m.end()

      self.add_match_to_hit_header_results( match = m )
      self.process_blocks( blocks = m.groupdict()[ "blocks" ] )

    remaining = hits[ start: ].strip()

    if remaining:
      raise ValueError("Uninterpretable tail: %s" % repr( remaining ))


  def process_blocks(self, blocks):

    matches = []
    start = 0

    for match in self.BLOCKS.finditer( blocks ):
      if match.start() != start:
        assert start < match.start()

        unknown = blocks[ start : match.start() ].strip()

        if unknown:
          raise ValueError("Uninterpretable: %s" % repr( unknown ))

      start = match.end()
      matches.append( match.groups() )

    remaining = blocks[ start: ].strip()

    if remaining:
      raise ValueError("Uninterpretable: %s" % repr( remaining ))

    if not matches:
      raise ValueError("Empty homology block")

    assert all([len( matches[0] ) == len( a ) for a in matches[1:]])

    return self.merge_and_process_block_hits( matches = matches )


  def merge_sequence_numbers(self, starts, ends, others):

    for ( s, e ) in zip( starts[1:], ends[:-1] ):
      if s != e + 1:
        raise ValueError("Incorrect sequence indices")

      if not all([others[0] == o for o in others[1:]]):
        raise ValueError("Incorrect sequence indices")

    return ( starts[0], ends[-1], others[0] )


class hhsearch_parser(hhpred_parser):
  """
  Specific for hhsearch output
  """

  HITS = re.compile(
    r"""
    No \s ( \d+) \s* (?: \n | \r\n | \r )
    >( [\w]{4} ) _  ( [\w]+ ) \s ( [^\n]* )(?: \n | \r\n | \r )
    Probab = ( [+-]? \d+ \. \d* ) \s+
    E-value = ( \d+ \.? \d* )( e[+-]? \d+ )? \s+
    Score = ( [+-]? \d+\.\d+ ) \s+
    (?: Aligned_cols | Aligned_columns ) = ( \d+ ) \s+
    Identities = ( \d+ ) %
    [^\n]*
    (?: \n | \r\n | \r )
    (?P<blocks> .*? )(?= (?:^No) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL | re.MULTILINE
    )
  BLOCKS = re.compile(
    r"""
    (?: Q .* (?: \n | \r\n | \r ) )*
    Q \s+ [\w:\.|-]* \s+ ( \d+ ) \s+ ( [\w\.-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* (?: \n | \r\n | \r )
    Q \s+ Consensus \s+ \d+ \s+ [\w\.~-]+ \s+ \d+ \s+ \( \d+ \) \s* (?: \n | \r\n | \r )
    \s* [ \.\-+|=]* (?: \n | \r\n | \r )
    T \s+ Consensus \s+ \d+ \s+ [\w\.~-]+ \s+ \d+ \s+ \( \d+ \) \s* (?: \n | \r\n | \r )
    T \s+ \w+ \s+ ( \d+ ) \s+ ( [\w\.-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* (?: \n | \r\n | \r )
    (?: T .* (?: \n | \r\n | \r ) )*
    (?: Confidence .* (?: \n | \r\n | \r ))?
    """,
    re.VERBOSE
    )

  def setup_data_arrays(self):

    self.indices = []
    self.pdbs = []
    self.chains = []
    self.annotations = []
    self.probabs = []
    self.e_values = []
    self.scores = []
    self.aligned_cols = []
    self.identities = []

    self.query_starts = []
    self.query_ends = []
    self.query_others = []
    self.query_alignments = []

    self.hit_starts = []
    self.hit_ends = []
    self.hit_others = []
    self.hit_alignments = []


  def restrict(self, max_count):

    self.indices = self.indices[:max_count]
    self.pdbs = self.pdbs[:max_count]
    self.chains = self.chains[:max_count]
    self.annotations = self.annotations[:max_count]
    self.probabs = self.probabs[:max_count]
    self.e_values = self.e_values[:max_count]
    self.scores = self.scores[:max_count]
    self.aligned_cols = self.aligned_cols[:max_count]
    self.identities = self.identities[:max_count]

    self.query_starts = self.query_starts[:max_count]
    self.query_ends = self.query_ends[:max_count]
    self.query_others = self.query_others[:max_count]
    self.query_alignments = self.query_alignments[:max_count]

    self.hit_starts = self.hit_starts[:max_count]
    self.hit_ends = self.hit_ends[:max_count]
    self.hit_others = self.hit_others[:max_count]
    self.hit_alignments = self.hit_alignments[:max_count]


  def add_match_to_hit_header_results(self, match):

    self.indices.append( int( match.group( 1 ) ) )
    self.pdbs.append( match.group( 2 ).upper() )
    self.chains.append( match.group( 3 ) )
    self.annotations.append( match.group( 4 ) )
    self.probabs.append( float( match.group( 5 ) ) )
    number = match.group( 6 )
    if match.group( 7 ): number += match.group( 7 )
    self.e_values.append( float( number ) )
    self.scores.append( float( match.group( 8 ) ) )
    self.aligned_cols.append( int( match.group( 9 ) ) )
    self.identities.append( float( match.group( 10 ) ) )


  def merge_and_process_block_hits(self, matches):

    data = list(zip( *matches ))
    assert len( data ) == 8

    sequences = [ data[1], data[5] ]

    for index in range( len( matches ) ):
      if len( sequences[0][ index ] ) != len( sequences[1][ index ] ):
        raise ValueError("Incorrect alignments")

    mergeds = [ "".join( a ) for a in sequences ]

    q_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[0] ],
      ends = [ int( d ) for d in data[2] ],
      others = [ int( d ) for d in data[3] ]
      )

    t_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[4] ],
      ends = [ int( d ) for d in data[6] ],
      others = [ int( d ) for d in data[7] ]
      )

    self.query_starts.append( q_indices[0] )
    self.query_ends.append( q_indices[1] )
    self.query_others.append( q_indices[2] )
    self.query_alignments.append( mergeds[0] )

    self.hit_starts.append( t_indices[0] )
    self.hit_ends.append( t_indices[1] )
    self.hit_others.append( t_indices[2] )
    self.hit_alignments.append( mergeds[1] )


  def hits(self):

    data = list(zip(
      self.pdbs,
      self.chains,
      self.annotations,
      self.query_alignments,
      self.hit_alignments,
      ))

    for ( pdb, chain, annotation, query, hit ) in data:
      alignment = clustal_alignment(
        names = [ "target", "%s_%s" % ( pdb, chain ) ],
        alignments = [ query, hit ],
        program = "HHPred"
        )

      yield homology_search_hit(
        identifier = pdb,
        chain = chain,
        annotation = annotation,
        alignment = alignment
        )


  def __len__(self):

    return len( self.pdbs )


class hhalign_parser(hhpred_parser):
  """
  Specific for hhalign output
  """

  HITS = re.compile(
    r"""
    No \s ( \d+) \s* (?: \n | \r\n | \r )
    > \s* ( [^\n]* ) (?: \n | \r\n | \r )
    Probab = ( [+-]? \d+ \. \d* ) \s+
    E-value = ( \d+ \.? \d* )( e[+-]? \d+ )? \s+
    Score = ( \d+\.\d+ ) \s+
    (?: Aligned_cols | Aligned_columns ) = ( \d+ ) \s+
    Identities = ( \d+ ) %
    [^\n]*
    (?: \n | \r\n | \r )
    (?P<blocks> .*? )(?= (?:^No) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL | re.MULTILINE
    )
  BLOCKS = re.compile(
    r"""
    Q \s+ ss_pred \s+ ( [\w-]+ ) \s* (?: \n | \r\n | \r )
    Q \s+ ss_conf \s+ ( [\w-]+ ) \s* (?: \n | \r\n | \r )
    Q \s+ [\w:\.]+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* (?: \n | \r\n | \r )
    Q \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* (?: \n | \r\n | \r )
    \s+ ( [ \.\-+|=]+ ) (?: \n | \r\n | \r )
    T \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* (?: \n | \r\n | \r )
    T \s+ [\w\.]+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* (?: \n | \r\n | \r )
    T \s+ ss_pred \s+ ( [\w-]+ ) \s* (?: \n | \r\n | \r )
    T \s+ ss_conf \s+ ( [\w-]+ ) \s* (?: \n | \r\n | \r )
    """,
    re.VERBOSE
    )

  def setup_data_arrays(self):

    self.indices = []
    self.annotations = []
    self.probabs = []
    self.e_values = []
    self.scores = []
    self.aligned_cols = []
    self.identities = []

    self.query_starts = []
    self.query_ends = []
    self.query_others = []
    self.query_alignments = []
    self.query_consensi = []
    self.query_ss_preds = []
    self.query_ss_confs = []

    self.midlines = []

    self.hit_starts = []
    self.hit_ends = []
    self.hit_others = []
    self.hit_alignments = []
    self.hit_consensi = []
    self.hit_ss_preds = []
    self.hit_ss_confs = []


  def add_match_to_hit_header_results(self, match):

    self.indices.append( int( match.group( 1 ) ) )
    self.annotations.append( match.group( 2 ) )
    self.probabs.append( float( match.group( 3 ) ) )
    number = match.group( 4 )
    if match.group( 5 ):
      number += match.group( 5 )
    self.e_values.append( float( number ) )
    self.scores.append( float( match.group( 6 ) ) )
    self.aligned_cols.append( int( match.group( 7 ) ) )
    self.identities.append( float( match.group( 8 ) ) )


  def merge_and_process_block_hits(self, matches):

    data = list(zip( *matches ))
    assert len( data ) == 21

    sequences = [
        data[0], data[1], data[3], data[7],
        data[12], data[16], data[19], data[20] ]
    midlines = []

    for ( index , alis) in enumerate( zip( *sequences ) ):
      count = len( alis[0] )

      if not all([count == len( c ) for c in alis[1:]]):
        raise ValueError("Incorrect alignments")

      midlines.append(
        " " * ( count - len( data[10][ index ] ) ) + data[10][ index ]
        )

    merged = [ reduce( operator.add, a ) for a in sequences ]
    assert len( midlines ) == len( matches )
    midline = reduce( operator.add, midlines )

    # Comment out consistency check
    # if data[2] != data[6] or data[4] != data[8] or data[5] != data[9]:
    #  raise ValueError, "Inconsistent query numbering"

    q_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[2] ],
      ends = [ int( d ) for d in data[4] ],
      others = [ int( d ) for d in data[5] ]
      )

    # Comment out consistency check
    # if data[11] != data[15] or data[13] != data[17] or data[14] != data[18]:
    #  raise ValueError, "Inconsistent target numbering"

    t_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[11] ],
      ends = [ int( d ) for d in data[13] ],
      others = [ int( d ) for d in data[18] ]
      )

    self.query_starts.append( q_indices[0] )
    self.query_ends.append( q_indices[1] )
    self.query_others.append( q_indices[2] )
    self.query_ss_preds.append( merged[0] )
    self.query_ss_confs.append( merged[1] )
    self.query_alignments.append( merged[2] )
    self.query_consensi.append( merged[3] )

    self.midlines.append( midline )

    self.hit_starts.append( t_indices[0] )
    self.hit_ends.append( t_indices[1] )
    self.hit_others.append( t_indices[2] )
    self.hit_consensi.append( merged[4] )
    self.hit_alignments.append( merged[5] )
    self.hit_ss_preds.append( merged[6] )
    self.hit_ss_confs.append( merged[7] )


  def alignment(self):

    assert len( self.annotations ) == 1
    return clustal_alignment(
      names = [ "target", self.annotations[0] ],
      alignments = [ self.query_alignments[0], self.hit_alignments[0] ],
      program = "HHPred"
      )

def any_hh_file(file_name):
  """Parse any HHalign file"""
  with open(file_name) as f:
    data = f.read()
  for parser in [hhalign_parser, hhsearch_parser] :
    try :
      p = parser(data)
    except KeyboardInterrupt : raise
    except Exception as e :
      pass
    else :
      return p
  raise RuntimeError("Not an HHpred/HHalign/HHsearch file!")

def any_a3m_file(file_name):
  """Parse any a3m file"""
  with open(file_name) as f:
    data = f.read()
  try:
    a3m_info = read_a3m(data)
  except KeyboardInterrupt : raise
  except Exception as e :
    pass
  else :
    return a3m_info
  raise RuntimeError("Not an a3m file")

def read_a3m(text):
  """ Read text as a3m format
      Ignore any lines starting with # at top
      Lines starting with > or blank lines are separators and are ignored
      Text between separators are sequence info.
        First line is base sequence.  It cannot contain any lowercase or "-"
        or "."
        All other lines contain upper case characters (sequence matches),
        lower case characters (insertions), and "-" or "." (gaps).
      The total number of characters in base sequence must equal the
        number of upper case characters plus number of gap characters, minus
        the number of insertion characters in each other line.
  """

  # Get text as utf-8
  if hasattr(text,'decode'):
    text = text.decode(encoding='utf-8')

  # Read in line representing each sequence
  sequences = []
  new_line = ""
  for line in text.splitlines():
    line = line.strip()
    if not sequences and line.startswith("#"): continue # skip leading # lines
    if not line or line.startswith(">"):  # separator
      if new_line:
        sequences.append(new_line)
        new_line = ""
      continue
    new_line += line.strip().replace(" ","")
  if new_line:
    sequences.append(new_line)


  # Check for illegal characters and length of lines
  base_sequence = sequences[0]
  n = len(base_sequence)
  import re
  for s in sequences:
    if not ok_a3m_sequence(s, n = n, base_sequence = base_sequence):
      return None

  from libtbx import group_args
  a3m_info = group_args(group_args_type = 'a3m_info',
    base_sequence = base_sequence,
    sequence_length = len(base_sequence),
    sequences = sequences,
     )
  return a3m_info

def ok_a3m_sequence(s, n = None, base_sequence = None):
  """ Check a sequence and make sure it matches expectations for an a3m line
  """
  # Remove blanks/linefeeds and convert . to -
  s = s.replace(" ","").replace("\r","").replace("\n","")
  s = s.replace(".","-")

  n_gap_chars = s.count("-")
  s_all = s

  # Get upper and lowercase
  s = s.replace("-","")
  # count lowercase/uppercase
  n_upper = 0
  n_lower = 0
  n_other = 0
  s_upper = ""
  for c in s:
    if c >="A" and c <= "Z":
      n_upper += 1
      s_upper += c
    elif c >="a" and c <= "z":
      n_lower += 1
    else:
      n_other += 1

  if n_other > 0:
    return False
  if n_upper + n_gap_chars == n:
    return True
  else:
    return False
def composition_from_sequence_file(file_name, log=None):
  """Get composition from a sequence file"""
  if (log is None):
    log = sys.stdout
  try :
    seq_file, non_compliant = any_sequence_format(file_name)
    if (seq_file is None):
      raise ValueError("Could not parse %s" % file_name)
  except Exception as e :
    print(str(e), file=log)
    return None
  else :
    if (len(non_compliant) > 0):
      print("Warning: non-compliant entries in sequence file", file=log)
      for nc in non_compliant :
        print("  " + str(nc), file=log)
    n_residues = n_bases = 0
    for seq_entry in seq_file :
      (n_res_seq, n_base_seq) = composition_from_sequence(seq_entry.sequence)
      n_residues += n_res_seq
      n_bases += n_base_seq
    return group_args(n_residues=n_residues, n_bases=n_bases)

def composition_from_sequence(sequence):
  """Get composition from a sequence"""
  seq = sequence.upper()
  n_residues = n_bases = 0
  n_valid = len(seq) - seq.count("X")
  n_na = seq.count("A") + seq.count("U") + seq.count("T") + \
    seq.count("G") + seq.count("C")
  if (n_na >= int(n_valid * 0.9)):
    n_bases += len(seq)
  else :
    n_residues += len(seq)
  return n_residues, n_bases

def duplicate_multiple_chains(text):
  """make chains that are marked with >7WZE_1|Chains A, B[Auth X]| or similar
  N times that many chains"""
  new_groups = []
  next_lines = []
  unused_lines = []
  previous_line = None
  for line in text.splitlines():
    line = line.strip()
    if line.startswith(">"):
      next_lines = [line]
      new_groups.append(next_lines)
    elif not line.strip():
      next_lines = [">"]
      new_groups.append(next_lines)
    elif next_lines:
      next_lines.append(line)
    else:
      unused_lines.append(line)
    previous_line = line
  if unused_lines and new_groups:
    return text # could not do anything with this

  elif unused_lines:
    return text # nothing to do

  else: # usual
    new_lines = []
    for new_group in new_groups:
      if not new_group or not new_group[0]:
        continue
      first_line = new_group[0]
      n = get_number_of_dups(first_line)
      lines_in_group= new_group[1:]
      for i in range(n):
        new_lines.append(first_line)
        new_lines.append("".join(lines_in_group))
    text = "\n".join(new_lines)
    text = text.replace(">\n\n>",">")
    text = text.replace(">\n\n>",">")
    return text

def remove_inside_brackets(text):
  """ remove everything enclosed in []"""
  if not text.find("[")>-1:
    return text
  new_text = ""
  found_lb = False
  for c in text:
    if c=="[":
      found_lb = True
    elif c == "]":
      found_lb = False
    elif found_lb:
      pass # skip it
    else:
      new_text += c
  return new_text

def get_number_of_dups(line):
  """ Get number of duplicate chains: looks like >7WZE_1|Chains A, B[Auth X]| or similar """
  if not line.startswith(">"):
    return 1
  spl = line.split("|")
  if len(spl) < 2:
    return 1
  text = spl[1].strip()
  text = remove_inside_brackets(text)
  if not text.startswith("Chain"):
    return 1
  text = text.replace("Chains","").replace("Chain","").replace("=",""
       ).replace(",","")
  values = text.split()
  return max(1, len(values))

def clear_empty_lines(text, apply_duplicate_multiple_chains = False,
    keep_labels = False):
  """First duplicate any multiple chains, then clear empty lines."""
  if apply_duplicate_multiple_chains:
    text = duplicate_multiple_chains(text)
  # make empty lines just a blank line.  Includes >>> etc.
  # If keep_labels, make the starting line for each group start with >

  new_lines=[]
  prev_line = ""
  label_line = ""
  for line in text.splitlines():
    if keep_labels and line.startswith(">"):
      label_line = line
    if not line.replace(">","").replace(" ",""):
       line=""
    elif line.startswith(">"):
       line=""
    line=line.replace("?","")
    if (not line) and (not prev_line):
      continue # skip blanks if dup or at beginning
    if line and label_line:
      new_lines.append(label_line)
      label_line = ""
    new_lines.append(line)
    prev_line = line
  return "\n".join(new_lines)+"\n"

def get_sequences(file_name=None,text=None,remove_duplicates=None,
     apply_duplicate_multiple_chains = False,
     remove_unknowns = False,
     return_sequences_with_labels = False):
  """return simple list of sequences in this file. duplicates included
  unless remove_duplicates=True.
  remove unknowns (X) if requested.
  If return_sequences_with_labels,  return sequence objects with labels"""
  if not text:
    if not file_name:
      raise Sorry("Missing file for get_sequences: %s" %(
        file_name))
    with open(file_name) as f:
      text = f.read()
  # clear any lines that have only > and nothing else
  text=clear_empty_lines(text, apply_duplicate_multiple_chains,
    keep_labels = return_sequences_with_labels)

  ( sequences, unknowns ) = parse_sequence( text )

  simple_sequence_list=[]
  sequence_object_list = []
  for sequence in sequences:
    if remove_duplicates and sequence.sequence in simple_sequence_list:
      continue # it is a duplicate
    elif remove_unknowns: # remove any X and take it
      sequence.sequence = sequence.sequence.upper().replace("X","")
      simple_sequence_list.append(sequence.sequence)
      sequence_object_list.append(sequence)
    else: # take it
      sequence.sequence = sequence.sequence.upper()
      simple_sequence_list.append(sequence.sequence)
      sequence_object_list.append(sequence)
  if return_sequences_with_labels:
    return sequence_object_list
  else:
    return simple_sequence_list

#####################################################################
####   Methods to try and guess chain types from sequences ##########
#####################################################################
def guess_chain_types_from_sequences(file_name=None,text=None,
    return_as_dict=False,minimum_fraction=None,
    likely_chain_types=None):
  """Guess what chain types are in this sequence file"""
  if not text:
    if not file_name:
      from libtbx.utils import Sorry
      raise Sorry("Missing file for guess_chain_types_from_sequences: %s" %(
        file_name))
    with open(file_name) as f:
      text = f.read()
  # clear any lines that have only > and nothing else
  text=clear_empty_lines(text)

  chain_types=[]
  ( sequences, unknowns ) = parse_sequence( text )
  dd={}
  dd_n={}
  total_residues=0
  for sequence in sequences:
    chain_type,n_residues=chain_type_and_residues(text=sequence.sequence,
      likely_chain_types=likely_chain_types)
    if chain_type is None and n_residues is None:
      continue
    if chain_type and not chain_type in chain_types:
      chain_types.append(chain_type)
      dd[chain_type]=[]
      dd_n[chain_type]=0
    if chain_type:
      dd[chain_type].append(sequence)
      dd_n[chain_type]+=len(sequence.sequence)
      total_residues+=len(sequence.sequence)
  if minimum_fraction and len(chain_types)>1 and total_residues>1:
    # remove anything < minimum_fraction
    new_chain_types=[]
    new_dd={}
    new_dd_n={}
    for chain_type in chain_types:
      if dd_n[chain_type]/total_residues >= minimum_fraction:
        new_chain_types.append(chain_type)
        new_dd[chain_type]=dd[chain_type]
        new_dd_n[chain_type]=dd_n[chain_type]
    chain_types=new_chain_types
    dd=new_dd
    dd_n=new_dd_n
  if return_as_dict:
    return dd # dict of chain_types and sequences for each chain_type
  else:
    chain_types.sort()
    return chain_types

def text_from_chains_matching_chain_type(file_name=None,text=None,
    chain_type=None,width=80):
  """Get sequence text from all the chains that match chain_type"""
  dd=guess_chain_types_from_sequences(file_name=file_name,
    text=text,return_as_dict=True,likely_chain_types=[chain_type])
  sequence_text=""
  for ct in dd.keys():
    if chain_type is None or ct==chain_type:
      for seq in dd[ct]:
        sequence_text+="""
%s
 """ %(seq.format(width=width))
  return sequence_text

def count_letters(letters="",text="",only_count_non_allowed=None):
  """Count letter in text"""
  n=0
  if only_count_non_allowed: # count the ones that are not there
    for t in text:
      if not t in letters:
        n+=1

  else: # usual
    for let in letters:
      n+=text.count(let)
  return n

def chain_type_and_residues(text=None,chain_type=None,likely_chain_types=None):
  """guess the type of chain from text string containing 1-letter codes
  and count residues.
  If chain_type is specified, just use it.
  If likely_chain_types are specified, use them if possible.

   Assumptions:
    1. few or no letters that are not part of the correct dict (there
      may be a few like X or other unknowns)
    2. if all letters are are A is is because poly-ala and poly-gly
       are common for unknowns.
    3. otherwise if it can be DNA or protein it is DNA (because chance is
       very high that if it were protein there would be a non-DNA letter)
   Method:  Choose the chain-type that matches the most letters in text.
    If a tie, take the chain type that has the fewest letters.
    If likely_chain_types are specified, use one from there first

  """

  text=text.replace(" ","").replace("\n","").lower()
  if not text:
    return None,None
  letter_dict={
    'PROTEIN':"acdefghiklmnpqrstvwy",
    'DNA':"gact",
    'RNA':"gacu",}
  if chain_type not in [None,'None']:
    for key in list(letter_dict.keys()):
      if key != chain_type:
        del letter_dict[key]
  # Get all allowed letters
  all_letters=[]
  for chain_type in letter_dict.keys():
    for let in letter_dict[chain_type]:
      if not let in all_letters: all_letters.append(let)

  # remove non-allowed letters from text
  new_text=""
  for let in text:
    if let in all_letters:
      new_text+=let
  text=new_text
  if not text:
    return None,None
  # See which chain_type matches best
  count_dict={}
  non_allowed_count_dict={}
  for chain_type in letter_dict.keys():
    count_dict[chain_type]=count_letters(letters=letter_dict[chain_type],
      text=text)
    non_allowed_count_dict[chain_type]=count_letters(
      letters=letter_dict[chain_type],
      text=text,only_count_non_allowed=True)
  # Take max count_dict. If tie, take minimum non-allowed. If tie take the one
  #  with fewer letters (i.e., DNA instead of protein if matches both except
  #  poly-ala not poly-A)
  score_list=[]
  for chain_type in letter_dict.keys():
    score_list.append([count_dict[chain_type],chain_type])
  score_list.sort(key=itemgetter(0))
  score_list.reverse()
  ok_list=[]
  best_score=score_list[0][0]
  for score,chain_type in score_list:
    if score==best_score:
      ok_list.append(chain_type)
  residues=best_score
  if len(ok_list)<1: return None,None
  if len(ok_list)==1: return ok_list[0],residues

  # take one from the likely list
  if likely_chain_types:
    likely_results=[]
    for lct in likely_chain_types:
      if lct in ok_list:
        likely_results.append(lct)

    if len(likely_results)==1:
      return likely_results[0],residues

  # decide which of the ones with the most matches is best..
  score_list=[]
  for chain_type in ok_list:
    score_list.append([non_allowed_count_dict[chain_type],chain_type])
  score_list.sort(key=itemgetter(0))
  ok_list=[]
  best_score=score_list[0][0]
  for score,chain_type in score_list:
    if score==best_score:
      ok_list.append(chain_type)
  if len(ok_list)<1: return None,None
  if len(ok_list)==1: return ok_list[0],residues

  if text.replace("g","a")==residues*"a" and "PROTEIN" in letter_dict.keys():
    # special case, all Adenine or Ala
    return "PROTEIN",residues

  score_list=[]
  for chain_type in ok_list:
    score_list.append([len(letter_dict[chain_type]),chain_type])
  score_list.sort(key=itemgetter(0))
  ok_list=[]
  best_score=score_list[0][0]
  for score,chain_type in score_list:
    if score==best_score:
      ok_list.append(chain_type)
  if len(ok_list)<1:
    return None,None
  else:
    return ok_list[0],residues

#####################################################################
####   END OF methods to try and guess chain types from sequences ###
#####################################################################

def random_sequence(n_residues=None,residue_basket=None,
   chain_type = 'PROTEIN'):
  """Return n_residues random residues"""
  assert n_residues and (residue_basket or chain_type)
  if not residue_basket:
    chain_type = chain_type.upper()
    if chain_type == "PROTEIN":
        # Approximate eukaryotic frequencies using W as basic unit
        residue_basket = "AAAAAACCCEEEEDDDDDGGGGGG"+\
          "FFFIIIHHKKKKKKMLLLLLLNNNQQQPPPPSSSSSSRRRTTTTTWVVVVVYYY"
    elif chain_type == "DNA":
        residue_basket = "GATC"
    elif chain_type == "RNA":
        residue_basket = "GAUC"
    else:
        raise Sorry("Chain type needs to be RNA/DNA/PROTEIN")

  import random
  s=""
  nn=len(residue_basket)-1
  for i in range(n_residues):
    id=random.randint(0,nn)
    s+=residue_basket[id]
  return s
