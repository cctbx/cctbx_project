import re
import operator
import os.path

# Wrap lines that are longer than 'width'
def wrap(text, width):

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
      raise ValueError, ( cls.format, str( data ) )

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
      raise ValueError, (
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


  def extend(self, sequences):

    if len( sequences ) != self.multiplicity():
      raise ValueError, "Inconsistent extension sequence set"

    pre_alis = []
    pre_gaps = []
    post_alis = []
    post_gaps = []

    for ( a, s ) in zip( self.alignments, sequences ):
      partial = "".join( [ c for c in a if c != self.gap ] )
      pos = s.find( partial )

      if pos == -1:
        raise ValueError, "Alignment does not match sequence"

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


  def _set_alignments(self, alignments):

    # All alignments should have the same length
    for o in alignments[1:]:
      if (len( alignments[0] ) != len( o )):
        raise ValueError, "'alignments' do not have the same length"

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


  def format(self, aln_width, caption_width, middle_line=None):

    if not middle_line:
      middle_line = midline().compare(
        alignments = self.alignments,
        gap = self.gap
        )

    elif len( middle_line ) != self.length():
      raise ValueError, "Incorrect midline length"

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
          self.type( **dict( kwargs.items() + match.groupdict().items() ) )
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


def partition_sequence_str(data, regex):
  """
  Partition sequence files
  """

  position = 0

  for match in regex.finditer( data ):
    current = match.start()
    unknown = data[ position : current ].strip()
    yield ( match.groups(), ( position, current, unknown ) )
    position = match.end()

  remaining = data[ position : ].strip()

  if remaining:
    yield ( (), ( position, len( data ), remaining ) )

  raise StopIteration


def parse_sequence_str(data, format):
  """
  Generic purpose sequence header parser
  """

  partition = partition_sequence_str( data = data, regex = format.regex )

  for ( groups, ( n_start, n_end, n_data ) ) in partition:
    if n_data:
      raise ValueError, (
        "Uninterpretable block from %s to %s:\nExpected: %s\nFound: %s" % (
          n_start,
          n_end,
          format.expected,
          n_data,
          )
        )

    body = "".join( [ c for c in groups[-1] if not c.isspace() ] )
    yield format.create( headers = groups[:-1], body = body )


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
    ^ > \s*
    (?P<name> [^\s^:]+ ) :? \s*
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
  Tom's format sequence parser
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

  ( name, extension ) = os.path.splitext( file_name )

  return _implemented_sequence_parsers.get( extension )


def sequence_parser_for_extension(extension):

  return _implemented_sequence_parsers.get( extension )


def known_sequence_formats():

  return _implemented_sequence_parsers.keys()

def any_sequence_format (file_name) :
  format_parser = sequence_parser_for(file_name)
  data = open(file_name, "r").read()
  seq_object = None
  if (format_parser is not None) :
    try :
      objects, non_compliant = format_parser.parse(data)
      assert (len(objects) > 0)
      assert (objects[0].sequence != "")
    except Exception, e :
      pass
    else :
      return objects, non_compliant
  for other_parser in [fasta_sequence_parse.parse, pir_sequence_parse.parse,
                       seq_sequence_parse.parse, tf_sequence_parse] :
    if (other_parser is not format_parser) :
      try :
        objects, non_compliant = other_parser(data)
        assert (len(objects) > 0)
        assert (objects[0].sequence != "")
      except Exception, e :
        pass
      else :
        return objects, non_compliant
  # fallback: unformatted
  data = re.sub("\s", "", data)
  if (re.search("[^a-zA-Z\*]", data) is None) :
    seq = sequence(data)
    return [seq], []
  return None, None

def merge_sequence_files (file_names, output_file,
    include_non_compliant=False) :
  assert (len(file_names) > 0)
  seq_out = open(output_file, "w")
  for seq_file in file_names :
    objects, non_compliant = any_sequence_format(seq_file)
    for seq_record in objects :
      name = str(seq_record.name)
      if (name == "") :
        name = "none"
      seq_out.write("> %s\n" % name)
      seq_out.write("%s\n" % seq_record.sequence)
  seq_out.close()

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


def alignment_parser_for_extension(extension):

  return _implemented_alignment_parsers.get( extension )


def known_alignment_formats():

  return _implemented_alignment_parsers.keys()


class homology_search_hit(object):
  """
  A hit from a homology search
  """

  def __init__(self, identifier, chain, annotation, alignment):

    self.identifier = identifier
    self.chain = chain
    self.annotation = annotation
    self.alignment = alignment


class hhpred_parser(object):
  """
  Parses .hhr files from HHPred
  """

  SPLIT = re.compile(
    r"""
    ( Query .*? ) \n\n
    ( \s+ No \s+ Hit .*? ) \n\n
    ( .*? )
    (?= (?:Done!) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL
    )

  HEADER = re.compile(
    r"""
    \A
    Query \s+ ( \S .* ) \n
    Match_columns \s+ ( \d+ ) \n
    No_of_seqs \s+ ( \d + ) \s out \s of \s ( \d+ ) \n
    Neff \s+ (\S[^\n]*) \n
    Searched_HMMs \s+ ( \d+ ) \n
    Date \s+ (\S[^\n]*) \n
    Command \s+ (\S.*)
    \Z
    """,
    re.VERBOSE
    )


  def __init__(self, output):

    res = self.SPLIT.search( output )

    if not res:
      raise ValueError, "Incorrect file format"

    ( header, summary, hits ) = res.groups()
    self.process_header( header = header )
    self.process_hits( hits = hits )


  def process_header(self, header):

    res = self.HEADER.search( header )

    if not res:
      raise ValueError, "Incorrect header format"

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
          raise ValueError, "Uninterpretable block: %s" % repr( unknown )

      start = m.end()

      self.add_match_to_hit_header_results( match = m )
      self.process_blocks( blocks = m.groupdict()[ "blocks" ] )

    remaining = hits[ start: ].strip()

    if remaining:
      raise ValueError, "Uninterpretable tail: %s" % repr( remaining )


  def process_blocks(self, blocks):

    matches = []
    start = 0

    for match in self.BLOCKS.finditer( blocks ):
      if match.start() != start:
        assert start < match.start()

        unknown = blocks[ start : match.start() ].strip()

        if unknown:
          raise ValueError, "Uninterpretable: %s" % repr( unknown )

      start = match.end()
      matches.append( match.groups() )

    remaining = blocks[ start: ].strip()

    if remaining:
      raise ValueError, "Uninterpretable: %s" % repr( remaining )

    if not matches:
      raise ValueError, "Empty homology block"

    assert all([len( matches[0] ) == len( a ) for a in matches[1:]])

    return self.merge_and_process_block_hits( matches = matches )


  def merge_sequence_numbers(self, starts, ends, others):

    for ( s, e ) in zip( starts[1:], ends[:-1] ):
      if s != e + 1:
        raise ValueError, "Incorrect sequence indices"

      if not all([others[0] == o for o in others[1:]]):
        raise ValueError, "Incorrect sequence indices"

    return ( starts[0], ends[-1], others[0] )


class hhsearch_parser(hhpred_parser):
  """
  Specific for hhsearch output
  """

  HITS = re.compile(
    r"""
    No \s ( \d+) \s* \n
    >( [\w]{4} ) _  ( [\w]? ) \s ( [^\n]* )\n
    Probab = ( [+-]? \d+ \. \d* ) \s+
    E-value = ( \d+ \.? \d* )( e[+-]? \d+ )? \s+
    Score = ( \d+\.\d+ ) \s+
    Aligned_cols = ( \d+ ) \s+
    Identities = ( \d+ ) % \s+
    Similarity = ( -? \d+ \. \d+ ) \s+
    Sum_probs = ( \d+ \. \d+ ) \n
    (?P<blocks> .*? )(?= (?:^No) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL | re.MULTILINE
    )
  BLOCKS = re.compile(
    r"""
    Q \s+ ss_pred \s+ ( [\w-]+ ) \s* \n
    Q \s+ [\w:\.]+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    Q \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    \s+ ( [ \.\-+|=]+ ) \n
    T \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    T \s+ \w+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    T \s+ ss_dssp \s+ ( [\w-]+ ) \s* \n
    T \s+ ss_pred \s+ ( [\w-]+ ) \s* \n
    (?: Confidence \s+ ( [\w ]+ ) \s* \n)?
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
    self.similarities = []
    self.sum_probs = []

    self.query_starts = []
    self.query_ends = []
    self.query_others = []
    self.query_alignments = []
    self.query_consensi = []
    self.query_ss_preds = []

    self.midlines = []

    self.hit_starts = []
    self.hit_ends = []
    self.hit_others = []
    self.hit_alignments = []
    self.hit_consensi = []
    self.hit_ss_preds = []
    self.hit_ss_dssps = []


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
    self.similarities = self.similarities[:max_count]
    self.sum_probs = self.sum_probs[:max_count]

    self.query_starts = self.query_starts[:max_count]
    self.query_ends = self.query_ends[:max_count]
    self.query_others = self.query_others[:max_count]
    self.query_alignments = self.query_alignments[:max_count]
    self.query_consensi = self.query_consensi[:max_count]
    self.query_ss_preds = self.query_ss_preds[:max_count]

    self.midlines = self.midlines[:max_count]

    self.hit_starts = self.hit_starts[:max_count]
    self.hit_ends = self.hit_ends[:max_count]
    self.hit_others = self.hit_others[:max_count]
    self.hit_alignments = self.hit_alignments[:max_count]
    self.hit_consensi = self.hit_consensi[:max_count]
    self.hit_ss_preds = self.hit_ss_preds[:max_count]
    self.hit_ss_dssps = self.hit_ss_dssps[:max_count]


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
    self.similarities.append( float( match.group( 11 ) ) )
    self.sum_probs.append( float( match.group( 12 ) ) )


  def merge_and_process_block_hits(self, matches):

    data = zip( *matches )
    assert len( data ) == 21

    sequences = [ data[0], data[2], data[6], data[11], data[15], data[18], data[19] ]
    midlines = []

    for ( index , alis) in enumerate( zip( *sequences ) ):
      count = len( alis[0] )

      if not all([count == len( c ) for c in alis[1:]]):
        raise ValueError, "Incorrect alignments"

      midlines.append(
        " " * ( count - len( data[9][ index ] ) ) + data[9][ index ]
        )

    merged = [ reduce( operator.add, a ) for a in sequences ]
    assert len( midlines ) == len( matches )
    midline = reduce( operator.add, midlines )

    if data[1] != data[5] or data[3] != data[7] or data[4] != data[8]:
      raise ValueError, "Inconsistent query numbering"

    q_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[1] ],
      ends = [ int( d ) for d in data[3] ],
      others = [ int( d ) for d in data[4] ]
      )

    if data[10] != data[14] or data[12] != data[16] or data[13] != data[17]:
      raise ValueError, "Inconsistent target numbering"

    t_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[10] ],
      ends = [ int( d ) for d in data[12] ],
      others = [ int( d ) for d in data[13] ]
      )

    self.query_starts.append( q_indices[0] )
    self.query_ends.append( q_indices[1] )
    self.query_others.append( q_indices[2] )
    self.query_ss_preds.append( merged[0] )
    self.query_alignments.append( merged[1] )
    self.query_consensi.append( merged[2] )

    self.midlines.append( midline )

    self.hit_starts.append( t_indices[0] )
    self.hit_ends.append( t_indices[1] )
    self.hit_others.append( t_indices[2] )
    self.hit_consensi.append( merged[3] )
    self.hit_alignments.append( merged[4] )
    self.hit_ss_dssps.append( merged[5] )
    self.hit_ss_preds.append( merged[6] )


  def hits(self):

    data = zip(
      self.pdbs,
      self.chains,
      self.annotations,
      self.query_alignments,
      self.hit_alignments,
      )

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


class hhalign_parser(hhpred_parser):
  """
  Specific for hhalign output
  """

  HITS = re.compile(
    r"""
    No \s ( \d+) \s* \n
    > \s* ( [^\n]* ) \n
    Probab = ( [+-]? \d+ \. \d* ) \s+
    E-value = ( \d+ \.? \d* )( e[+-]? \d+ )? \s+
    Score = ( \d+\.\d+ ) \s+
    Aligned_columns = ( \d+ ) \s+
    Identities = ( \d+ ) % \n
    (?P<blocks> .*? )(?= (?:^No) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL | re.MULTILINE
    )
  BLOCKS = re.compile(
    r"""
    Q \s+ ss_pred \s+ ( [\w-]+ ) \s* \n
    Q \s+ ss_conf \s+ ( [\w-]+ ) \s* \n
    Q \s+ [\w:\.]+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    Q \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    \s+ ( [ \.\-+|=]+ ) \n
    T \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    T \s+ [\w\.]+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    T \s+ ss_pred \s+ ( [\w-]+ ) \s* \n
    T \s+ ss_conf \s+ ( [\w-]+ ) \s* \n
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

    data = zip( *matches )
    assert len( data ) == 21

    sequences = [
        data[0], data[1], data[3], data[7],
        data[12], data[16], data[19], data[20] ]
    midlines = []

    for ( index , alis) in enumerate( zip( *sequences ) ):
      count = len( alis[0] )

      if not all([count == len( c ) for c in alis[1:]]):
        raise ValueError, "Incorrect alignments"

      midlines.append(
        " " * ( count - len( data[10][ index ] ) ) + data[10][ index ]
        )

    merged = [ reduce( operator.add, a ) for a in sequences ]
    assert len( midlines ) == len( matches )
    midline = reduce( operator.add, midlines )

    if data[2] != data[6] or data[4] != data[8] or data[5] != data[9]:
      raise ValueError, "Inconsistent query numbering"

    q_indices = self.merge_sequence_numbers(
      starts = [ int( d ) for d in data[2] ],
      ends = [ int( d ) for d in data[4] ],
      others = [ int( d ) for d in data[5] ]
      )

    if data[11] != data[15] or data[13] != data[17] or data[14] != data[18]:
      raise ValueError, "Inconsistent target numbering"

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
