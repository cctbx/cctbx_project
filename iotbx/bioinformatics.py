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

    self._set_alignments( alignments = alignments )

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
  ".fa": fasta_sequence_parse,
  ".faa" : fasta_sequence_parse,
  ".pir": pir_sequence_parse,
  ".seq": seq_sequence_parse,
  ".dat": seq_sequence_parse,
  }

def sequence_parser_for(file_name):

  ( name, extension ) = os.path.splitext( file_name )

  return _implemented_sequence_parsers.get( extension )


def known_sequence_formats():

  return _implemented_sequence_parsers.keys()

def any_sequence_format (file_name) :
  format_parser = sequence_parser_for(file_name)
  data = open(file_name, "r").read()
  seq_object = None
  if (format_parser is not None) :
    try :
      seq_object = format_parser.parse(data)
      assert seq_object is not None
    except Exception, e :
      pass
    else :
      return seq_object
  for other_parser in [fasta_sequence_parse, pir_sequence_parse,
                       seq_sequence_parse] :
    if (other_parser is not format_parser) :
      try :
        seq_object = other_parser.parse(data)
        assert seq_object is not None
      except Exception, e :
        pass
      else :
        return seq_object
  return None

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


class hhpred_homology_search(object):
  """
  Parses .hhr files from HHPred
  """

  SPLIT = re.compile(
    r"""
    ( Query .*? ) \n\n
    ( \s+ No \s+ Hit .*? ) \n\n
    ( .* )
    Done!
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

  HIT = re.compile(
    r"""
    No \s ( \d+) \s* \n
    >( [\w]{4} ) _  ( [\w]? ) \s ([^\n]*)\n
    Probab = ( [+-]? \d+ \. \d* ) \s+
    E-value = ( \d+ \. \d+ )( e[+-]? \d+ )? \s+
    Score = ( \d+\.\d+ ) \s+
    Aligned_cols = ( \d+ ) \s+
    Identities = ( \d+ ) % \s+
    Similarity = ( \d+ \. \d+ ) \s+
    Sum_probs = ( \d+ \. \d+ ) \n
    ( .*? )(?= (?:^No) | (?:\Z) )
    """,
    re.VERBOSE | re.DOTALL | re.MULTILINE
    )

  BLOCK = re.compile(
    r"""
    Q \s+ ss_pred \s+ ( [\w-]+ ) \s* \n
    Q \s+ [\w:]+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    Q \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    \s+ ( [ \.\-+|]+ ) \n
    T \s+ Consensus \s+ ( \d+ ) \s+ ( [\w~-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    T \s+ \w+ \s+ ( \d+ ) \s+ ( [\w-]+ ) \s+ ( \d+ ) \s+ \( ( \d+ ) \) \s* \n
    T \s+ ss_dssp \s+ ( [\w-]+ ) \s* \n
    T \s+ ss_pred \s+ ( [\w-]+ ) \s* \n
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


    for m in self.HIT.finditer( hits ):
      self.indices.append( int( m.group( 1 ) ) )
      self.pdbs.append( m.group( 2 ) )
      self.chains.append( m.group( 3 ) )
      self.annotations.append( m.group( 4 ) )
      self.probabs.append( float( m.group( 5 ) ) )
      self.e_values.append( float( m.group( 6 ) + m.group( 7 ) ) )
      self.scores.append( float( m.group( 8 ) ) )
      self.aligned_cols.append( int( m.group( 9 ) ) )
      self.identities.append( float( m.group( 10 ) ) )
      self.similarities.append( float( m.group( 11 ) ) )
      self.sum_probs.append( float( m.group( 12 ) ) )

      ( q_indices, h_indices, alignments, midline ) = self.process_blocks(
        blocks = m.group( 13 )
        )
      self.query_starts.append( q_indices[0] )
      self.query_ends.append( q_indices[1] )
      self.query_others.append( q_indices[2] )
      self.query_ss_preds.append( alignments[0] )
      self.query_alignments.append( alignments[1] )
      self.query_consensi.append( alignments[2] )

      self.midlines.append( midline )

      self.hit_starts.append( h_indices[0] )
      self.hit_ends.append( h_indices[1] )
      self.hit_others.append( h_indices[2] )
      self.hit_consensi.append( alignments[3] )
      self.hit_alignments.append( alignments[4] )
      self.hit_ss_dssps.append( alignments[5] )
      self.hit_ss_preds.append( alignments[6] )


  @classmethod
  def process_blocks(cls, blocks):

    data = [ [] for i in range( 20 ) ]

    for b in cls.BLOCK.finditer( blocks ):
      data[0].append( b.group( 1 ) ) # q_ss_pred

      data[1].append( int( b.group( 2 ) ) ) # q_sequence_start_indices
      data[2].append( b.group( 3 ) ) # q_sequence
      data[3].append( int( b.group( 4 ) ) ) # q_sequence_end_indices
      data[4].append( int( b.group( 5 ) ) ) # q_sequence_end_other_indices

      data[5].append( int( b.group( 6 ) ) ) # q_consensus_start_indides
      data[6].append( b.group( 7 ) ) # q_consensus
      data[7].append( int( b.group( 8 ) ) ) # q_consensus_end_indices
      data[8].append( int( b.group( 9 ) ) ) # q_consensus_end_other_indices

      data[9].append( b.group( 10 ) ) # midline

      data[10].append( int( b.group( 11 ) ) ) # t_consensus_start_indides
      data[11].append( b.group( 12 ) ) # t_consensus
      data[12].append( int( b.group( 13 ) ) ) # t_consensus_end_indices
      data[13].append( int( b.group( 14 ) ) ) # t_consensus_end_other_indices

      data[14].append( int( b.group( 15 ) ) ) # t_sequence_start_indices
      data[15].append( b.group( 16 ) ) # t_sequence
      data[16].append( int( b.group( 17 ) ) ) # t_sequence_end_indices
      data[17].append( int( b.group( 18 ) ) ) # t_sequence_end_other_indices

      data[18].append( b.group( 19 ) ) # t_ss_dssp
      data[19].append( b.group( 20 ) ) # t_ss_pred


    if not data[0]:
      raise ValueError, "Incorrect homology block format"

    if not all( len( data[0] ) == len( a ) for a in data[1:] ):
      raise ValueError, "Inconsistent alignments"

    sequences = [ data[0], data[2], data[6], data[11], data[15], data[18], data[19] ]

    for ( index , alis) in enumerate( zip( *sequences ) ):
      count = len( alis[0] )

      if not all( count == len( c ) for c in alis[1:] ):
        raise ValueError, "Incorrect alignments"

      data[9][index] = " " * ( count - len( data[9][ index ] ) ) + data[9][ index ]

    merged = [ reduce( operator.add, a ) for a in sequences ]

    if data[1] != data[5] or data[3] != data[7] or data[4] != data[8]:
      raise ValueError, "Inconsistent query numbering"

    q_indices = cls.merge_sequence_numbers(
      starts = data[1],
      ends = data[3],
      others = data[4]
      )

    if data[10] != data[14] or data[12] != data[16] or data[13] != data[17]:
      raise ValueError, "Inconsistent target numbering"

    t_indices = cls.merge_sequence_numbers(
      starts = data[10],
      ends = data[12],
      others = data[13]
      )

    return ( q_indices, t_indices, merged, reduce( operator.add, data[9] ) )


  @classmethod
  def merge_sequence_numbers(cls, starts, ends, others):

    for ( s, e ) in zip( starts[1:], ends[:-1] ):
      if s != e + 1:
        raise ValueError, "Incorrect sequence indices"

      if not all( others[0] == o for o in others[1:] ):
        raise ValueError, "Incorrect sequence indices"

    return ( starts[0], ends[-1], others[0] )
