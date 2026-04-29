"""
  Homology search by BLAST, XML output
"""
from __future__ import absolute_import, division, print_function

from iotbx.bioinformatics.xmlbuild import Single, Multiple
from iotbx.bioinformatics.xmlbuild import Text, TextWithDefault, Group
from iotbx.bioinformatics.xmlbuild import Parser

class Result(object):
  """
  Homology search by BLAST, XML output
  """

  def __init__(self, root):

    self.root = root


  def restrict(self, max_count):

    if self.root.iterations:
      del self.root.iterations[-1].hits[max_count:]


  def hits(self):

    if not self.root.iterations:
      raise StopIteration

    from iotbx import bioinformatics

    for h in self.root.iterations[-1].hits:
      pieces = h.accession.split( "_" )
      assert 0 < len( pieces )

      if len( pieces ) == 1:
        pdb = pieces[0]
        chain = ""

      else:
        ( pdb, chain ) = pieces[:2]

      alignment = bioinformatics.clustal_alignment(
        names = [ "target", "%s_%s" % ( pdb, chain ) ],
        alignments = [ h.hsps[0].query.seq, h.hsps[0].hit.seq ],
        program = "NCBI-BLAST"
        )
      yield bioinformatics.homology_search_hit(
        identifier = pdb,
        chain = chain,
        annotation = h.annotation,
        alignment = alignment
        )


  def __len__(self):

    if not self.root.iterations:
      return 0

    return len( self.root.iterations[-1].hits )


_alternatives = frozenset( "-" )

def substitute_underscore(data):

  return "".join( "_" if c in _alternatives else c for c in data )


BUILDER = Single(
  extractions = [
    Text( attribute = "program", tag = "BlastOutput_program" ),
    Text( attribute = "version", tag = "BlastOutput_version" ),
    Text( attribute = "reference", tag = "BlastOutput_reference" ),
    Text( attribute = "db", tag = "BlastOutput_db" ),
    Text( attribute = "query_id", tag = "BlastOutput_query-ID" ),
    Text( attribute = "query_def", tag = "BlastOutput_query-def" ),
    Text( attribute = "query_len", tag = "BlastOutput_query-len", conversion = int ),
    ],
  child_data_tagged = {
    "BlastOutput_param": (
      "parameters",
      Single(
        extractions = [
          Text( attribute = "matrix", tag = "Parameters/Parameters_matrix" ),
          Text( attribute = "expect", tag = "Parameters/Parameters_expect", conversion = float ),
          Text( attribute = "gap_open", tag = "Parameters/Parameters_gap-open", conversion = float ),
          Text( attribute = "gap_extend", tag = "Parameters/Parameters_gap-extend", conversion = float ),
          Text( attribute = "filter", tag = "Parameters/Parameters_filter" ),
          ],
        ),
      ),
    "BlastOutput_iterations": (
      "iterations",
      Multiple(
        tag = "Iteration",
        processor = Single(
          extractions = [
            Text( attribute = "num", tag = "Iteration_iter-num", conversion = int ),
            TextWithDefault( attribute = "query_id", tag = "Iteration_query-ID", default = "" ),
            TextWithDefault( attribute = "query_def", tag = "Iteration_query-def", default = "" ),
            TextWithDefault(
              attribute = "query_len",
              tag = "Iteration_query-len",
              conversion = int,
              default = 0,
              ),
            ],
          child_data_tagged = {
            "Iteration_hits": (
              "hits",
              Multiple(
                tag = "Hit",
                processor = Single(
                  extractions = [
                    Text( attribute = "num", tag = "Hit_num", conversion = int ),
                    Text( attribute = "identifier", tag = "Hit_id" ),
                    Text( attribute = "annotation", tag = "Hit_def" ),
                    Text( attribute = "accession", tag = "Hit_accession", conversion = substitute_underscore ),
                    Text( attribute = "length", tag = "Hit_len", conversion = int ),
                    ],
                  child_data_tagged = {
                    "Hit_hsps": (
                      "hsps",
                      Multiple(
                        tag = "Hsp",
                        processor = Single(
                          extractions = [
                            Text( attribute = "num", tag = "Hsp_num", conversion = int ),
                            Text( attribute = "bit_score", tag = "Hsp_bit-score", conversion = float ),
                            Text( attribute = "score", tag = "Hsp_score", conversion = float ),
                            Text( attribute = "evalue", tag = "Hsp_evalue", conversion = float ),
                            Group(
                              attribute = "query",
                              extractions = [
                                Text( attribute = "start", tag = "Hsp_query-from", conversion = int ),
                                Text( attribute = "end", tag = "Hsp_query-to", conversion = int ),
                                TextWithDefault(
                                  attribute = "frame",
                                  tag = "Hsp_query-frame",
                                  conversion = int,
                                  default = 0,
                                  ),
                                Text( attribute = "seq", tag = "Hsp_qseq" ),
                                ],
                              ),
                            Group(
                              attribute = "hit",
                              extractions = [
                                Text( attribute = "start", tag = "Hsp_hit-from", conversion = int ),
                                Text( attribute = "end", tag = "Hsp_hit-to", conversion = int ),
                                TextWithDefault(
                                  attribute = "frame",
                                  tag = "Hsp_hit-frame",
                                  conversion = int,
                                  default = 0,
                                  ),
                                Text( attribute = "seq", tag = "Hsp_hseq" ),
                                ],
                              ),
                            Text( attribute = "identity", tag = "Hsp_identity", conversion = int ),
                            Text( attribute = "positive", tag = "Hsp_positive", conversion = int ),
                            TextWithDefault( attribute = "gaps", tag = "Hsp_gaps", conversion = int, default = 0 ),
                            Text( attribute = "length", tag = "Hsp_align-len", conversion = int ),
                            Text( attribute = "midline", tag = "Hsp_midline" ),
                            ],
                          ),
                        ),
                      ),
                    },
                  ),
                ),
              ),
            "Iteration_stat": (
              "statistics",
              Single(
                extractions = [
                  Text( attribute = "db_num", tag = "Statistics/Statistics_db-num", conversion = int ),
                  Text( attribute = "db_len", tag = "Statistics/Statistics_db-len", conversion = int ),
                  Text( attribute = "hsp_len", tag = "Statistics/Statistics_hsp-len", conversion = int ),
                  Text(
                    attribute = "eff_space",
                    tag = "Statistics/Statistics_eff-space",
                    conversion = float,
                    ),
                  Text( attribute = "kappa", tag = "Statistics/Statistics_kappa", conversion = float ),
                  Text( attribute = "lambdav", tag = "Statistics/Statistics_lambda", conversion = float ),
                  Text( attribute = "entropy", tag = "Statistics/Statistics_entropy", conversion = float ),
                  ],
                ),
              ),
            },
          ),
        ),
      ),
    },
  )

parse = Parser( tag = "BlastOutput", builder = BUILDER, restype = Result, cElementTree = True )
