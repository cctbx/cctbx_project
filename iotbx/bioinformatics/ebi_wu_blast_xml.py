"""Carry out homology searach by WU-BLAST at EBI"""

from __future__ import absolute_import, division, print_function

from iotbx.bioinformatics.xmlbuild import Single, Multiple, DataAttribute
from iotbx.bioinformatics.xmlbuild import Text, TextWithDefault, Attribute, AttributeWithDefault, Value
from iotbx.bioinformatics.xmlbuild import Parser

class Result(object):
  """
  Homology search by WU-BLAST at EBI
  """

  def __init__(self, root):

    self.root = root


  def restrict(self, max_count):

    del self.root.hits[max_count:]


  def hits(self):

    from iotbx import bioinformatics

    for h in self.root.hits:
      pieces = h.identifier.split( "_" )
      assert 0 < len( pieces )

      if len( pieces ) == 1:
        pdb = pieces[0]
        chain = ""

      else:
        ( pdb, chain ) = pieces[:2]

      alignment = bioinformatics.clustal_alignment(
        names = [ "target", "%s_%s" % ( pdb, chain ) ],
        alignments = [ h.alignments[0].query.seq, h.alignments[0].match.seq ],
        program = "WU-BLAST"
        )
      yield bioinformatics.homology_search_hit(
        identifier = pdb,
        chain = chain,
        annotation = h.description,
        alignment = alignment
        )

  def __len__(self):

    return len( self.root.hits )


from xml.etree.cElementTree import QName

def ebi_qualified(tag):

  return QName( "http://www.ebi.ac.uk/schema", tag ).text


BUILDER = Single(
  child_data_tagged = {
    ebi_qualified( tag = "Header" ): (
      "header",
      Single(
        child_data_tagged = {
          ebi_qualified( tag = "program" ): (
            "program",
            Single(
              extractions = [
                Attribute( attribute = "name", name = "name" ),
                Attribute( attribute = "version", name = "version" ),
                AttributeWithDefault(
                  attribute = "citation",
                  name = "citation",
                  default = "",
                  ),
                ],
              ),
            ),
          ebi_qualified( tag = "commandLine" ): ( "command", DataAttribute( name = "command" ) ),
          ebi_qualified( tag = "parameters" ): (
            "parameters",
            Single(
              extractions = [
                TextWithDefault(
                  attribute = "scores",
                  tag = ebi_qualified( "scores" ),
                  conversion = float,
                  default = 0.0,
                  ),
                TextWithDefault(
                  attribute = "alignments",
                  tag = ebi_qualified( "alignments" ),
                  conversion = int,
                  default = 0,
                  ),
                TextWithDefault(
                  attribute = "matrix",
                  tag = ebi_qualified( "matrix" ),
                  default = "",
                  ),
                TextWithDefault(
                  attribute = "expectation_upper",
                  tag = ebi_qualified( "expectationUpper" ),
                  conversion = float,
                  default = 0.0,
                  ),
                TextWithDefault(
                  attribute = "statistics",
                  tag = ebi_qualified( "statistics" ),
                  default = "",
                  ),
                TextWithDefault(
                  attribute = "filter",
                  tag = ebi_qualified( "filter" ),
                  default = "",
                  ),
                ],
              child_data_tagged = {
                ebi_qualified( tag = "sequences" ): (
                  "sequences",
                  Multiple(
                    tag = ebi_qualified( tag = "sequence" ),
                    processor = Single(
                      extractions = [
                        Attribute( attribute = "number", name = "number", conversion = int ),
                        Attribute( attribute = "name", name = "name" ),
                        Attribute( attribute = "type", name = "type" ),
                        Attribute( attribute = "length", name = "length", conversion = int ),
                        ]
                      ),
                    ),
                  ),
                ebi_qualified( tag = "databases" ): (
                  "databases",
                  Multiple(
                    tag = ebi_qualified( tag = "database" ),
                    processor = Single(
                      extractions = [
                        Attribute( attribute = "number", name = "number", conversion = int ),
                        Attribute( attribute = "name", name = "name" ),
                        Attribute( attribute = "type", name = "type" ),
                        AttributeWithDefault( attribute = "created", name = "created", default = "" ),
                        ]
                      ),
                    ),
                  ),
                },
              ),
            ),
          ebi_qualified( tag = "timeInfo" ): (
            "time",
            Single(
              extractions = [
                Attribute( attribute = "start", name = "start" ),
                Attribute( attribute = "end", name = "end" ),
                Attribute( attribute = "search", name = "search" ),
                ],
              ),
            ),
          },
        )
      ),
    ebi_qualified( "hits" ): (
      "hits",
      Multiple(
        tag = ebi_qualified( tag = "hit" ),
        processor = Single(
          extractions = [
            Attribute( attribute = "number", name = "number", conversion = int ),
            Attribute( attribute = "database", name = "database" ),
            Attribute( attribute = "identifier", name = "id" ),
            Attribute( attribute = "length", name = "length", conversion = int ),
            Attribute( attribute = "description", name = "description" ),
            ],
          child_data_tagged = {
            ebi_qualified( tag = "alignments" ): (
              "alignments",
              Multiple(
                tag = ebi_qualified( tag = "alignment" ),
                processor = Single(
                  extractions = [
                    Attribute( attribute = "number", name = "number", conversion = int ),
                    Text(
                      attribute = "score",
                      tag = ebi_qualified( tag = "score" ),
                      conversion = float,
                      ),
                    Text(
                      attribute = "bits",
                      tag = ebi_qualified( tag = "bits" ),
                      conversion = float,
                      ),
                    Text(
                      attribute = "expectation",
                      tag = ebi_qualified( tag = "expectation" ),
                      conversion = float,
                      ),
                    TextWithDefault(
                      attribute = "probability",
                      tag = ebi_qualified( tag = "probability" ),
                      conversion = float,
                      default = 0.0,
                      ),
                    Text(
                      attribute = "identity",
                      tag = ebi_qualified( tag = "identity" ),
                      conversion = float,
                      ),
                    Text(
                      attribute = "positives",
                      tag = ebi_qualified( tag = "positives" ),
                      conversion = float,
                      ),
                    Text( attribute = "pattern", tag = ebi_qualified( tag = "pattern" ) ),
                    ],
                  child_data_tagged = {
                    ebi_qualified( tag = "querySeq" ): (
                      "query",
                      Single(
                        extractions = [
                          Attribute( attribute = "start", name = "start", conversion = int ),
                          Attribute( attribute = "end", name = "end", conversion = int ),
                          Value( attribute = "seq" ),
                          ]
                        )
                      ),
                    ebi_qualified( tag = "matchSeq" ): (
                      "match",
                      Single(
                        extractions = [
                          Attribute( attribute = "start", name = "start", conversion = int ),
                          Attribute( attribute = "end", name = "end", conversion = int ),
                          Value( attribute = "seq" ),
                          ]
                        )
                      ),
                    },
                  ),
                ),
              ),
            }
          ),
        ),
      ),
    },
  )

parse = Parser(
  tag = ebi_qualified( tag = "EBIApplicationResult" ),
  builder = BUILDER,
  restype = Result,
  cElementTree = True,
  )

