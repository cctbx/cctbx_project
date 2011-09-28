#ifndef UCIF_PARSER_H
#define UCIF_PARSER_H

#include <istream>
#include <string>

#include <ucif/cifLexer.h>
#include <ucif/cifParser.h>
#include <ucif/cifWalker.h>
#include <ucif/utils.h>

namespace ucif {

class parser
{
  private:  // emphasize the following members are private
    parser(const parser&);
    const parser& operator=(const parser&);

  public:

    parser() {}

    parser(builder_base* builder,
           std::string input_string,
           std::string filename="memory", bool strict=true)
    {
      input = antlr3NewAsciiStringInPlaceStream(
        pANTLR3_UINT8(input_string.c_str()),
        input_string.size(),
        pANTLR3_UINT8(filename.c_str()));
      lxr = cifLexerNew(input);
      tstream = antlr3CommonTokenStreamSourceNew(
        ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));
      psr = cifParserNew(tstream);
      psr->pParser->rec->displayRecognitionError = parser_displayRecognitionError;
      psr->errors = builder->new_array();
      lxr->pLexer->rec->displayRecognitionError = lexer_displayRecognitionError;
      lxr->errors = builder->new_array();
      cif_AST = psr->parse(psr, strict);
      if (lxr->errors->size() == 0 && psr->errors->size() == 0 && cif_AST.tree != NULL) {
        nodes   = antlr3CommonTreeNodeStreamNewTree(cif_AST.tree, ANTLR3_SIZE_HINT);
        tree_psr        = cifWalkerNew(nodes);
        tree_psr->errors = builder->new_array();
        tree_psr->parse(tree_psr, builder);
      }
      else {
        tree_psr = NULL;
      }
      fflush(stderr);
    }

    ~parser()
    {
      // Essential to clean up after ourselves (in reverse order)
      if (tree_psr != NULL) {
        nodes->free(nodes);
        nodes = NULL;
        delete tree_psr->errors;
        tree_psr->free(tree_psr);
        tree_psr = NULL;
      }
      delete psr->errors;
      delete lxr->errors;
      psr->free(psr);
      tstream->free(tstream);
      lxr->free(lxr);
      input->close(input);
    }

    pcifLexer lxr;
    pcifParser psr;
    pcifWalker tree_psr;

  private:
    pANTLR3_COMMON_TOKEN_STREAM tstream;
    pANTLR3_INPUT_STREAM input;
    pANTLR3_COMMON_TREE_NODE_STREAM     nodes;
    cifParser_parse_return cif_AST;

};

} // namespace ucif

#endif //GUARD
