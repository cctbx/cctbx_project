#ifndef IOTBX_CIF_PARSER_H
#define IOTBX_CIF_PARSER_H

#include <istream>
#include <string>

#include <iotbx/cif/cifLexer.h>
#include <iotbx/cif/cifParser.h>
#include <iotbx/cif/cifWalker.h>
#include <iotbx/cif/utils.h>

#include <boost/noncopyable.hpp>

namespace iotbx { namespace cif {

class parser : private boost::noncopyable
{

  public:

    parser() {}

    parser(std::string input_string, builder_base* builder, bool strict=true)
    {
      input = antlr3NewAsciiStringInPlaceStream(pANTLR3_UINT8(
        input_string.c_str()), input_string.size(), pANTLR3_UINT8("memory"));
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

}} // namespace iotbx::cif

#endif //GUARD
