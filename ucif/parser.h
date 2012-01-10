#ifndef UCIF_PARSER_H
#define UCIF_PARSER_H

#include <istream>
#include <string>

#include <ucif/cifLexer.h>
#include <ucif/cifParser.h>
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
      input = antlr3StringStreamNew(
        pANTLR3_UINT8(input_string.c_str()),
        ANTLR3_ENC_8BIT,
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
      psr->parse(psr, builder, strict);
      fflush(stderr);
    }

    ~parser()
    {
      // Essential to clean up after ourselves (in reverse order)
      delete psr->errors;
      delete lxr->errors;
      psr->free(psr);
      tstream->free(tstream);
      lxr->free(lxr);
      input->close(input);
    }

    pcifLexer lxr;
    pcifParser psr;

  private:
    pANTLR3_COMMON_TOKEN_STREAM tstream;
    pANTLR3_INPUT_STREAM input;

};

} // namespace ucif

#endif //GUARD
