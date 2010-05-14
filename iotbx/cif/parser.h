#ifndef IOTBX_CIF_PARSER_H
#define IOTBX_CIF_PARSER_H

#include <istream>
#include <string>

#include <iotbx/cif/cif2Lexer.h>
#include <iotbx/cif/cif2Parser.h>
#include <boost/python/object.hpp>

namespace iotbx { namespace cif {

class parser
{


  public:

    parser(std::string filename, boost::python::object& builder)
    {
      input = antlr3AsciiFileStreamNew(pANTLR3_UINT8(filename.c_str()));
      if (input == NULL)
      {
        input = antlr3NewAsciiStringInPlaceStream(pANTLR3_UINT8(filename.c_str()), filename.size(), pANTLR3_UINT8("memory"));
      }
      lxr = cif2LexerNew(input);
      tstream = antlr3CommonTokenStreamSourceNew(ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));
      psr = cif2ParserNew(tstream);
      psr->parse(psr, builder);
    }

    ~parser()
    {
      // Essential to clean up after ourselves (in reverse order)
      psr->free(psr);
      tstream->free(tstream);
      lxr->free(lxr);
      input->close(input);
    }

  private:
      pANTLR3_COMMON_TOKEN_STREAM       tstream;
      pANTLR3_INPUT_STREAM input;
      pcif2Lexer lxr;
      pcif2Parser psr;

};

}} // namespace iotbx::cif

#endif //GUARD
