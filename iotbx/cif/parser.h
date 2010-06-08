#ifndef IOTBX_CIF_PARSER_H
#define IOTBX_CIF_PARSER_H

#include <istream>
#include <string>

#include <iotbx/cif/cifLexer.h>
#include <iotbx/cif/cifParser.h>
#include <boost/python/object.hpp>
#include <boost/noncopyable.hpp>

namespace iotbx { namespace cif {

class parser : private boost::noncopyable
{

  public:

    parser() {}

    parser(std::string filename, boost::python::object& builder)
    {
      input = antlr3AsciiFileStreamNew(pANTLR3_UINT8(filename.c_str()));
      if (input == NULL)
      {
        input = antlr3NewAsciiStringInPlaceStream(pANTLR3_UINT8(
          filename.c_str()), filename.size(), pANTLR3_UINT8("memory"));
      }
      lxr = cifLexerNew(input);
      tstream = antlr3CommonTokenStreamSourceNew(
        ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));
      psr = cifParserNew(tstream);
      psr->parse(psr, builder);
      fflush(stderr);
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
      pcifLexer lxr;
      pcifParser psr;

};

}} // namespace iotbx::cif

#endif //GUARD
