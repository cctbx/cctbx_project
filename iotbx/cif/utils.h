#ifndef IOTBX_CIF_UTILS_H
#define IOTBX_CIF_UTILS_H

#include <string>

#include <boost/format.hpp>

#include <iotbx/cif/cifParser.h>

// Workaround conflict with min/max macros defined by Windef.h
// http://support.microsoft.com/kb/143208
#if defined(min) && defined(max)
  #define min_redefined min
  #define max_redefined max
  #undef min
  #undef max
#endif
#include <scitbx/array_family/shared.h>
#ifdef min_redefined
  #define max max_redefined
  #define min min_redefined
#endif

namespace iotbx { namespace cif {

/// Error display method taken from antlr3baserecognizer.c
/// Modified to print output to C++ iostream (to enable Python redirection).

/// Standard/Example error display method.
/// No generic error message display funciton coudl possibly do everything correctly
/// for all possible parsers. Hence you are provided with this example routine, which
/// you should override in your parser/tree parser to do as you will.
///
/// Here we depart somewhat from the Java runtime as that has now split up a lot
/// of the error display routines into spearate units. However, ther is little advantage
/// to this in the C version as you will probably implement all such routines as a
/// separate translation unit, rather than install them all as pointers to functions
/// in the base recognizer.
///
static void
parser_displayRecognitionError (pANTLR3_BASE_RECOGNIZER recognizer, pANTLR3_UINT8 * tokenNames)

{
  pANTLR3_PARSER              parser;
  pANTLR3_TREE_PARSER         tparser;
  pANTLR3_INT_STREAM          is;
  pANTLR3_STRING              ttext;
  pANTLR3_STRING              ftext;
  pANTLR3_EXCEPTION           ex;
  pANTLR3_COMMON_TOKEN        theToken;
  pANTLR3_BASE_TREE           theBaseTree;
  pANTLR3_COMMON_TREE         theCommonTree;

  std::string message;
    // Retrieve some info for easy reading.
    //
    ex      = recognizer->state->exception;
    ttext   = NULL;

    // See if there is a 'filename' we can use
    //
    if (ex->streamName == NULL)
    {
      if (((pANTLR3_COMMON_TOKEN)(ex->token))->type == ANTLR3_TOKEN_EOF)
      {
        message += "-end of input-(";
      }
      else
      {
        message += "-unknown source-(";
      }
    }
    else
    {
      ftext = ex->streamName->to8(ex->streamName);
      message += str(boost::format("%s(") %ftext->chars);
    }

    // Next comes the line number
    //

    message += str(boost::format("line %d) ") %recognizer->state->exception->line);
    message += str(boost::format(" : error %d : %s")
      %recognizer->state->exception->type
      %(pANTLR3_UINT8) (recognizer->state->exception->message));


    // How we determine the next piece is dependent on which thing raised the
    // error.
    //
    switch (recognizer->type)
    {
    case ANTLR3_TYPE_PARSER:

      // Prepare the knowledge we know we have
      //
      parser      = (pANTLR3_PARSER) (recognizer->super);
      tparser     = NULL;
      is          = parser->tstream->istream;
      theToken    = (pANTLR3_COMMON_TOKEN)(recognizer->state->exception->token);
      //ttext      = theToken->toString(theToken);
      ttext       =  theToken->getText(theToken);

      message += str(boost::format(", at offset %d") %recognizer->state->exception->charPositionInLine);
      if  (theToken != NULL)
      {
        if (theToken->type == ANTLR3_TOKEN_EOF)
        {
          message += ", at <EOF>";
        }
        else
        {
          // Guard against null text in a token
          //
          message += str(boost::format("\n    near %s\n    ") %(ttext == NULL ? (pANTLR3_UINT8)"<no text for the token>" : ttext->chars));
        }
      }
      break;

    case ANTLR3_TYPE_TREE_PARSER:

      tparser         = (pANTLR3_TREE_PARSER) (recognizer->super);
      parser          = NULL;
      is                      = tparser->ctnstream->tnstream->istream;
      theBaseTree     = (pANTLR3_BASE_TREE)(recognizer->state->exception->token);
      ttext           = theBaseTree->toStringTree(theBaseTree);

      if (theBaseTree != NULL)
      {
        theCommonTree = (pANTLR3_COMMON_TREE) theBaseTree->super;

        if (theCommonTree != NULL)
        {
          theToken = (pANTLR3_COMMON_TOKEN) theBaseTree->getToken(theBaseTree);
        }
        message += str(boost::format(", at offset %d") %theBaseTree->getCharPositionInLine(theBaseTree));
        message += str(boost::format(", near %s") %ttext->chars);
      }
      break;

    default:

      message += "Base recognizer function displayRecognitionError called by unknown parser type - provide override for this function\n";
      return;
      break;
    }

    // Although this function should generally be provided by the implementation, this one
    // should be as helpful as possible for grammar developers and serve as an example
    // of what you can do with each exception type. In general, when you make up your
    // 'real' handler, you should debug the routine with all possible errors you expect
    // which will then let you be as specific as possible about all circumstances.
    //
    // Note that in the general case, errors thrown by tree parsers indicate a problem
    // with the output of the parser or with the tree grammar itself. The job of the parser
    // is to produce a perfect (in traversal terms) syntactically correct tree, so errors
    // at that stage should really be semantic errors that your own code determines and handles
    // in whatever way is appropriate.
    //
    switch  (ex->type)
    {
    case ANTLR3_UNWANTED_TOKEN_EXCEPTION:

      // Indicates that the recognizer was fed a token which seesm to be
      // spurious input. We can detect this when the token that follows
      // this unwanted token would normally be part of the syntactically
      // correct stream. Then we can see that the token we are looking at
      // is just something that should not be there and throw this exception.
      //
      if (tokenNames == NULL)
      {
        message += " : Extraneous input...";
      }
      else
      {
        if (ex->expecting == ANTLR3_TOKEN_EOF)
        {
          message += " : Extraneous input - expected <EOF>\n";
        }
        else
        {
          message += str(boost::format(" : Extraneous input - expected %s ...\n") %tokenNames[ex->expecting]);
        }
      }
      break;

    case ANTLR3_MISSING_TOKEN_EXCEPTION:

      // Indicates that the recognizer detected that the token we just
      // hit would be valid syntactically if preceeded by a particular
      // token. Perhaps a missing ';' at line end or a missing ',' in an
      // expression list, and such like.
      //
      if (tokenNames == NULL)
      {
        message += str(boost::format(" : Missing token (%d)...\n") %ex->expecting);
      }
      else
      {
        if (ex->expecting == ANTLR3_TOKEN_EOF)
        {
          message += " : Missing <EOF>\n";
        }
        else
        {
          message += str(boost::format(" : Missing %s \n") %tokenNames[ex->expecting]);
        }
      }
      break;

    case ANTLR3_RECOGNITION_EXCEPTION:

      // Indicates that the recognizer received a token
      // in the input that was not predicted. This is the basic exception type
      // from which all others are derived. So we assume it was a syntax error.
      // You may get this if there are not more tokens and more are needed
      // to complete a parse for instance.
      //
      message += " : syntax error...\n";
      break;

    case ANTLR3_MISMATCHED_TOKEN_EXCEPTION:

      // We were expecting to see one thing and got another. This is the
      // most common error if we coudl not detect a missing or unwanted token.
      // Here you can spend your efforts to
      // derive more useful error messages based on the expected
      // token set and the last token and so on. The error following
      // bitmaps do a good job of reducing the set that we were looking
      // for down to something small. Knowing what you are parsing may be
      // able to allow you to be even more specific about an error.
      //
      if (tokenNames == NULL)
      {
        message += " : syntax error...\n";
      }
      else
      {
        if (ex->expecting == ANTLR3_TOKEN_EOF)
        {
          message += " : expected <EOF>\n";
        }
        else
        {
          message += str(boost::format(" : expected %s ...\n") %tokenNames[ex->expecting]);
        }
      }
      break;

    case ANTLR3_NO_VIABLE_ALT_EXCEPTION:

      // We could not pick any alt decision from the input given
      // so god knows what happened - however when you examine your grammar,
      // you should. It means that at the point where the current token occurred
      // that the DFA indicates nowhere to go from here.
      //
      message += " : cannot match to any predicted input...\n";

      break;

    case ANTLR3_MISMATCHED_SET_EXCEPTION:

      {
        ANTLR3_UINT32     count;
        ANTLR3_UINT32     bit;
        ANTLR3_UINT32     size;
        ANTLR3_UINT32     numbits;
        pANTLR3_BITSET    errBits;

        // This means we were able to deal with one of a set of
        // possible tokens at this point, but we did not see any
        // member of that set.
        //
        message += " : unexpected input...\n  expected one of : ";

        // What tokens could we have accepted at this point in the
        // parse?
        //
        count   = 0;
        errBits = antlr3BitsetLoad              (ex->expectingSet);
        numbits = errBits->numBits              (errBits);
        size    = errBits->size                 (errBits);

        if (size > 0)
        {
          // However many tokens we could have dealt with here, it is usually
          // not useful to print ALL of the set here. I arbitrarily chose 8
          // here, but you should do whatever makes sense for you of course.
          // No token number 0, so look for bit 1 and on.
          //
          for (bit = 1; bit < numbits && count < 8 && count < size; bit++)
          {
            // TODO: This doesn;t look right - should be asking if the bit is set!!
            //
            if (tokenNames[bit])
            {
              message += str(boost::format("%s%s") %(count > 0 ? ", " : "")
                                                   %tokenNames[bit]);
              count++;
            }
          }
          message += "\n";
        }
        else
        {
          message += "Actually dude, we didn't seem to be expecting anything here, or at least\n";
          message += "I could not work out what I was expecting, like so many of us these days!\n";
        }
      }
      break;

    case ANTLR3_EARLY_EXIT_EXCEPTION:

      // We entered a loop requiring a number of token sequences
      // but found a token that ended that sequence earlier than
      // we should have done.
      //
      message += " : missing elements...\n";
      break;

    default:

      // We don't handle any other exceptions here, but you can
      // if you wish. If we get an exception that hits this point
      // then we are just going to report what we know about the
      // token.
      //
      message += " : syntax not recognized...\n";
      break;
    }

    // Here you have the token that was in error which if this is
    // the standard implementation will tell you the line and offset
    // and also record the address of the start of the line in the
    // input stream. You could therefore print the source line and so on.
    // Generally though, I would expect that your lexer/parser will keep
    // its own map of lines and source pointers or whatever as there
    // are a lot of specific things you need to know about the input
    // to do something like that.
    // Here is where you do it though :-).
    //

    // How we determine the next piece is dependent on which thing raised the
    // error.
    //
    switch (recognizer->type)
    {
    case ANTLR3_TYPE_PARSER:

      pANTLR3_PARSER parser;
      pcifParser generated;

      parser = (pANTLR3_PARSER) (recognizer->super);
      generated  = (pcifParser)(parser->super);

      generated->errors->push_back(message);
      break;

    default:
      std::cerr << message;
      break;
    }

} // displayRecognitionError()


/** Default lexer error handler (works for 8 bit streams only!!!)
 */
static void
lexer_displayRecognitionError (pANTLR3_BASE_RECOGNIZER recognizer, pANTLR3_UINT8 * tokenNames)
{
  pANTLR3_LEXER               lexer;
  pANTLR3_EXCEPTION           ex;
  pANTLR3_STRING              ftext;

  lexer = (pANTLR3_LEXER)(recognizer->super);
  ex = lexer->rec->state->exception;


  std::string message;

  // See if there is a 'filename' we can use
  //
  if (ex->name == NULL)
  {
    message += "-unknown source-(";
  }
  else
  {
    ftext = ex->streamName->to8(ex->streamName);
    message += str(boost::format("%s(") %ftext->chars);
  }

  message += str(boost::format("line %d) ") %recognizer->state->exception->line);
  message += str(boost::format(": lexer error %d :\n\t%s at offset %d, ")
                                              %ex->type
                                              %((pANTLR3_UINT8)(ex->message))
                                              %(ex->charPositionInLine+1))
  ;
  {
    ANTLR3_INT32 width;

    width = ANTLR3_UINT32_CAST(( (pANTLR3_UINT8)(lexer->input->data) + (lexer->input->size(lexer->input) )) - (pANTLR3_UINT8)(ex->index));

    if (width >= 1)
    {
      if (isprint(ex->c))
      {
        message += str(boost::format("near '%c' :\n") %ex->c);
      }
      else
      {
        message += str(boost::format("near char(%#02X) :\n") %(ANTLR3_UINT8)(ex->c));
      }
      boost::format fmt(str(boost::format("\t%%.%is\n") % std::min(width, 20)));
      message += str(fmt % (pANTLR3_UINT8)(ex->index));
    }
    else
    {
      message += "(end of input).\n\t This indicates a poorly specified lexer RULE\n\t or unterminated input element such as: \"STRING[\"]\n";
      message += str(boost::format("\t The lexer was matching from line %d, offset %d, which\n\t ")
                                  %((ANTLR3_UINT32)(lexer->rec->state->tokenStartLine))
                                  %((ANTLR3_UINT32)(lexer->rec->state->tokenStartCharPositionInLine))
                                  );
      width = ANTLR3_UINT32_CAST(((pANTLR3_UINT8)(lexer->input->data)+(lexer->input->size(lexer->input))) - (pANTLR3_UINT8)(lexer->rec->state->tokenStartCharIndex));

      if (width >= 1)
      {
        boost::format fmt(str(boost::format("looks like this:\n\t\t%%.%is\n")
                              % std::min(width, 20)));
        message += str(
          fmt % (pANTLR3_UINT8)(lexer->rec->state->tokenStartCharIndex));
      }
      else
      {
        message += "is also the end of the line, so you must check your lexer rules\n";
      }
    }
  }

  pcifLexer generated;
  lexer = (pANTLR3_LEXER) (recognizer->super);
  generated = (pcifLexer)(lexer->super);

  generated->errors->push_back(message);
}


}} // namespace iotbx::cif

#endif // GUARD
