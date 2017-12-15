#ifndef UCIF_UTILS_H
#define UCIF_UTILS_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <ucif/cifParser.h>

namespace ucif {

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

  std::ostringstream message;
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
        message << "-end of input-(";
      }
      else
      {
        message << "-unknown source-(";
      }
    }
    else
    {
      ftext = ex->streamName->to8(ex->streamName);
      message << ftext->chars << "(";
    }

    // Next comes the line number
    //

    message << "line " << recognizer->state->exception->line << ")";
    message << " : error " << recognizer->state->exception->type << " : "
      << (pANTLR3_UINT8) (recognizer->state->exception->message);


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
      ttext       =  theToken->getText(theToken);

      message << ", at offset " << recognizer->state->exception->charPositionInLine;
      if  (theToken != NULL)
      {
        if (theToken->type == ANTLR3_TOKEN_EOF)
        {
          message << ", at <EOF>";
        }
        else
        {
          // Guard against null text in a token
          //
          message << "\n    near "
            << (ttext == NULL ? (pANTLR3_UINT8)"<no text for the token>" : ttext->chars)
            << "\n    ";
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
        message << ", at offset " << theBaseTree->getCharPositionInLine(theBaseTree);
        message << ", near " << ttext->chars;
      }
      break;

    default:

      message << "Base recognizer function displayRecognitionError called by unknown parser type - provide override for this function\n";
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
        message << " : Extraneous input...";
      }
      else
      {
        if (ex->expecting == ANTLR3_TOKEN_EOF)
        {
          message << " : Extraneous input - expected <EOF>\n";
        }
        else
        {
          message << " : Extraneous input - expected " << tokenNames[ex->expecting] << " ...\n";
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
        message << " : Missing token (" << ex->expecting << ")...\n";
      }
      else
      {
        if (ex->expecting == ANTLR3_TOKEN_EOF)
        {
          message <<" : Missing <EOF>\n";
        }
        else
        {
          message << " : Missing " << tokenNames[ex->expecting] << " \n";
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
      message << " : syntax error...\n";
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
        message << " : syntax error...\n";
      }
      else
      {
        if (ex->expecting == ANTLR3_TOKEN_EOF)
        {
          message << " : expected <EOF>\n";
        }
        else
        {
          message << " : expected " << tokenNames[ex->expecting] << " ...\n";
        }
      }
      break;

    case ANTLR3_NO_VIABLE_ALT_EXCEPTION:

      // We could not pick any alt decision from the input given
      // so god knows what happened - however when you examine your grammar,
      // you should. It means that at the point where the current token occurred
      // that the DFA indicates nowhere to go from here.
      //
      message << " : cannot match to any predicted input...\n";

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
        message << " : unexpected input...\n  expected one of : ";

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
              message << (count > 0 ? ", " : "") << tokenNames[bit];
              count++;
            }
          }
          message << "\n";
        }
        else
        {
          message << "Unable to process inputs\n";
        }
      }
      break;

    case ANTLR3_EARLY_EXIT_EXCEPTION:

      // We entered a loop requiring a number of token sequences
      // but found a token that ended that sequence earlier than
      // we should have done.
      //
      message << " : missing elements...\n";
      break;

    default:

      // We don't handle any other exceptions here, but you can
      // if you wish. If we get an exception that hits this point
      // then we are just going to report what we know about the
      // token.
      //
      message << " : syntax not recognized...\n";
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

      generated->errors->push_back(message.str());
      break;

    default:
      std::cerr << message.str();
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


  std::ostringstream message;

  // See if there is a 'filename' we can use
  //
  if (ex->name == NULL)
  {
    message << "-unknown source-(";
  }
  else
  {
    ftext = ex->streamName->to8(ex->streamName);
    message << ftext->chars << "(";
  }

  message << "line " << recognizer->state->exception->line << ") ";
  message << ": lexer error " << ex->type << " :\n\t"
    << (pANTLR3_UINT8)(ex->message) << " at offset " << (ex->charPositionInLine+1) << ", "
  ;
  {
    ANTLR3_INT32 width;

    width = ANTLR3_UINT32_CAST(( (pANTLR3_UINT8)(lexer->input->data) + (lexer->input->size(lexer->input) )) - (pANTLR3_UINT8)(ex->index));

    if (width >= 1)
    {
      if (isprint(ex->c))
      {
        message << "near '" << ex->c << "' :\n";
      }
      else
      {
        message << "near char(" << std::hex << (ANTLR3_UINT8)(ex->c) << ") :\n";
      }
      message << "\t" << std::string(reinterpret_cast<const char*>(ex->index), std::min(width, 20)) << "\n";
    }
    else
    {
      message << "(end of input).\n\t This indicates a poorly specified lexer RULE\n\t or unterminated input element such as: \"STRING[\"]\n";
      message << "\t The lexer was matching from line "
              << (ANTLR3_UINT32)(lexer->rec->state->tokenStartLine)
              << ", offset "
              << (ANTLR3_UINT32)(lexer->rec->state->tokenStartCharPositionInLine)
              << ", which\n\t ";
      width = ANTLR3_UINT32_CAST(((pANTLR3_UINT8)(lexer->input->data)+(lexer->input->size(lexer->input))) - (pANTLR3_UINT8)(lexer->rec->state->tokenStartCharIndex));

      if (width >= 1)
      {
        message << "looks like this:\n\t\t"
                << std::string(reinterpret_cast<const char*>(lexer->rec->state->tokenStartCharIndex), std::min(width, 20))
                << "\n";
      }
      else
      {
        message << "is also the end of the line, so you must check your lexer rules\n";
      }
    }
  }

  pcifLexer generated;
  lexer = (pANTLR3_LEXER) (recognizer->super);
  generated = (pcifLexer)(lexer->super);

  generated->errors->push_back(message.str());
}


} // namespace ucif

#endif // GUARD
