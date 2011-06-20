/** CIF Version 1.1 Working specification grammar

Translated from the grammar defined at

http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#bnf

A compiled version of the parser, with C language target, but contains
C++ code in the actions, therefore the output files must be renamed to *.cpp

Richard Gildea
April 2010
*/

grammar cif;

options {
    language=C;
}

@lexer::includes{
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

}

@includes{
#include <string>
#include <vector>

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

#include <boost/python/object.hpp>
}

@parser::context
{
    scitbx::af::shared<std::string> *errors;
}

@lexer::context
{
    scitbx::af::shared<std::string> *errors;
}

@parser::apifuncs
{
  PARSER->super = (void *)ctx;
}

@lexer::apifuncs
{
  LEXER->super = (void *)ctx;
}
@members {
std::string to_std_string(pANTLR3_COMMON_TOKEN token) {
	ANTLR3_MARKER start = token->getStartIndex(token);
	ANTLR3_MARKER stop = token->getStopIndex(token);
	std::string str((const char*)start, stop-start+1);
	if ((str[0] == '\'' && str[str.size()-1] == '\'') ||
	    (str[0] == '"' && str[str.size()-1] == '"'))
	  { str = str.substr(1, str.size()-2); }
	return str;
}

}
/*------------------------------------------------------------------
 * PARSER RULES
 *------------------------------------------------------------------*/

// The start rule
parse[boost::python::object & builder_]
scope { boost::python::object *builder; }
@init { $parse::builder = new boost::python::object(builder_); }
@after { delete $parse::builder; }

	: cif (EOF | '\u001a' /*Ctrl-Z*/) ;
/*------------------------------------------------------------------
 * BASIC STRUCTURE OF A CIF
 *------------------------------------------------------------------*/

cif
	:	(COMMENTS)? (WHITESPACE)* ( data_block ( WHITESPACE* data_block )* (WHITESPACE)* )?
	;

loop_body
	:	v1=value
{ ($data_items::curr_loop_values)->push_back(to_std_string($v1.start)); }
	      ( WHITESPACE+
	        v2=value
{ ($data_items::curr_loop_values)->push_back(to_std_string($v2.start)); }
	       )*
	;

save_frame
	:	SAVE_FRAME_HEADING
{ ($parse::builder)->attr("start_save_frame")(to_std_string($SAVE_FRAME_HEADING)); }
	      ( WHITESPACE+ data_items )+ WHITESPACE+ SAVE
{ ($parse::builder)->attr("end_save_frame")(); }
	;

data_items
scope { scitbx::af::shared<std::string> *curr_loop_values;
        scitbx::af::shared<std::string> *curr_loop_headers;
}
@init { $data_items::curr_loop_values = new scitbx::af::shared<std::string>();
	$data_items::curr_loop_headers = new scitbx::af::shared<std::string>();
}
@after { delete $data_items::curr_loop_values;
 	 delete $data_items::curr_loop_headers;
}
	:	TAG WHITESPACE* value
{
  ($parse::builder)->attr("add_data_item")(
  to_std_string($TAG),
  to_std_string($value.start));
}
	      | loop_header WHITESPACE* loop_body
{
  scitbx::af::shared<std::string> &values = *($data_items::curr_loop_values);
  int n_cols = $data_items::curr_loop_headers->size();
  if (values.size() \% n_cols != 0) {
    std::string msg = "Wrong number of data items for loop containing ";
    msg += (*$data_items::curr_loop_headers)[0];
    CTX->errors->push_back(msg);
  }
  else {
    try {
      ($parse::builder)->attr("add_loop")($data_items::curr_loop_headers, values);
    }
    catch (boost::python::error_already_set&) {
      PyErr_Print();
      PyErr_Clear();
    }
  }
}
	;

data_block
	:	DATA_BLOCK_HEADING
{ ($parse::builder)->attr("add_data_block")(to_std_string($DATA_BLOCK_HEADING)); }
	      ( WHITESPACE+ ( data_items | save_frame ) )*
	;

loop_header
	:	LOOP_ ( WHITESPACE+ TAG
{ ($data_items::curr_loop_headers)->push_back(to_std_string($TAG)); }
		 )+ WHITESPACE
	;

/*------------------------------------------------------------------
 * TAGS AND VALUES
 *------------------------------------------------------------------*/

inapplicable
	:	'.' ;

unknown	:	'?' ;

value 	:	inapplicable | unknown | '-' | char_string  | NUMERIC | text_field
	;
	catch [RecognitionException re] {
	}

char_string
	:	CHAR_STRING ;

text_field
	:	SEMI_COLON_TEXT_FIELD ;




/*------------------------------------------------------------------
 * LEXER RULES
 *------------------------------------------------------------------*/

/*------------------------------------------------------------------
 * CHARACTER SETS
 *------------------------------------------------------------------*/

fragment EOL
	:	( '\n' | '\r' | '\r\n' ) ;

fragment DOUBLE_QUOTE
	:	'"' ;

fragment SINGLE_QUOTE
	:	'\'' ;

fragment ORDINARY_CHAR
	: 	'!' | '%' | '&' | '(' | ')' | '*' | '+' | ',' | '-' | '.' | '/' |
	( '0'.. '9' ) | ':' | '<' | '=' | '>' | '?' | '@' | ('A'..'Z') | ('a'..'z') |
	'\\' | '^' | '`' | '{' | '|' | '}' | '~'
	
	;

fragment NON_BLANK_CHAR_
	:	ORDINARY_CHAR | DOUBLE_QUOTE | SINGLE_QUOTE | '#' | '$' | '_' | '[' | ']' | ';' ;

fragment TEXT_LEAD_CHAR
	:	ORDINARY_CHAR | DOUBLE_QUOTE | SINGLE_QUOTE | '#' | '$' | '_' | '[' | ']' | ' ' | '\t' ;

fragment ANY_PRINT_CHAR
	:	ORDINARY_CHAR | '#' | '$' | '_' | '[' | ']' | ' ' | '\t' | ';'
// temporarily disable quotes from any print char until
// I figure out how to do it properly...
//	      | DOUBLE_QUOTE | SINGLE_QUOTE
	;

TAG	:	'_' (NON_BLANK_CHAR_)+ ;


/*------------------------------------------------------------------
 * RESERVED WORDS - define these after semicolon text field
 *------------------------------------------------------------------*/

fragment
DATA_	:	( 'D' | 'd' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'A' | 'a' ) '_' ;

fragment
SAVE_	:	( 'S' | 's' ) ( 'A' | 'a' ) ( 'V' | 'v' ) ( 'E' | 'e' ) '_' ;

LOOP_ 	:	( 'L' | 'l' ) ( 'O' | 'o' ) ( 'O' | 'o' ) ( 'P' | 'p' ) '_' ;

GLOBAL_ :	( 'G' | 'g' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'B' | 'b' ) ( 'A' | 'a' ) ( 'L' | 'l' ) '_' ;

STOP_	:	( 'S' | 's' ) ( 'T' | 't' ) ( 'O' | 'o' ) ( 'P' | 'p' ) '_' ;

/*------------------------------------------------------------------
 * SPECIAL KEY WORDS
 *------------------------------------------------------------------*/

//VERSION	:	'#\\#CIF_' (DIGIT)+ '.' (DIGIT)+ ;

DATA_BLOCK_HEADING
	:	DATA_ (NON_BLANK_CHAR)+ ;

SAVE_FRAME_HEADING
	:	SAVE_ (NON_BLANK_CHAR)+ ;

SAVE	:	SAVE_ ;

/*------------------------------------------------------------------
 * NUMERICS
 *------------------------------------------------------------------*/

fragment DIGIT	: '0'..'9' ;

fragment EXPONENT
        : 	( ( 'e' | 'E') | ( 'e' | 'E')( '+' | '-' ) ) (DIGIT)+ ;

fragment INTEGER
	: 	( '+' | '-' )? (DIGIT)+ ;

fragment FLOAT
	: 	INTEGER EXPONENT | ( ( '+' | '-' )? ( (DIGIT)* '.' (DIGIT)+) | (DIGIT)+ '.' ) (EXPONENT)? ;

fragment UNSIGNED_INTEGER
	:	(DIGIT)+ ;

fragment NUMBER
	:	INTEGER | FLOAT;

NUMERIC	:	NUMBER | ( NUMBER '(' (UNSIGNED_INTEGER)+ ')' ) ;

/*------------------------------------------------------------------
 * CHARACTER STRINGS AND FIELDS
 *------------------------------------------------------------------*/

fragment UNQUOTED_STRING

	:	(({ GETCHARPOSITIONINLINE() > 0 }?=> ';')
		  | ORDINARY_CHAR ) (NON_BLANK_CHAR_)* ;

// apparently a single quoted string such as 'a dog's life' is legal...
fragment SINGLE_QUOTED_STRING
	:	SINGLE_QUOTE
		( ( (SINGLE_QUOTE NON_BLANK_CHAR_)=>SINGLE_QUOTE ) | ANY_PRINT_CHAR | DOUBLE_QUOTE )*
		SINGLE_QUOTE
	;

fragment DOUBLE_QUOTED_STRING
	:	DOUBLE_QUOTE
		( ( (DOUBLE_QUOTE NON_BLANK_CHAR_)=>DOUBLE_QUOTE ) | ANY_PRINT_CHAR | SINGLE_QUOTE )*
	        DOUBLE_QUOTE
	;

CHAR_STRING
	:	UNQUOTED_STRING | SINGLE_QUOTED_STRING | DOUBLE_QUOTED_STRING;

SEMI_COLON_TEXT_FIELD
	:	( { GETCHARPOSITIONINLINE() == 0 }?=> ';')
		( ( ANY_PRINT_CHAR | SINGLE_QUOTE | DOUBLE_QUOTE )* EOL
		( (TEXT_LEAD_CHAR ( ANY_PRINT_CHAR | SINGLE_QUOTE | DOUBLE_QUOTE )* )? EOL)* )
		
		{ $start += 1; EMIT(); } // strip semicolons
		
		';'
	;

/*------------------------------------------------------------------
 * WHITE SPACE AND COMMENTS
 *------------------------------------------------------------------*/

COMMENTS
	:	( ( '#' (ANY_PRINT_CHAR | SINGLE_QUOTE | DOUBLE_QUOTE )*
  		  ( EOL | { LA(1) == EOF }? ) )+ )
	        { $channel = HIDDEN; }
	;



//TOKENIZED_COMMENTS
//	:	( ' ' | '\t' | EOL )+ COMMENTS_
//	        { $channel = HIDDEN; }
//	;

// Redefine this as non-fragment so can be seen by the parser
NON_BLANK_CHAR
	:	NON_BLANK_CHAR_ ;

WHITESPACE
	: 	( '\t' | ' ' | EOL | '\u000C' )+
		//{ $channel = HIDDEN; }
	;
