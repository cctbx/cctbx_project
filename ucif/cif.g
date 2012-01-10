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
#include <ucif/builder.h>
}

@includes{
#include <string>
#include <ucif/builder.h>
}

@parser::context
{
  ucif::array_wrapper_base* errors;
}

@lexer::context
{
  ucif::array_wrapper_base* errors;
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
  ANTLR3_MARKER start = token->start;
  ANTLR3_MARKER stop = token->stop;
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
parse[ucif::builder_base* builder_, bool strict_]
scope {
  ucif::builder_base* builder;
  bool strict;
}
@init {
  $parse::builder = builder_;
  $parse::strict = strict_;
}
  : cif (EOF | '\u001a' /*Ctrl-Z*/) ;
/*------------------------------------------------------------------
  * BASIC STRUCTURE OF A CIF
  *------------------------------------------------------------------*/

cif
  :	COMMENTS? data_block*
  ;

loop_body
  :	v1=value
{ ($data_items::curr_loop_values)->push_back(to_std_string($v1.start)); }
  ( v2=value
{ ($data_items::curr_loop_values)->push_back(to_std_string($v2.start)); }
  )*
;

save_frame
  :	SAVE_FRAME_HEADING
{ ($parse::builder)->start_save_frame(to_std_string($SAVE_FRAME_HEADING)); }
  ( data_items )+ SAVE
{ ($parse::builder)->end_save_frame(); }
  ;

data_items
scope {
  ucif::array_wrapper_base* curr_loop_values;
  ucif::array_wrapper_base* curr_loop_headers;
}
@init {
  $data_items::curr_loop_values = ($parse::builder)->new_array();
  $data_items::curr_loop_headers = ($parse::builder)->new_array();
}
@after {
  delete $data_items::curr_loop_values;
  delete $data_items::curr_loop_headers;
}
  :	TAG value
{
  ($parse::builder)->add_data_item(
    to_std_string($TAG),
    to_std_string($value.start));
}
  | loop_header loop_body
{
  ucif::array_wrapper_base* values = $data_items::curr_loop_values;
  int n_cols = $data_items::curr_loop_headers->size();
  if (values->size() \% n_cols != 0) {
    std::string msg = "Wrong number of data items for loop containing ";
    msg += (*$data_items::curr_loop_headers)[0];
    CTX->errors->push_back(msg);
  }
  else {
    ($parse::builder)->add_loop(*$data_items::curr_loop_headers, *values);
  }
}
  ;

data_block
  :	( DATA_BLOCK_HEADING
{ ($parse::builder)->add_data_block(to_std_string($DATA_BLOCK_HEADING)); }
  ( ( data_items | save_frame ) )*
  )
  | ( {!$parse::strict}?=>GLOBAL_ ( ( data_items | save_frame ) )* ) // global blocks are ignored
;

loop_header
  :	LOOP_ ( TAG
{ ($data_items::curr_loop_headers)->push_back(to_std_string($TAG)); }
  )+
;

/*------------------------------------------------------------------
  * TAGS AND VALUES
  *------------------------------------------------------------------*/


value
  :	INAPPLICABLE | UNKNOWN | NUMERIC | char_string | text_field
  ;

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

INAPPLICABLE
  : '.' ;

UNKNOWN : '?' ;

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
    { SKIP(); }
  ;
