/** CIF Version 1.1 Working specification grammar

Translated from the grammar defined at

http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#bnf

This is a prototype version of the parser using a Python target,
and is hence much slower than the compiled C++ target version.

Richard Gildea
April 2010
*/

grammar cifProto;

options {
    language=Python;
}

@header {
import iotbx.cif
from cctbx.array_family import flex
}

@members {
paraphrases = []
def getErrorMessage(self, e, tokenNames):
  msg = Parser.getErrorMessage(self, e, tokenNames)
  if len(self.paraphrases) > 0:
    paraphrase = self.paraphrases[-1]
    msg += " " + paraphrase
  return msg
}

/*------------------------------------------------------------------
 * PARSER RULES
 *------------------------------------------------------------------*/

// The start rule
parse[builder]
@init { self.builder = builder }
	: cif EOF ;
/*------------------------------------------------------------------
 * BASIC STRUCTURE OF A CIF
 *------------------------------------------------------------------*/

cif
	:	(COMMENTS)? (WHITESPACE)* ( data_block ( WHITESPACE* data_block )* (WHITESPACE)* )?
	;

loop_body
@init { self.curr_loop_values = flex.std_string()}
	:	v1=value
{ self.curr_loop_values.append(str($v1.text)) }
	      ( WHITESPACE+
	        v2=value
{ self.curr_loop_values.append(str($v2.text)) }
	       )*
	       { self.paraphrases.pop() } // remove "in loop_" from stack
	;

save_frame
	:	SAVE_FRAME_HEADING
{ self.builder.start_save_frame($SAVE_FRAME_HEADING.text) }
	 ( WHITESPACE+ data_items )+ WHITESPACE+ SAVE
{ self.builder.end_save_frame() }
	;

data_items
@init  { self.paraphrases.append("in data items") }
@after { self.paraphrases.pop() }
	: 	( TAG WHITESPACE value
{ self.builder.add_data_item($TAG.text, $value.text) }
		)
	      | ( loop_header WHITESPACE* loop_body
{
try:
  self.builder.add_loop(self.curr_loop_headers, data=self.curr_loop_values)
except AssertionError, e:
  print e
}
		)
	;

data_block
	:	DATA_BLOCK_HEADING
{ self.builder.add_data_block($DATA_BLOCK_HEADING.text) }
	      ( WHITESPACE+ ( data_items | save_frame ) )*
	;

loop_header
@init { self.curr_loop_headers = flex.std_string()}
	:	LOOP_ { self.paraphrases.append("in loop_") }
		( WHITESPACE+ TAG
{ self.curr_loop_headers.append(str($TAG.text)) }
		 )+ WHITESPACE
		;

/*------------------------------------------------------------------
 * TAGS AND VALUES
 *------------------------------------------------------------------*/

inapplicable
	:	'.' ;

unknown	:	'?' ;

value
@init  { self.paraphrases.append("in value") }
@after { self.paraphrases.pop() }
 	:	inapplicable | unknown | '-' | char_string  | numeric| text_field
 	;
 	catch [RecognitionException, re] {
	  raise re
	}

unsigned_integer
	:	(DIGIT)+ ;

integer	: 	( '+' | '-' )? unsigned_integer ;

float_	: 	integer EXPONENT | ( ( '+' | '-' )? ( (DIGIT)* '.' unsigned_integer) | (DIGIT)+ '.' ) (EXPONENT)? ;

number  :	integer | float_ ;

numeric	:	number | ( number '(' (DIGIT)+ ')' ) ;

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

TAG	:	'_' ( 'A'..'Z' | 'a'..'z' ) (NON_BLANK_CHAR_)* ;

/*------------------------------------------------------------------
 * CHARACTER STRINGS AND FIELDS
 *------------------------------------------------------------------*/

SEMI_COLON_TEXT_FIELD
	:	';'
		( ( ANY_PRINT_CHAR | SINGLE_QUOTE | DOUBLE_QUOTE )* EOL
		( (TEXT_LEAD_CHAR ( ANY_PRINT_CHAR | SINGLE_QUOTE | DOUBLE_QUOTE )* )? EOL)* )
		';'
	;

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

VERSION	:	'#\\#CIF_' (DIGIT)+ '.' (DIGIT)+ ;

DATA_BLOCK_HEADING
	:	DATA_ (NON_BLANK_CHAR)+ ;

SAVE_FRAME_HEADING
	:	SAVE_ (NON_BLANK_CHAR)+ ;

SAVE	:	SAVE_ ;

// apparently a single quoted string such as 'a dog's life' is legal...
fragment SINGLE_QUOTED_STRING
	:	SINGLE_QUOTE
		( ( (SINGLE_QUOTE NON_BLANK_CHAR_)=>SINGLE_QUOTE ) | ANY_PRINT_CHAR | DOUBLE_QUOTE )*
		SINGLE_QUOTE
		{ self.setText(self.getText()[1:-1]) }
	;

fragment DOUBLE_QUOTED_STRING
	:	DOUBLE_QUOTE
		( ( (DOUBLE_QUOTE NON_BLANK_CHAR_)=>DOUBLE_QUOTE ) | ANY_PRINT_CHAR | SINGLE_QUOTE )*
	        DOUBLE_QUOTE
		{ self.setText(self.getText()[1:-1]) }
	;

/*------------------------------------------------------------------
 * NUMERICS
 *------------------------------------------------------------------*/

DIGIT	: '0'..'9' ;

EXPONENT: 	( ( 'e' | 'E') | ( 'e' | 'E')( '+' | '-' ) ) (DIGIT)+ ;

// UNQUOTED_STRING must be defined after the digits if we want to catch digits first
fragment UNQUOTED_STRING

	:	( ORDINARY_CHAR | ';' ) (NON_BLANK_CHAR_)* ;

CHAR_STRING
	:	SINGLE_QUOTED_STRING | DOUBLE_QUOTED_STRING | UNQUOTED_STRING;

/*------------------------------------------------------------------
 * WHITE SPACE AND COMMENTS
 *------------------------------------------------------------------*/

COMMENTS
	:	( ( '#' (ANY_PRINT_CHAR | SINGLE_QUOTE | DOUBLE_QUOTE )*
	          ( EOL | { self.input.LA(1) == EOF }? ) )+ )
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
