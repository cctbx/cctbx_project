/** CIF Version 1.1 Working specification grammar

Translated from the grammar defined at

http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#bnf

A compiled version of the parser, with C language target, but contains
C++ code in the actions, therefore the output files must be renamed to *.cpp

Richard Gildea
April 2010
*/

tree grammar cifWalker;

options {
  tokenVocab=cif;
  ASTLabelType=pANTLR3_BASE_TREE;
  language=C;
}

@includes{
#include <string>
#include <ucif/builder.h>
}

@context
{
  ucif::array_wrapper_base* errors;
}

@members {
  std::string to_std_string(pANTLR3_BASE_TREE node) {
    pANTLR3_COMMON_TOKEN token = ((pANTLR3_COMMON_TREE)(node->super))->token;
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
parse[ucif::builder_base* builder_]
scope { ucif::builder_base* builder; }
@init { $parse::builder = builder_; }

  : cif ;
/*------------------------------------------------------------------
  * BASIC STRUCTURE OF A CIF
  *------------------------------------------------------------------*/

cif
  :	data_block*
  ;

loop_body
  :
    ( value
{ ($data_items::curr_loop_values)->push_back(to_std_string($value.start)); }
  )+
;

save_frame
  :	^(SAVE SAVE_FRAME_HEADING
{ ($parse::builder)->start_save_frame(to_std_string($SAVE_FRAME_HEADING)); }
  ( data_items )+
{ ($parse::builder)->end_save_frame(); }
  )
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
  :	 ^(TAG_VALUE_PAIR TAG value)
{
  ($parse::builder)->add_data_item(
    to_std_string($TAG),
    to_std_string($value.start));
}
  | ^(LOOP loop_header loop_body)
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

data_block_heading
  : DATA_BLOCK_HEADING | GLOBAL_
  ;

data_block
  :	^(DATA_BLOCK data_block_heading
{ ($parse::builder)->add_data_block(to_std_string($data_block_heading.start)); }
  data_items* save_frame*
)
;

loop_header
  :	^(LOOP_ ( TAG
{ ($data_items::curr_loop_headers)->push_back(to_std_string($TAG)); }
  )+ )
;

/*------------------------------------------------------------------
  * TAGS AND VALUES
  *------------------------------------------------------------------*/

inapplicable
  : '.' ;

unknown : '?' ;

value
  : INAPPLICABLE | UNKNOWN | NUMERIC | char_string | text_field
  ;

char_string
  : CHAR_STRING ;

text_field
  : SEMI_COLON_TEXT_FIELD ;
