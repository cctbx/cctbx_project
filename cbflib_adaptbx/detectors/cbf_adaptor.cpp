#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <fstream>
#include <iostream>
#include <exception>
#include <scitbx/array_family/flex_types.h>
#include <examples/img.h>
#include <cbflib_adaptbx/detectors/cbf_adaptor.h>

namespace af = scitbx::af;
namespace ide = iotbx::detectors;

bool
ide::CBFAdaptor::file_is_transposed(){
    std::string elem,precedence;
    const char *temp_value;
    cbf_failnez ( cbf_find_category(cbf_h,"array_structure_list") );
    elem = "ELEMENT_X";

    cbf_failnez ( cbf_find_column(cbf_h,"axis_set_id") );
    cbf_failnez ( cbf_rewind_row(cbf_h) );
    cbf_failnez ( cbf_find_row(cbf_h,elem.c_str()) );
    cbf_failnez ( cbf_rewind_column(cbf_h) );
    cbf_failnez ( cbf_find_column(cbf_h,"precedence") );
    cbf_failnez ( cbf_get_value(cbf_h,&temp_value) );
    precedence = std::string(temp_value);

    if (precedence=="1") return false;
    if (precedence=="2") return true;
    throw iotbx::detectors::Error ("Unable to determine precedence of ELEMENT_X");
}

std::string
ide::CBFAdaptor::raster_description(){
    std::string precedence,direction;
    std::string result="";
    const char *temp_value;
    cbf_failnez ( cbf_find_category(cbf_h,"array_structure_list") );
    typedef std::vector<std::string> el_lst;
    el_lst elements;
    elements.push_back("ELEMENT_X");elements.push_back("ELEMENT_Y");

    for (el_lst::const_iterator e = elements.begin(); e!=elements.end(); ++e) {
      cbf_failnez ( cbf_find_column(cbf_h,"axis_set_id") );
      cbf_failnez ( cbf_rewind_row(cbf_h) );
      cbf_failnez ( cbf_find_row(cbf_h,e->c_str()) );
      cbf_failnez ( cbf_rewind_column(cbf_h) );
      cbf_failnez ( cbf_find_column(cbf_h,"direction") );
      cbf_failnez ( cbf_get_value(cbf_h,&temp_value) );
      direction = std::string(temp_value);
      cbf_failnez ( cbf_rewind_column(cbf_h) );
      cbf_failnez ( cbf_find_column(cbf_h,"precedence") );
      cbf_failnez ( cbf_get_value(cbf_h,&temp_value) );
      precedence = std::string(temp_value);
      result = result+std::string(*e)+" "+direction+" precedence="+precedence+'\n';
    }

    return result;
}
