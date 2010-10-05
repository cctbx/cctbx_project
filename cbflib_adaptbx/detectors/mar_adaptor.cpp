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
#include <cbflib_adaptbx/detectors/basic.h>
#include <cbflib_adaptbx/detectors/mar_adaptor.h>

namespace af = scitbx::af;

namespace iotbx {
  namespace detectors {

AnyImgAdaptor::AnyImgAdaptor(const std::string& filename):
    filename(filename),
    img(img_make_handle()),
    read_header_ok(false),
    read_data_ok(false){
    img_set_tags (img, 0);
    img_set_dimensions (img, 0, 0);
  }
void AnyImgAdaptor::read_header() {
  if (read_header_ok) {return;}
  if (!img) throw Error("img_BAD_ARGUMENT");
  private_file = std::fopen(filename.c_str(),"rb");
  if (!private_file) throw Error("img_BAD_OPEN");
  read_header_ok = header_read_specific(img,private_file,img_org_data);
  data_read_position = std::ftell(private_file);
  std::fclose(private_file);
  if (!read_header_ok) throw Error("Bad header read");
}
void AnyImgAdaptor::read_data() {
  if (read_data_ok) {return;}
  read_header();
  if (!img) throw Error("img_BAD_ARGUMENT from read_data");
  private_file = std::fopen(filename.c_str(),"rb");
  if (!private_file) throw Error("img_BAD_OPEN from read_data");
  std::fseek(private_file,data_read_position,0);
  read_data_ok = data_read_specific(img,private_file,img_org_data);
  std::fclose(private_file);
  if (!read_data_ok) throw Error("Bad data read");
}
std::string AnyImgAdaptor::get_field(const std::string& field_name){
  try {
    const char * found_value;
    found_value = img_get_field(img,field_name.c_str());
    if (!found_value) throw Error("No value found for field "+field_name);
    char stack_value [64];
    std::strncpy (stack_value, found_value, 63);
    stack_value [63] = 0;
    std::string return_candidate(stack_value);
    //while (*(return_candidate.rbegin())==' ') {
    //  return_candidate.erase(return_candidate.size()-1,1);
    //}
    return return_candidate;
  } catch (...) {
    throw Error("problem parsing the field "+field_name);
  }
}
af::flex_int AnyImgAdaptor::rawdata() {
  read_data();
  af::flex_int z(af::flex_grid<>((long)size1(),(long)size2()));
  int* begin = z.begin();
  std::size_t sz = z.size();
  std::size_t side = static_cast<std::size_t>(
    std::sqrt(static_cast<double>(sz)));
  int conv = 0;
  //for (std::size_t i = 0; i < sz; i++) {
  //    begin[i] = img->image[i];
  //}
  if (1){
  for (std::size_t i = 0; i < side; i++) {
    for (std::size_t j = 0; j < side; j++) {
      if (conv==0)
        {begin[i*side+j] = img->image[i*side+j]; }
      else if (conv==1)
        {begin[i*side+j] = img->image[j*side+i]; }
      else if (conv==2)
        {begin[i*side+j] = img->image[(side-1-j)*side+i]; }
      else if (conv==3)
        {begin[i*side+j] = img->image[j*side+(side-1-i)]; }
      else if (conv==4)
        {begin[i*side+j] = img->image[(side-1-j)*side+(side-1-i)]; }
      else if (conv==5)
        {begin[i*side+j] = img->image[i*side+(side-1-j)]; }
      else if (conv==6)
        {begin[i*side+j] = img->image[(side-1-i)*side+j]; }
      else if (conv==7)
        {begin[i*side+j] = img->image[(side-1-i)*side+(side-1-j)]; }

    }
  }
  }

  return z;
}


  }
}

bool
iotbx::detectors::Mar345Adaptor::header_read_specific
  (img_handle i, FILE* f, int* g){
  int status = img_read_mar345header (i, f, g);
  std::string detector = get_field("DETECTOR");
  for (std::string::iterator i = detector.begin(); i!=detector.end(); ++i){
    *i=std::tolower(*i);
  }
  if (detector.find("mar")==std::string::npos ||
      detector.find("345")==std::string::npos){
    throw Error("detector type other than mar345 from mar345 reader");
  }
  return (status==0);
}
bool
iotbx::detectors::Mar345Adaptor::data_read_specific
  (img_handle i, FILE* f, int* g){
  int status = img_read_mar345data (i, f, g);
  return (status==0);
}
