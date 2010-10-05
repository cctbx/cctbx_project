#ifndef CBFLIB_MAR_AD_H
#define CBFLIB_MAR_AD_H
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <exception>
#include <scitbx/array_family/flex_types.h>
#include <examples/img.h>
#include <cbflib_adaptbx/detectors/basic.h>

extern "C" {
int img_set_tags (img_handle img, int tags);
int img_read_mar345header (img_handle img, FILE *file, int *org_data);
int img_read_mar345data (img_handle img, FILE *file, int *org_data);
}

namespace af = scitbx::af;

namespace iotbx {
  namespace detectors {

class AnyImgAdaptor {
 private:
  std::string filename;
  img_handle img;
  int img_status;
  bool read_header_ok;
  long data_read_position;
  bool read_data_ok;
  FILE* private_file;
  int img_org_data[4];

 public:
  AnyImgAdaptor(const std::string& filename);
  void read_header();
  inline virtual bool header_read_specific(img_handle, FILE*, int*){return false;}
  //can't make it abstract or Boost python barfs
  void read_data();
  inline virtual bool data_read_specific(img_handle, FILE*, int*){return false;}
  inline ~AnyImgAdaptor() {
    img_free_handle (img);
  }
  std::string get_field(const std::string& field_name);
  inline double get_number(const std::string& field_name){
    return img_get_number(img,field_name.c_str());
  }
  inline int columns(){return img_columns(img);}
  inline int rows(){return img_rows(img);}
  inline int size1() { read_header(); return img_org_data[0]; }
  inline int size2() { read_header();return img_org_data[1]; }
  af::flex_int rawdata();
};

class Mar345Adaptor: public AnyImgAdaptor {
 public:
  inline Mar345Adaptor(const std::string& filename):
    AnyImgAdaptor(filename){}

  bool header_read_specific(img_handle i, FILE* f, int* g);
  bool data_read_specific(img_handle i, FILE* f, int* g);
  inline double pixel_size() { read_header();
    return get_number("PIXEL SIZE");
  }
  inline double wavelength() { read_header();
    return get_number("WAVELENGTH");
  }
  inline double distance() { read_header();
    return get_number("DISTANCE");
  }
  inline double gain() { read_header();
    return 1.55;
  }
  inline double overload() { read_header();
    return 249862; //determined from a single image.
  }
  inline double osc_range() { read_header();
    return get_number("OSCILLATION RANGE");
  }
  inline double osc_start() { read_header();
    return get_number("PHI");
  }
  inline double twotheta() { read_header();// XXX Unknown whether 0.001 factor is needed
    return get_number("TWOTHETA");
  }
  inline std::string test(std::string item) {read_header(); return get_field(item);}
};

  }
}

#endif
