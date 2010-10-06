#ifndef CBFLIB_CBF_AD_H
#define CBFLIB_CBF_AD_H
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <exception>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/versa_matrix.h>
#include <include/cbf.h>
#include <include/cbf_simple.h>
#include <cbflib_adaptbx/detectors/basic.h>

#undef cbf_failnez
#define cbf_failnez(x) { int err; err = (x); if (err) { \
  std::cout<<"error code "<<err<<std::endl; throw iotbx::detectors::Error ( \
  "CBFlib error in " #x " "); }}
//#include <boost/timer.hpp>

namespace iotbx {
  namespace detectors {

struct transform_flags {
  /* Operations necessary to put the data in the order
     ELEMENT_X increasing precedence=1
     ELEMENT_Y increasing precedence=2

     Order of operations:
       reverse_slow and reverse_fast (if true) applied first to assure
         each direction==increasing.
       then transpose (if true) to assure ELEMENT_X precedence==1
  */
  bool transpose, reverse_slow, reverse_fast;
};

class CBFAdaptor {
 protected:
  std::string filename;
  FILE *private_file;
  bool read_header_already;
  const char *array_id;
  int dimension [2];
  int id, index;
  std::size_t nelem_read;
  std::size_t i_size1,i_size2;
  int i_rows,i_columns;
  double d_overload, d_wavelength, d_detector_distance, d_pixel_size;
  double d_osc_start, d_osc_range;
  bool data_transposed;
 public:
  cbf_handle cbf_h;
  double beam_index1,beam_index2,beam_center1, beam_center2;

 public:
  inline CBFAdaptor(const std::string& filename):
    filename(filename),read_header_already(false),id(0),data_transposed(false){
    /* Create the cbf */
    cbf_failnez (cbf_make_handle (&cbf_h))
  }

  inline ~CBFAdaptor(){
    /* The CBF manual promises that the CBF library will close the file;
       CBFlib 0.7.8.1 assures this */
    /* Free the cbf */
    cbf_failnez (cbf_free_handle (cbf_h))
  }

  inline double overload(){ read_header();
    cbf_failnez ( cbf_get_overload(cbf_h,0,&d_overload) );
    return d_overload;}

  inline double wavelength(){ read_header();
    cbf_failnez ( cbf_get_wavelength(cbf_h,&d_wavelength) );
    return d_wavelength;}

  inline double distance(){ read_header();
    try {
      cbf_detector detector1;
      cbf_failnez ( cbf_construct_detector(cbf_h,&detector1,0) );
      cbf_failnez ( cbf_get_detector_distance(detector1,&d_detector_distance) );
      cbf_failnez ( cbf_free_detector(detector1) );
    } catch (...) {
      cbf_failnez (cbf_rewind_datablock (cbf_h))
      cbf_failnez (cbf_find_category    (cbf_h, "diffrn_measurement"))
      cbf_failnez (cbf_find_column      (cbf_h, "sample_detector_distance"))
      cbf_failnez (cbf_get_doublevalue        (cbf_h, &d_detector_distance))
      d_detector_distance *= 1000.; //conversion from m to mm
    }
    return d_detector_distance;
  }

  inline double pixel_size(){ read_header();
    try {
      cbf_failnez ( cbf_get_pixel_size(cbf_h,0,1,&d_pixel_size) );
      //to get a handle on negative beam_centers, need to read the pixel size
      //  of both axes to see any that are negative.
    } catch (...) {
      try {
      cbf_detector detector1;
      cbf_failnez ( cbf_construct_detector(cbf_h,&detector1,0) );
      cbf_failnez ( cbf_get_inferred_pixel_size(detector1,1,&d_pixel_size) );
      cbf_failnez ( cbf_free_detector(detector1) );
      } catch (...) {

      }
    }
    return d_pixel_size;
  }

  inline double osc_range(){ read_header();
    cbf_goniometer goniometer1;
    cbf_failnez ( cbf_construct_goniometer(cbf_h,&goniometer1) );
    cbf_failnez (
      cbf_get_rotation_range(goniometer1,0,&d_osc_start,&d_osc_range) );
    cbf_failnez ( cbf_free_goniometer(goniometer1) );
    return d_osc_range;
  }

  inline double osc_start(){ read_header();
    cbf_goniometer goniometer1;
    cbf_failnez ( cbf_construct_goniometer(cbf_h,&goniometer1) );
    cbf_failnez (
      cbf_get_rotation_range(goniometer1,0,&d_osc_start,&d_osc_range) );
    cbf_failnez ( cbf_free_goniometer(goniometer1) );
    return d_osc_start;
  }

  inline void read_header(){
    /* assumptions
       only one detector
       square pixels
       for now, all metadata provided by separate functions
       therefore keep file open until the CBFAdaptor goes out of scope
    */
    if (read_header_already) {return;}
    if (!cbf_h) throw Error("bad CBF handle");
    private_file = std::fopen(filename.c_str(),"rb");
    if (!private_file) throw Error("cbf file BAD_OPEN");
    cbf_failnez (cbf_read_file (cbf_h, private_file, MSG_DIGEST))
    //file must be left open & is closed by the cbf library.

    try {
      //WARNING--To avoid memory leak the cbf_detector
      //  class will have to be wrapped in its own c++ class with defined
      //  destructor.  Same comment for cbf_handle, cbf_goniometer
      cbf_detector detector1;
      cbf_failnez ( cbf_construct_detector(cbf_h,&detector1,0) );
      cbf_failnez (
      cbf_get_beam_center(detector1,&beam_index1,&beam_index2,&beam_center1,&beam_center2) );
      //SCITBX_EXAMINE(beam_index1);
      //SCITBX_EXAMINE(beam_index2);
      //SCITBX_EXAMINE(beam_center1);
      //SCITBX_EXAMINE(beam_center2);
      //Center can be negative; index is positive
      cbf_failnez ( cbf_free_detector(detector1) );
    } catch (iotbx::detectors::Error& e) {throw e;}

    read_header_already = true;
  }

  inline scitbx::af::flex_int read_data(){

    /* Get the image identifier */
    cbf_failnez (cbf_rewind_datablock (cbf_h))
    cbf_failnez (cbf_find_category    (cbf_h, "diffrn_data_frame"))
    cbf_failnez (cbf_find_column      (cbf_h, "array_id"))
    cbf_failnez (cbf_get_value        (cbf_h, &array_id))

    /* Get the image dimensions (second dimension = fast, first = slow) */
    cbf_failnez (cbf_find_category    (cbf_h, "array_structure_list"))
    cbf_failnez (cbf_rewind_row       (cbf_h))
    cbf_failnez (cbf_find_column      (cbf_h, "array_id"))

    dimension [0] = dimension [1] = 0;

    while (cbf_find_nextrow (cbf_h, array_id) == 0){
      cbf_failnez (cbf_find_column      (cbf_h, "precedence"))
      cbf_failnez (cbf_get_integervalue (cbf_h, &index))
      if (index >= 1 && index <= 2){
        cbf_failnez (cbf_find_column (cbf_h, "dimension"))
        cbf_failnez (cbf_get_integervalue (cbf_h, &dimension [2 - index]))
      } else {
        throw Error ("thrown out of cbf_find_nextrow loop");
      }
      cbf_failnez (cbf_find_column (cbf_h, "array_id"))
    }

    if (dimension [0] == 0 || dimension [1] == 0) {
      throw Error ("thrown out of dimension test");
    }

    /* Find the binary data */
    cbf_failnez (cbf_find_category (cbf_h, "array_data"))
    cbf_failnez (cbf_find_column   (cbf_h, "array_id"))
    cbf_failnez (cbf_find_row      (cbf_h, array_id))
    cbf_failnez (cbf_find_column   (cbf_h, "data"))

    //C++ weirdness
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(dimension[0],dimension[1])),
      scitbx::af::init_functor_null<int>() ); //don't waste 0.04 seconds initializing to zero.
    int* begin = z.begin();
    std::size_t sz = z.size();

    /* Read the binary data */
    data_transposed=false;
    //boost::timer T = boost::timer();
    cbf_failnez (cbf_get_integerarray (cbf_h, //cbf handle
                 &id,                         //ptr to binary section identifier
                 begin,                       //array ptr
                 sizeof (int),                //element size
                 1,                           //flag of signed data type
                 sz,                          //elements requested
                 &nelem_read                  //elements read
                 ))
    //std::cout<<"time elapsed for get_integerarray: " <<T.elapsed()<<std::endl;
    SCITBX_ASSERT(sz==nelem_read);

    if (file_is_transposed()) {
      //std::cout<<"transposing data...";
      scitbx::af::matrix_transpose_in_place(z);
      //std::cout<<"done"<<std::endl;
      data_transposed=true;
    }

    return z;
  }

  inline int size1() { read_header();
    cbf_failnez ( cbf_get_image_size(cbf_h,0,0,&i_size1,&i_size2) );
    if (file_is_transposed()) {return i_size2;}
    else {return i_size1;} }

  inline int size2() { read_header();
    cbf_failnez ( cbf_get_image_size(cbf_h,0,0,&i_size1,&i_size2) );
    if (file_is_transposed()) {return i_size1;}
    else {return i_size2;} }

  inline double twotheta(){ // 2-theta will be supported in the future; not yet
    return 0.0;}

  //! True when ELEMENT_X is slow.
  bool file_is_transposed() const;

  std::string raster_description();

  iotbx::detectors::transform_flags
  transform_flags() const;

};

class Goniometer {
  private:
    double d_osc_start, d_osc_range;
    cbf_goniometer goniometer1;
  public:
    Goniometer(const cbf_handle& cbf_h){
      cbf_construct_goniometer(cbf_h,&goniometer1);
      cbf_get_rotation_range(goniometer1,0,&d_osc_start,&d_osc_range);
    }
    double osc_start(){ return d_osc_start; }
    double osc_range(){ return d_osc_range; }
    ~Goniometer(){
      cbf_free_goniometer(goniometer1);
    }
};

  }//namespace detectors
}//namespace iotbx

#endif
