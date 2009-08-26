#include <scitbx/array_family/flex_types.h>
#include <iotbx/detectors/scanbox.h>

namespace Distl {

  //! Class for breaking a detector image into component modules
  /*! For example, the ADSC Quantum 315 detector consists of a 3 x 3 array
      of fiber-optic taper modules.  For the unbinned image, each square
      module is surrounded by a 4-pixel strip of null pixels (value==0).
  */
  class image_divider {
    int nullvalue;
    scitbx::af::flex_int data;
    interval_list slow_tiles,fast_tiles;

    public:
      //! Constructor
      /*!
          data:  an integer*4 array containing the original image
          nullvalue: the integer value assigned to an inactive pixel
          (0 for ADSC; -1 for Pilatus, e.g.)

          The constructor analyzes all the data and determines slow & fast
          strips of values that have a null value.  This allows rectangluar
          module boundaries to be determined on-the-fly, without hard-coding
          a pre-defined configuration from the manufacturer
       */
      image_divider(scitbx::af::flex_int, const int&);
      //! module_count: the total number of modules in this image
      int module_count()const;
      //! module data
      /*!
          This function returns a new rectangular array consisting of
          the raw data from a single module.
          Modules are numbered 0,1,2,... in the same slow+fast sense
          as the data laid out in the raw data array.
       */
      scitbx::af::flex_int tile_data(const int&)const;
      //! inclusive pixel boundaries of the slow module interval
      Distl::interval tile_slow_interval(const int&)const;
      //! inclusive pixel boundaries of the fast module interval
      Distl::interval tile_fast_interval(const int&)const;
  };

}//namespace Distl
