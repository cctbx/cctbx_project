
#include <simtbx/nanoBragg/nanoBragg.h>

namespace simtbx {
namespace nanoBragg {

class derivative_manager{
  public:
    derivative_manager();
    void initialize(int sdim, int fdim);

    double* floatimage;
    af::flex_double raw_pixels;
    double value;
};

class rot_manager: public derivative_manager{
  public:
    rot_manager();
}; // end of rot_manager

class diffBragg: public nanoBragg{
  public:
  diffBragg(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
            int verbose, int panel_id);

  void initialize_managers();
  void add_diffBragg_spots();
  mat3 RX, RY, RZ, RXYZ;

  /* derivative managers */
  rot_manager rotX_man;
  rot_manager rotY_man;
  rot_manager rotZ_man;

}; // end of diffBragg

} // end of namespace nanoBragg
} // end of namespace simtbx
