#ifndef DISTL_WRAPPER_H
#define DISTL_WRAPPER_H

#include <string>
#include <cmath>

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/constants.h>

#include <spotfinder/core_toolbox/libdistl.h>

namespace af = scitbx::af;
namespace spotfinder {
namespace distltbx {

typedef af::tiny<int, 2> intxy;

typedef af::shared<Distl::spot> spot_list_t;
typedef Distl::spot w_spot;

class SpotError : public std::exception {
private:
  std::string s;
  const char* file;
public:
  inline SpotError(std::string s):s(s){}
  virtual const char* what() const throw();
  virtual ~SpotError() throw();
};
inline const char* SpotError::what() const throw() {
  char m[120];
  sprintf (m,"%s\n",s.c_str());
  std::string mess(m);
  return mess.c_str();
}
inline SpotError::~SpotError() throw() {}

inline af::flex_int linear_char_scratchpad(int rows,int cols){
  // int array of zeroes; user-defined purpose
  af::flex_int z(af::flex_grid<>(cols,rows));
  return z;
}

class w_Distl {
public:
  w_Distl(std::string,bool);
  Distl::diffimage finder;
  void set_resolution_outer(const double&);
  void setspotimg(const double&, const double&, const double&,
                  const double&, const double&, af::flex_int const&,
                  const int&,
                  const double&);
  void set_tiling(const std::string&);
  af::flex_double Z_data();//DISTL's pixelwise Z-score array
  af::flex_int mod_data();//user-defined modified dataset

  inline void get_underload(){
    finder.image_geometry = Distl::get_image_geometry(finder.pixelvalue);
    finder.underloadvalue = finder.get_underload();
  }

  inline void set_minimum_spot_area(const int& A){
    finder.spotarealowcut = A;
    finder.spotbasesize = A;}

  inline void set_minimum_signal_height(const double& A){
    finder.bgupperint[0] = A;
    finder.bgupperint[1] = A;
    finder.bgupperint[2] = A;}

  inline void set_minimum_spot_height(const double& A){
    finder.difflowerint = A;}

  inline void set_spot_area_maximum_factor(const double& A){
    finder.spotareamaxfactor = A;}

  inline void set_scanbox_windows(af::shared<int> A){
    finder.scanboxsize[0] = A[0];
    finder.scanboxsize[1] = A[1];
    finder.scanboxsize[2] = A[2];
    }

  inline void parameter_guarantees(){
    // difflowerint (minimum_spot_height) >= bgupperint (minimum_signal_height)
    double* max_bg = std::max_element(finder.bgupperint, finder.bgupperint+2);
    if (finder.difflowerint < *max_bg) { finder.difflowerint = *max_bg; }

    if (finder.spotareamaxfactor < 1.0) { finder.spotareamaxfactor = 1.0; }
  }

  inline void pxlclassify(){ finder.pxlclassify(); }

  inline void search_icerings(){ finder.search_icerings(); }

  inline void search_maximas(){ finder.search_maximas(); }

  inline void search_spots(){ finder.search_spots(); }

  inline void search_overloadpatches(){
    finder.search_overloadpatches();
    finder.imgresolution();
  }

  void finish_analysis();
  scitbx::af::shared<w_spot> spots;
  inline int nicerings() { return finder.icerings.size(); }

  bool isIsolated(const w_spot&, const double&) const;
  scitbx::af::shared<Distl::icering> icerings;
  inline double imgresol() { return finder.imgresol; }

  inline
  af::shared<double>
  background_resolutions(){
    return finder.scanbox_background_resolutions; }

  inline
  af::shared<double>
  background_means() const {
    return finder.scanbox_background_means; }

  inline
  af::shared<double>
  background_wndw_sz() const {
    return finder.scanbox_background_wndw_sz; }
};

} //namespace

} //namespace

#endif //DISTL_WRAPPER_H
