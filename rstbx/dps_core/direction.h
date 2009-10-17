#ifndef RSTBX_DIRECTION_H
#define RSTBX_DIRECTION_H

#include <boost/shared_ptr.hpp>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/shared.h>

namespace af = scitbx::af;
namespace rstbx {

typedef scitbx::vec3<double>                 point;
typedef scitbx::mat3<double>                 matrix;

struct Directional_FFT; //forward declaration permitting extract
typedef boost::shared_ptr<Directional_FFT>   fftptr;

struct Direction {
  typedef scitbx::vec3<double>                 point;
  typedef std::size_t                          sztype;
  double psi,phi;
  point dvec;
  sztype kmax;
  double kval;
  double kval0;
  double kval2;
  double kval3;
  af::shared<double> ff;
  sztype m;
  double delta_p;
  double pmin;
  double uc_length;
  Direction(const double &, const double &);
  Direction(const point &);
  Direction();
  void initialize();
  bool is_nearly_collinear(const Direction &)const;
  inline af::shared<double> getff() { return ff; }
  inline point bvec(){return point(uc_length*dvec[0],uc_length*dvec[1],uc_length*dvec[2]);}
  void extract_directional_properties(fftptr, const bool PS = true);
  // fill the original angle with kmax,kval,power spectrum,m,delta_p,pmin,
  // and unit cell length
};

struct kvalcmp {
  bool operator()(const Direction&, const Direction&);
};

} //namespace

#endif //RSTBX_DIRECTION_H
