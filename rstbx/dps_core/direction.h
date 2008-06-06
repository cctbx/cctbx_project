#ifndef RSTBX_DIRECTION_H
#define RSTBX_DIRECTION_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/shared.h>

namespace af = scitbx::af;
namespace labelit {
namespace dptbx {

typedef scitbx::vec3<double>                 point;
typedef scitbx::mat3<double>                 matrix;

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
  bool is_nearly_collinear(const Direction &)const;
  inline af::shared<double> getff() { return ff; }
};

struct kvalcmp {
  bool operator()(const Direction&, const Direction&);
};

} //namespace

} //namespace

#endif //RSTBX_DIRECTION_H
