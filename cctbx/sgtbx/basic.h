#ifndef CCTBX_SGTBX_BASIC_H
#define CCTBX_SGTBX_BASIC_H

#include <scitbx/mat3.h>
#include <scitbx/rational.h>
#include <cctbx/error.h>

namespace cctbx {

  //! Shorthand for default vec3 type in space group toolbox.
  typedef scitbx::vec3<int> sg_vec3;
  //! Shorthand for default mat3 type in space group toolbox.
  typedef scitbx::mat3<int> sg_mat3;

  //! Space Group Toolbox namespace.
  namespace sgtbx {

  //! Special class for rational vector exceptions.
  /*! Used by sgtbx::tr_vec, sgtbx::rot_mx
   */
  class error_rational_vector : public error
  {
    public:
      //! Constructor.
      error_rational_vector(const char* file, long line,
                            std::string const& msg = "") throw()
      : error(file, line, msg) {}

      //! Virtual destructor.
      virtual ~error_rational_vector() throw() {}
  };

  /*! \brief Maximum number of representative rotation matrices for
      3-dimensional crystallographic space groups.
   */
  static const std::size_t n_max_repr_rot_mx = 24;

  //! Default space_group translation vector denominator.
  const int sg_t_den = 12;
  //! Default change_of_basis rotation matrix denominator.
  const int cb_r_den = 12;
  //! Default change_of_basis translation vector denominator.
  const int cb_t_den = 144;

  //! Checks default rotation and translation denominators for consistency.
  inline void sanity_check()
  {
    CCTBX_ASSERT(sg_t_den % 12 == 0);
    CCTBX_ASSERT(cb_t_den >= 2 * sg_t_den);
    CCTBX_ASSERT(cb_t_den % sg_t_den == 0);
  }

  typedef boost::rational<int> rat;
  typedef scitbx::vec3<rat> vec3_rat;

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_BASIC_H
