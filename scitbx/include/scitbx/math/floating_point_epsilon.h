#ifndef SCITBX_MATH_FLOATING_POINT_EPSILON_H
#define SCITBX_MATH_FLOATING_POINT_EPSILON_H

#include <scitbx/serialization/base_256.h>

namespace scitbx { namespace math {

  //! Returns value trimmed to the true precision of FloatType.
  /*! Some optimizers maintain a higher precision for some
      calculations. To avoid this and trim value to the true
      precision of FloatType, value is converted to a base-256
      string and then back to FloatType.
   */
  template <typename FloatType>
  FloatType
  trim_cast(FloatType const& value)
  {
    namespace base_256 = scitbx::serialization::base_256;
    char buf[4*sizeof(FloatType)];
    base_256::floating_point::to_string(buf, value);
    base_256::floating_point::from_string<FloatType> fs(buf);
    return fs.value;
  }

  /*! \brief Dynamic determination of the smallest floating point
      number such that 1 + floating_point_epsilon() differs from 1.
   */
  template <typename FloatType>
  struct floating_point_epsilon
  {
    //! Perform the dynamic determination.
    /*! Original comments from the FORTRAN code (Lbfgsb.2.1, function dpmeps):

      This subroutine computes the machine precision parameter
      dpmeps as the smallest floating point number such that
      1 + dpmeps differs from 1.

      This subroutine is based on the subroutine machar described in

      W. J. Cody,
      MACHAR: A subroutine to dynamically determine machine parameters,
      ACM Transactions on Mathematical Software, 14, 1988, pages 303-311.

      The subroutine statement is:

        subroutine dpeps(dpmeps)

      where

        dpmeps is a double precision variable.
          On entry dpmeps need not be specified.
          On exit dpmeps is the machine precision.

      MINPACK-2 Project. February 1991.
      Argonne National Laboratory and University of Minnesota.
      Brett M. Averick.
     */
    static FloatType
    get()
    {
      long i,ibeta,irnd,it,itemp,negep;
      FloatType a,b,beta,betain,betah,temp,tempa,temp1;
      FloatType zero = 0;
      FloatType one = 1;
      FloatType two = 2;
      // determine ibeta, beta ala malcolm.
      a = one;
      b = one;
      while (true) {
        a = trim_cast(a + a);
        temp = trim_cast(a + one);
        temp1 = trim_cast(temp - a);
        if (trim_cast(temp1 - one) != zero) break;
      }
      while (true) {
        b = trim_cast(b + b);
        temp = trim_cast(a + b);
        itemp = static_cast<long>(temp - a);
        if (itemp != 0) break;
      }
      ibeta = itemp;
      beta = static_cast<FloatType>(ibeta);
      // determine it, irnd.
      it = 0;
      b = one;
      while (true) {
        it = it + 1;
        b = trim_cast(b * beta);
        temp = trim_cast(b + one);
        temp1 = trim_cast(temp - b);
        if (trim_cast(temp1 - one) != zero) break;
      }
      irnd = 0;
      betah = trim_cast(beta/two);
      temp = trim_cast(a + betah);
      if (trim_cast(temp - a) != zero) irnd = 1;
      tempa = trim_cast(a + beta);
      temp = trim_cast(tempa + betah);
      if ((irnd == 0) && (trim_cast(temp - tempa) != zero)) irnd = 2;
      // determine dpmeps.
      negep = it + 3;
      betain = trim_cast(one/beta);
      a = one;
      for(i=0;i<negep;i++) {
        a = trim_cast(a*betain);
      }
      while (true) {
        temp = trim_cast(one + a);
        if (trim_cast(temp - one) != zero) break;
        a = a*beta;
      }
      FloatType result = a;
      if ((ibeta == 2) || (irnd == 0)) return result;
      a = trim_cast((a*(one + a))/two);
      temp = trim_cast(one + a);
      if (trim_cast(temp - one) != zero) result = a;
      return result;
    }
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_FLOATING_POINT_EPSILON_H
