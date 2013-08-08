#ifndef SCITBX_MATH_GOLAY_H
#define SCITBX_MATH_GOLAY_H

#include <scitbx/array_family/loops.h>
#include <scitbx/array_family/tiny.h>


namespace scitbx { namespace math {

  //! Generator matrix for the Golay (24,12) code.
  /* http://web.usna.navy.mil/~wdj/codes.mpl
     (Note that this is a symmetric matrix.)
   */
  static const int golay_24_12_generator_matrix[] = {
    0,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,0,1,1,1,0,0,0,1,0,
    1,1,0,1,1,1,0,0,0,1,0,1,
    1,0,1,1,1,0,0,0,1,0,1,1,
    1,1,1,1,0,0,0,1,0,1,1,0,
    1,1,1,0,0,0,1,0,1,1,0,1,
    1,1,0,0,0,1,0,1,1,0,1,1,
    1,0,0,0,1,0,1,1,0,1,1,1,
    1,0,0,1,0,1,1,0,1,1,1,0,
    1,0,1,0,1,1,0,1,1,1,0,0,
    1,1,0,1,1,0,1,1,1,0,0,0,
    1,0,1,1,0,1,1,1,0,0,0,1};

  //! Generator for the 4096 binary words of the Golay (24,12) code.
  /*! Joe Fields (fields@math.uic.edu) writes:
      The [24,12,8] extended binary Golay code is a well-known and
      remarkable combinatorial object. It consists of 4096 binary words of
      length 24. Besides the zero word and the word consisting of 24 ones,
      there are 759 words of weight 8. (The weight, or more accurately,
      Hamming weight of a binary word is simply the number of ones in it.)
      There are also 759 words of weight 16, and the remainder of the words
      (2576) are of weight 12. The weight 8 words are called octads and the
      weight 12 words are called dodecads.
      <p>
      http://www2.math.uic.edu/~fields/DecodingGolayHTML/introduction.html
   */
  template <typename BitType=int>
  class golay_24_12_generator
  {
    public:
      //! Initialization of the sequence of words.
      golay_24_12_generator() : loop_(2) {}

      //! True after exactly 4096 calls to next().
      bool
      at_end() const { return loop_.over() != 0; }

      /*! The word returned on the first call is the zero word.
          The last word is returned on the 4096th call and
          consists of 24 ones. An exception is thrown on
          the 4097th call. Use at_end() to avoid the exception.
       */
      af::tiny<BitType, 24>
      next()
      {
        if (at_end()) {
          throw error("golay_24_12_generator is exhausted.");
        }
        af::tiny<int, 12> const& v = loop_();
        af::tiny<BitType, 24> result;
        std::copy(v.begin(), v.end(), result.begin());
        for(int i=0;i<12;i++) {
          int bit = 0;
          for(int j=0;j<12;j++) {
            bit += golay_24_12_generator_matrix[i*12+j] * v[j];
          }
          result[i+12] = static_cast<BitType>(bit % 2);
        }
        loop_.incr();
        return result;
      }

    protected:
      af::nested_loop<af::tiny<int, 12> > loop_;
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_GOLAY_H
