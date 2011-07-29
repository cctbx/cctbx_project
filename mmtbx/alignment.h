#include <string>

namespace mmtbx {
//! Sequence alignment tools.
namespace alignment {

  namespace detail {

  template <class T> class Matrix {
   public:
    Matrix( int n1, int n2 ) { N1 = n1; N2 = n2; data = new T[N1*N2]; }
    ~Matrix() { delete[] data; }
    T& operator() ( int i1, int i2 ) { return data[i1*N2+i2]; }
   private:
    int N1, N2;
    T* data;
  };

  }

  //! Based on implementation by Kevin Cowtan.
  /*! Date: Tue, 28 Aug 2007 07:49:27 +0100
      From: Kevin Cowtan
      To: "Ralf W. Grosse-Kunstleve"
      Subject: Sequence alignment

      I hacked up a little sequence alignment program on the way home - feel
      free to use this code however you wish.

      Scores are +1 for match, 0 for mismatch, -0.5 for insertion/deletion.
    */
  struct pairwise_global
  {
    std::string result1;
    std::string result2;

    pairwise_global() {}

    pairwise_global(
      std::string seq1,
      std::string seq2 )
    {
      enum dirn { NUL, U, L, UL };  // directions: null, up, left, diag
      const float match = 1.0;
      const float mismatch = 0.0;
      const float insdel = -0.5;

      // pad sequences at start to allow first symbol to be aligned
      std::string s1 = " " + seq1;
      std::string s2 = " " + seq2;
      int n1 = s1.length();
      int n2 = s2.length();

      // initilize matrices.
      detail::Matrix<float> scores(n1,n2);
      detail::Matrix<char> dirns(n1,n2);
      for ( int i1 = 0; i1 < n1; i1++ )
        for ( int i2 = 0; i2 < n2; i2++ ) {
          scores(i1,i2) = -1.0e6;
          dirns(i1,i2) = NUL;
        }

      // now fill first row/col
      for ( int i1 = 1; i1 < n1; i1++ ) {
        scores(i1,0) = scores(i1-1,0) + insdel;
        dirns(i1,0) = U;
      }
      for ( int i2 = 1; i2 < n2; i2++ ) {
        scores(0,i2) = scores(0,i2-1) + insdel;
        dirns(0,i2) = L;
      }

      // fill the rest of the matrix
      for ( int i1 = 1; i1 < n1; i1++ )
        for ( int i2 = 1; i2 < n2; i2++ ) {
          // calc bonus for a match at this position
          int s = ( s1[i1] == s2[i2] ) ? 1 : 0;
          if (s1[i1] == 'X' && s2[i2] == 'X') {
            s = 0;
          }
          // calc best score obtainable for this position
          float sul = scores(i1-1,i2-1) + s;
          float su = scores(i1-1,i2) + s + insdel;
          float sl = scores(i1,i2-1) + s + insdel;
          // and select
          if ( sul >= su && sul >= sl ) {
            scores(i1,i2) = sul;
            dirns(i1,i2) = UL;
          } else if ( su > sl ) {
            scores(i1,i2) = su;
            dirns(i1,i2) = U;
          } else {
            scores(i1,i2) = sl;
            dirns(i1,i2) = L;
          }
        }

      // now trace backwards to build up the best sequence alignment
      std::string r1, r2;
      { // scope to avoid Visual C++ 7.1 warnings
      int i1, i2;
      i1 = n1 - 1;
      i2 = n2 - 1;
      while ( dirns(i1,i2) != NUL ) {
        if ( dirns(i1,i2) == UL ) {
          r1 = s1[i1] + r1;
          r2 = s2[i2] + r2;
          i1 = i1 - 1;
          i2 = i2 - 1;
        } else if ( dirns(i1,i2) == U ) {
          r1 = s1[i1] + r1;
          r2 = '-' + r2;
          i1 = i1 - 1;
        } else {
          r1 = '-' + r1;
          r2 = s2[i2] + r2;
          i2 = i2 - 1;
        }
      }
      }

      result1 = r1;
      result2 = r2;
    }
  };

}} // namespace mmtbx::alignment
