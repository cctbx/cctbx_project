#include <string>
#include <mmtbx/alignment/tables.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace mmtbx {
//! Sequence alignment tools.
namespace alignment {

    namespace detail {

    template<typename T>
    inline T maximum( T a, T b ) { return a < b ? b : a ; }

    template<typename T>
    inline T maximum( T a, T b, T c ) { return maximum( maximum(a,b), c ) ; }

}

namespace af=scitbx::af;


inline
int identity_f(char a, char b) {
  return int(a == b);
}

inline
int identity_f(std::string a, std::string b) {
  return int(a == b);
}

inline
int dayhoff_f(char a, char b) {
  std::string AAs = "ACDEFGHIKLMNPQRSTVWY";
  int i = AAs.find_first_of(a);
  int j = AAs.find_first_of(b);
  if (i<0 || j<0) return 0;
  return dayhoff_mdm78_similarity_scores(i,j);
}

inline
int blosum50_f(char a, char b) {
  std::string AAs = "ABCDEFGHIKLMNPQRSTVWXYZ";
  int i = AAs.find_first_of(a);
  int j = AAs.find_first_of(b);
  if (i<0 || j<0) return 0;
  return blosum50_similarity_scores(i,j);
}

inline
int dayhoff_f(std::string a, std::string b) {
  return 0;
}

inline
int blosum50_f(std::string a, std::string b) {
  return 0;
}

template<typename T>
int sim_fun_select(T a, T b, std::string descr) {
  // int (*sim_fun)(T, T);
  if (descr == "identity") return identity_f(a,b);
  else if (descr == "dayhoff") return dayhoff_f(a,b);
  else if (descr == "blosum50") return blosum50_f(a,b);
  return 0;
}


class align {

private:
  float gop;
  float gep;

  inline
  float gap_cost(float width) {return gop + width*gep;}

public:

  af::versa<float, af::c_grid<2> > M;
  af::versa<float, af::c_grid<2> > D;
  af::versa<float, af::c_grid<2> > I;
  af::versa<float, af::c_grid<2> > E;

  template <typename DataType>
  // we going to need std::string for sequence alignment
  // and array of strings for atom names alignment af::const_ref<std::string> const&
  align(
      DataType const& seq_a,
      DataType const& seq_b,
      af::const_ref<float> const& masking,
      std::string const& style = "global",
      float gap_opening_penalty=1,
      float gap_extension_penalty=1,
      std::string const& similarity_function="identity")
  {
    using namespace std;
    gop = gap_opening_penalty;
    gep = gap_extension_penalty;
    int m = seq_a.size();
    int n = seq_b.size();

    M.resize(af::c_grid<2>(m+1, n+1), 0);
    D.resize(af::c_grid<2>(m+1, n+1), 0);
    I.resize(af::c_grid<2>(m+1, n+1), 0);
    E.resize(af::c_grid<2>(m+1, n+1), 0);
    af::shared<float> G;
    G.resize(m+1); // scale on gap penalty at position i for sequence a.
                   // If all 1 it is standard gap penalty

    for (int i=0; i<m+1; i++){
      if (i==0 || i==m) {
        G[i]=gap_cost(1);
      } else {
       float scale = masking[i];
       G[i] = scale*gap_cost(1); // cost of gap at i i_scale bigger than default
      }
    }

    if (style == "global")
      for (int i=0; i<m+1; i++) M(i,0) = -gap_cost(i);
      for (int i=0; i<n+1; i++) M(0,i) = -gap_cost(i);

    for (int i=1; i<=m; i++)
      for (int j=1; j<=n; j++) {

        if (i==1) D(i,j) = M(i-1,j)-gap_cost(1);
        else D(i,j) = std::max( M(i-1,j)-G[i], D(i-1,j)-gep);

        if (j==1) I(i,j) = M(i,j-1)-gap_cost(1);
        else I(i,j) = std::max(M(i,j-1)-G[i], I(i,j-1)-gep);

        M(i,j) = detail::maximum( M(i-1,j-1)+sim_fun_select(seq_a[i-1],seq_b[j-1],similarity_function), D(i,j), I(i,j));
        if (style=="local") M(i,j) = std::max(M(i,j), float(0));

        if ( M(i,j)==D(i,j) ) E(i,j) = 1;    // deletion, i.e. of A[i]
        else if (M(i,j)==I(i,j)) E(i,j) = -1; // insertion, i.e. of B[j]
        else E(i,j) = 0;
      }
  }

}; // class align

}} // namespace mmtbx::alignment
