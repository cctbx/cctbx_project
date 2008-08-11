#ifndef SCITBX_LBFGSB_FORTRAN_INTERFACE_H
#define SCITBX_LBFGSB_FORTRAN_INTERFACE_H

#include <scitbx/lbfgsb/raw.h>

extern "C" {

  // This signature is correct for Tru64 Unix with f77 and IA32 Linux with g77.
  // Other platforms may require different signatures.
  void
  setulb_(
    const int* n,
    const int* m,
    double* x,
    double* l,
    double* u,
    int* nbd,
    double* f,
    double* g,
    const double* factr,
    const double* pgtol,
    double* wa,
    const int* iwa,
    char* task,
    const int* iprint,
    char* csave,
    int* lsave,
    int* isave,
    double* dsave,
    int task_len,
    int csave_len);
}

namespace scitbx { namespace lbfgsb { namespace fortran_interface {

  void
  to_character_buf(
    std::string const& s,
    char* b,
    int b_len)
  {
    std::size_t i = 0;
    for(;i<s.size();i++) b[i] = s[i];
    for(;i<b_len;i++) b[i] = ' ';
  }

  std::string
  from_character_buf(
    const char* b,
    int b_len)
  {
    std::string result;
    int l = b_len-1;
    for(;l>=0;l--) {
      if (b[l] != ' ') break;
    }
    return std::string(b, l+1);
  }

  void
  to_logical_buf(
    const bool* b,
    int* l,
    int l_len)
  {
    for(std::size_t i=0;i<l_len;i++) l[i] = (b[i] ? 1 : 0);
  }

  void
  from_logical_buf(
    const int* l,
    int l_len,
    bool* b)
  {
    for(std::size_t i=0;i<l_len;i++) b[i] = (l[i] ? true : false);
  }

  using namespace raw;

  void
  setulb(
    int const& n,
    int const& m,
    ref1<double> const& x,
    ref1<double> const& l,
    ref1<double> const& u,
    ref1<int> const& nbd,
    double& f,
    ref1<double> const& g,
    double const& factr,
    double const& pgtol,
    ref1<double> const& wa,
    ref1<int> const& iwa,
    std::string& task,
    int const& iprint,
    std::string& csave,
    ref1<bool> const& lsave,
    ref1<int> const& isave,
    ref1<double> const& dsave)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(x.size() == n);
    SCITBX_ASSERT(l.size() == n);
    SCITBX_ASSERT(u.size() == n);
    SCITBX_ASSERT(nbd.size() == n);
    SCITBX_ASSERT(g.size() == n);
    SCITBX_ASSERT(wa.size() == 2*m*n+4*n+12*m*m+12*m);
    SCITBX_ASSERT(iwa.size() == 3*n);
    SCITBX_ASSERT(lsave.size() == 4);
    SCITBX_ASSERT(isave.size() == 44);
    SCITBX_ASSERT(dsave.size() == 29);
#endif
    static double epsmch = math::floating_point_epsilon<double>::get();
    const int task_len = 60;
    const int csave_len = 60;
    const int lsave_len = 4;
    char task_buf[task_len+1];
    char csave_buf[csave_len+1];
    int lsave_buf[lsave_len];
    SCITBX_ASSERT(task.size() <= task_len);
    SCITBX_ASSERT(csave.size() <= csave_len);
    to_character_buf(task, task_buf, task_len);
    to_character_buf(csave, csave_buf, csave_len);
    to_logical_buf(lsave.begin(), lsave_buf, lsave_len);
    ::setulb_(
      &n, &m, x.begin(), l.begin(), u.begin(), nbd.begin(),
      &f, g.begin(), &factr, &pgtol, wa.begin(), iwa.begin(),
      task_buf, &iprint, csave_buf, lsave_buf, isave.begin(),
      dsave.begin(), task_len, csave_len);
    task = from_character_buf(task_buf, task_len);
    csave = from_character_buf(csave_buf, csave_len);
    from_logical_buf(lsave_buf, lsave_len, lsave.begin());
    // function dpmeps may produce incorrect results.
    // overwrite with the correct values.
    dsave(3) = factr*epsmch;
    dsave(5) = epsmch;
  }

}}} // namspace scitbx::lbfgsb::fortran_interface

#endif // SCITBX_LBFGSB_FORTRAN_INTERFACE_H
