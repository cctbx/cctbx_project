#ifndef FEM_DO_HPP
#define FEM_DO_HPP

namespace fem {

  class do_
  {
    int* i_;
    int l_;
    int s_;

    public:

      do_(int* i, int f, int l, int s=1) : i_(i), l_(l), s_(s) { (*i) = f; }

      operator bool()
      {
        if (s_ < 0) return (*i_) >= l_;
        return (*i_) <= l_;
      }

      void
      incr() { (*i_) += s_; }
  };

} // namespace fem

#define FEM_DO(i, f, l) \
  for(i=f; i<=l; i++)

#define FEM_DOSTEP(i, f, l, s) \
  for(fem::do_ fem_do(&(i), f, l, s);fem_do;fem_do.incr())

#endif // GUARD
