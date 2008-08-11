#ifndef CCTBX_SGTBX_UTILS_H
#define CCTBX_SGTBX_UTILS_H

namespace cctbx { namespace sgtbx { namespace utils {

  int change_denominator(
    const int *old_num, int old_den,
          int *new_num, int new_den, int n);

  class cmp_i_vec
  {
    public:
      cmp_i_vec(std::size_t n) : n_(n) {}

      bool operator()(const int *a, const int *b) const;

    private:
      std::size_t n_;
  };

}}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_UTILS_H
