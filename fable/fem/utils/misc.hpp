#ifndef FEM_UTILS_MISC_HPP
#define FEM_UTILS_MISC_HPP

namespace fem { namespace utils {

  struct one_time_flag
  {
    private:
      bool state;
      public:

    one_time_flag() : state(true) {}

    bool
    operator()()
    {
      if (!state) return false;
      state = false;
      return true;
    }
  };

  template <typename T>
  struct hide
  {
    T hidden;
  };

}} // namespace fem::utils

#endif // GUARD
