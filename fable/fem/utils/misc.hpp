#ifndef FEM_UTILS_MISC_HPP
#define FEM_UTILS_MISC_HPP

#include <fem/size_t.hpp>

namespace fem { namespace utils {

  template <typename T>
  struct hide
  {
    T hidden;
  };

  struct noncopyable // clone of boost::noncopyable
  {
    protected:
      noncopyable() {}

    private:
      noncopyable(noncopyable const&);
      noncopyable const& operator=(noncopyable const&);
  };

  template <typename T, size_t SmallSize=256>
  struct simple_buffer : noncopyable
  {
    T small_space[SmallSize];
    T* space;

    explicit
    simple_buffer(
      size_t size)
    :
      space(size <= SmallSize ? small_space : new T[size])
    {}

    ~simple_buffer()
    {
      if (space != small_space) delete[] space;
    }
  };

}} // namespace fem::utils

#endif // GUARD
