#ifndef FEM_UTILS_MISC_HPP
#define FEM_UTILS_MISC_HPP

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

}} // namespace fem::utils

#endif // GUARD
