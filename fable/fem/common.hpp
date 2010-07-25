#ifndef FEM_COMMON_HPP
#define FEM_COMMON_HPP

#include <fem/close_chain.hpp>
#include <fem/file_positioning_chain.hpp>
#include <fem/inquire_chain.hpp>
#include <fem/open_chain.hpp>
#include <fem/utils/misc.hpp>

namespace fem {

  struct common
  {
    fem::io io;
  };

  struct cmn_sve
  {
    private:
      // simplified boost::any, with pointer holder instead of value holder
      struct placeholder
      {
        virtual ~placeholder() {}
      };

      template <typename T>
      struct holder : placeholder
      {
        T* ptr;

        holder(T* ptr_) : ptr(ptr_) {}

        ~holder() { delete ptr; }

        private:
          holder& operator=(holder const&);
          holder(holder const&);
      };

      placeholder* content;

      cmn_sve(cmn_sve const&);
      cmn_sve& operator=(cmn_sve const&);

      public:

    cmn_sve() : content(0) {}

    ~cmn_sve() { delete content; }

    bool
    is_called_first_time() { return (content == 0); }

    template <typename T>
    void
    construct()
    {
      content = new holder<T>(new T);
    }

    template <typename T, typename D>
    void
    construct(
      D const& dynamic_parameters)
    {
      content = new holder<T>(new T(dynamic_parameters));
    }

    template <typename T>
    T&
    get()
    {
      return *(static_cast<holder<T>*>(content)->ptr);
    }
  };

  template <typename T>
  void
  no_operation_to_avoid_unused_variable_warning(
    const T&)
  {}

} // namespace fem

#define FEM_CMN_SVE(FUNC) \
  bool is_called_first_time = cmn.FUNC##_sve.is_called_first_time(); \
  if (is_called_first_time) { \
    cmn.FUNC##_sve.construct<FUNC##_save>(); \
  } \
  FUNC##_save& sve = cmn.FUNC##_sve.get<FUNC##_save>(); \
  fem::no_operation_to_avoid_unused_variable_warning(sve)

#define FEM_CMN_SVE_DYNAMIC_PARAMETERS(FUNC) \
  bool is_called_first_time = cmn.FUNC##_sve.is_called_first_time(); \
  if (is_called_first_time) { \
    cmn.FUNC##_sve.construct<FUNC##_save>(cmn.dynamic_parameters); \
  } \
  FUNC##_save& sve = cmn.FUNC##_sve.get<FUNC##_save>()

#endif // GUARD
