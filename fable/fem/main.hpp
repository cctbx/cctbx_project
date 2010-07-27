#ifndef FEM_MAIN_HPP
#define FEM_MAIN_HPP

#include <fem/intrinsics_extra.hpp>
#include <fem/stop.hpp>
#include <fem/utils/int_types.hpp>
#include <cstdio>

namespace fem {

  inline
  bool
  check_fem_utils_int_types()
  {
    bool result = true;
#define FABLE_LOC(bits, expected_sz) \
    { \
      int sz = static_cast<int>(sizeof(fem::utils::int##bits##_t)); \
      if (sz != expected_sz) { \
        std::fprintf(stderr, \
          "FATAL: sizeof(fem::utils::int%d_t) is %d but should be %d.\n", \
            bits, sz, expected_sz); \
        result = false; \
      } \
    }
    FABLE_LOC( 8, 1)
    FABLE_LOC(16, 2)
    FABLE_LOC(32, 4)
    FABLE_LOC(64, 8)
#undef FABLE_LOC
    if (!result) {
      std::fprintf(stderr,
        "NOTE: fem/utils/int_types.hpp"
        " needs to be adjusted for this platform.\n");
    }
    return result;
  }

  inline
  int
  main_with_catch(
    int argc,
    char const* argv[],
    void (*callable)(int argc, char const* argv[]))
  {
    user_plus_system_time();
    if (!check_fem_utils_int_types()) {
      return 255;
    }
    try {
      callable(argc, argv);
    }
    catch (fem::stop_info const& info) {
      std::fflush(stdout);
      std::fprintf(stderr, "%s\n", info.what());
      std::fflush(stderr);
    }
    catch (std::exception const& e) {
      std::fflush(stdout);
      char const* what = e.what();
      if (what == 0) what = "null";
      std::fprintf(stderr, "std::exception what(): %s\n", what);
      std::fflush(stderr);
      return 1;
    }
    catch (...) {
      std::fflush(stdout);
      std::fprintf(stderr, "Terminated by unknown C++ exception.\n");
      std::fflush(stderr);
      return 2;
    }
    return 0;
  }

  struct dynamic_parameters_from_argv
  {
    int argc;
    char const** argv;
    int i_arg;

    dynamic_parameters_from_argv(
      int argc_,
      char const* argv_[],
      int max_number_of_args)
    :
      argc(argc_),
      argv(argv_),
      i_arg(1)
    {
      if (argc-1 > max_number_of_args) {
        std::ostringstream o;
        o << "Too many command-line arguments (given: "
          << (argc-1)
          << ", max. expected: "
          << max_number_of_args
          << ")";
        throw std::runtime_error(o.str());
      }
    }

    template <typename T>
    dynamic_parameters_from_argv&
    reset_if_given(
      T& value)
    {
      if (i_arg < argc) {
        char const* arg = argv[i_arg];
        std::istringstream i(arg);
        i >> value;
        if (i.fail()) {
          std::ostringstream o;
          o << "Invalid command-line argument (field "
            << i_arg
            << "): \""
            << arg
            << "\"";
          throw std::runtime_error(o.str());
        }
        i_arg++;
      }
      return *this;
    }
  };

} // namespace fem

#endif // GUARD
