#ifndef FEM_MAIN_HPP
#define FEM_MAIN_HPP

#include <fem/intrinsics_extra.hpp>
#include <fem/stop.hpp>
#include <cstdio>

namespace fem {

  inline
  int
  main_with_catch(
    int argc,
    char const* argv[],
    void (*callable)(int argc, char const* argv[]))
  {
    user_plus_system_time();
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
