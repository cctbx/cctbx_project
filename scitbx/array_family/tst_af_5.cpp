#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/shared_algebra.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

  template <typename ElementType>
  struct exercise_complex_special
  {
    static void run()
    {
      af::shared<int> i;
      i.assign(af::tiny<int, 3>(0, 1, 2));
      af::shared<ElementType> e;
      e.assign(af::tiny<ElementType, 3>(0.1, 0.2, 0.3));
      af::shared<std::complex<ElementType> > c;
      c.assign(af::tiny<ElementType, 3>(0.4, 0.5, 2.0));
      { af::shared<ElementType> r = af::real(c);
        check_true(__LINE__, r[2] == std::real(c[2])); }
      { af::shared<ElementType> r = af::imag(c);
        check_true(__LINE__, r[2] == std::imag(c[2])); }
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
      { af::shared<ElementType> r = af::abs(c);
        check_true(__LINE__, r[2] == std::abs(c[2])); }
#endif
      { af::shared<ElementType> r = af::arg(c);
        check_true(__LINE__, r[2] == std::arg(c[2])); }
      { af::shared<ElementType> r = af::norm(c);
        check_true(__LINE__, r[2] == std::norm(c[2])); }
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
      { af::shared<std::complex<ElementType> > r = af::pow(c, i);
        check_true(__LINE__, r[2] == std::pow(c[2], i[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, i[0]);
        check_true(__LINE__, r[2] == std::pow(c[2], i[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c[0], i);
        check_true(__LINE__, r[2] == std::pow(c[0], i[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, e);
        check_true(__LINE__, r[2] == std::pow(c[2], e[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, e[0]);
        check_true(__LINE__, r[2] == std::pow(c[2], e[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c[0], e);
        check_true(__LINE__, r[2] == std::pow(c[0], e[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, c);
        check_true(__LINE__, r[2] == std::pow(c[2], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, c[0]);
        check_true(__LINE__, r[2] == std::pow(c[2], c[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c[0], c);
        check_true(__LINE__, r[2] == std::pow(c[0], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(e, c);
        check_true(__LINE__, r[2] == std::pow(e[2], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(e, c[0]);
        check_true(__LINE__, r[2] == std::pow(e[2], c[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(e[0], c);
        check_true(__LINE__, r[2] == std::pow(e[0], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::polar(e, e);
        check_true(__LINE__, approx_equal(r[2].real(),
                                          std::polar(e[2], e[2]).real()));
        check_true(__LINE__, approx_equal(r[2].imag(),
                                          std::polar(e[2], e[2]).imag())); }
      { af::shared<std::complex<ElementType> > r = af::polar(e, e[0]);
        check_true(__LINE__, approx_equal(r[2].real(),
                                          std::polar(e[2], e[0]).real()));
        check_true(__LINE__, approx_equal(r[2].imag(),
                                          std::polar(e[2], e[0]).imag())); }
      { af::shared<std::complex<ElementType> > r = af::polar(e[0], e);
        check_true(__LINE__, approx_equal(r[2].real(),
                                          std::polar(e[0], e[2]).real()));
        check_true(__LINE__, approx_equal(r[2].imag(),
                                          std::polar(e[0], e[2]).imag())); }
#endif
    }
  };

} // namespace <anonymous>

int main(int argc, char* argv[])
{
  for(;;) {
    exercise_complex_special<double>::run();
    if (argc == 1) break;
  }
  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }
  return 0;
}
