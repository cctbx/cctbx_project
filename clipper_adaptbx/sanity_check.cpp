#include <clipper/core/clipper_util.h>
#include <iostream>
#include <stdio.h>

namespace clipper { namespace {

template <typename FloatType>
void
show_as_unsigned(const char* label, FloatType const& x)
{
  std::cout << label << std::endl;
  std::cout << "  value: " << x << std::endl;
  if (sizeof(FloatType) == sizeof(unsigned int)) {
    printf("  int(%u):  0x%x",
      static_cast<unsigned>(sizeof(unsigned int)),
      *((unsigned int *)(&x)));
  }
  else if (sizeof(FloatType) == sizeof(unsigned long)) {
    printf("  int(%u):  0x%lxL",
      static_cast<unsigned>(sizeof(unsigned long)),
      *((unsigned long *)(&x)));
  }
#if !defined(__osf__)
  else if (sizeof(FloatType) == sizeof(unsigned long long)) {
#if defined(_MSC_VER)
    printf("  int(%u):  0x%I64xLL",
      sizeof(unsigned long long),
#else
    printf("  int(%u):  0x%llxLL",
      static_cast<unsigned>(sizeof(unsigned long long)),
#endif
      *((unsigned long long *)(&x)));
  }
#endif
  else {
    printf("  no matching integer type");
  }
  printf("\n");
  std::cout << "  is_nan:  " << Util::is_nan(x) << std::endl;;
  std::cout << "  isnan:   " << Util::isnan(x) << std::endl;;
  std::cout << "  is_null: " << Util::is_null(x) << std::endl;;
  std::cout << "  x!=x:    " << (x!=x) << std::endl;;
#if defined(FP_NAN) && !defined(__APPLE__)
  std::cout << "  isnanf() system macro: " << isnanf(x) << std::endl;
  std::cout << "  isnan()  system macro: " << isnan(x) << std::endl;
#endif
}

void
sanity_check()
{
  std::cout << "sizeof(itype32):  " << sizeof(itype32) << std::endl;
  std::cout << "sizeof(uitype32): " << sizeof(uitype32) << std::endl;
  std::cout << "sizeof(ftype32):  " << sizeof(ftype32) << std::endl;
  std::cout << "sizeof(itype64):  " << sizeof(itype64) << std::endl;
  std::cout << "sizeof(uitype64): " << sizeof(uitype64) << std::endl;
  std::cout << "sizeof(ftype64):  " << sizeof(ftype64) << std::endl;
  show_as_unsigned("nan float", Util::nanf());
  show_as_unsigned("nan double", Util::nand());
  {
    float x;
    Util::set_null(x);
    show_as_unsigned("null float", x);
  }
  {
    double x;
    Util::set_null(x);
    show_as_unsigned("null double", x);
  }
#if !defined(__osf__) || defined(_IEEE_FP)
  show_as_unsigned("sqrt(-1) double", sqrt(-1.));
#endif
}

}} // namespace clipper::<anonymous>

int
main()
{
  clipper::sanity_check();
  return 0;
}
