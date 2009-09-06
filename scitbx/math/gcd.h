#ifndef SCITBX_MATH_GCD_H
#define SCITBX_MATH_GCD_H

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
# define SCITBX_MATH_GCD_USING_ASM
#endif

namespace scitbx { namespace math {

  inline
  int
  gcd_int_simple(
    int a,
    int b)
  {
    for(;;) {
      if (b == 0) return (a < 0 ? -a : a);
      int edx = a % b;
      a = b;
      b = edx;
    }
  }

#if defined(SCITBX_MATH_GCD_USING_ASM)
  /* gcc 3.2 32-bit: ca. 30% faster than gcd_int_simple()
     gcc 4.0 32-bit: ca. 10% faster
     gcc 4.1 64-bit: ca. 15% faster
     gcc 4.4, icc 9.1 64-bit: practically same speed
  */
  inline
  int
  gcd_int_asm(
    int a,
    int b)
  {
    int result;
    __asm__(
      "  movl %1, %%eax\n"
      "  movl %2, %%ecx\n"
      ".L_for_%=:\n"
      "  cmpl $0, %%ecx\n" // if b == 0
      "  je .L_return_a_or_minus_a_%=\n"
      "  movl %%eax, %%edx\n"
      "  sarl $31, %%edx\n"
      "  idivl %%ecx\n" // edx = a % b
      "  movl %%ecx, %%eax\n" // a = b
      "  movl %%edx, %%ecx\n" // b = edx
      "  jmp .L_for_%=\n"
      ".L_return_a_or_minus_a_%=:\n"
         // next five lines: slightly faster than compare & jump
      "  movl %%eax, %%ecx\n"
      "  cltd\n"
      "  xorl %%edx, %%ecx\n"
      "  subl %%edx, %%ecx\n"
      "  movl %%ecx, %%eax\n"
        : "=a"(result) /* output */
        : "r"(a), "r"(b) /* input */
        : "%ecx", "%edx"); /* clobbered registers */
    return result;
  }
#endif

  inline
  int
  gcd_int(
    int a,
    int b)
  {
#if defined(SCITBX_MATH_GCD_USING_ASM)
    return gcd_int_asm(a, b);
#else
    return gcd_int_simple(a, b);
#endif
  }

}} // namespace scitbx::math

#endif // GUARD
