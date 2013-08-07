#ifndef SCITBX_MATH_GCD_H
#define SCITBX_MATH_GCD_H

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__)) \
 && (!defined(__clang__)) \
 && (defined(__APPLE_CC__) && __APPLE_CC__ >= 5465) // OS X 10.5
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
      int next_b = a % b;
      a = b;
      b = next_b;
    }
  }

  inline
  long
  gcd_long_simple(
    long a,
    long b)
  {
    for(;;) {
      if (b == 0) return (a < 0 ? -a : a);
      long next_b = a % b;
      a = b;
      b = next_b;
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
  gcd_int32_asm(
    int a,
    int b)
  {
    __asm__(
      "  movl %2, %%ecx\n"
      ".L_for_%=:\n"
      "  cmpl $0, %%ecx\n" // if b == 0
      "  je .L_return_a_or_minus_a_%=\n"
      "  movl %0, %%edx\n"
      "  sarl $31, %%edx\n"
      "  idivl %%ecx\n" // next_b = a % b
      "  movl %%ecx, %0\n" // a = b
      "  movl %%edx, %%ecx\n" // b = next_b
      "  jmp .L_for_%=\n"
      ".L_return_a_or_minus_a_%=:\n"
         // next five lines: slightly faster than compare & jump
      "  movl %0, %%ecx\n"
      "  cltd\n"
      "  xorl %%edx, %%ecx\n"
      "  subl %%edx, %%ecx\n"
      "  movl %%ecx, %0\n"
        : "=a"(a) /* output */
        : "0"(a), "r"(b) /* input */
        : "cc", "%ecx", "%edx"); /* clobbered registers */
    return a;
  }
#endif

#if defined(SCITBX_MATH_GCD_USING_ASM) && defined(__x86_64__)
  inline
  long
  gcd_int64_asm(
    long a,
    long b)
  {
    __asm__(
      "  movq %2, %%rcx\n"
      ".L_for_%=:\n"
      "  cmpq $0, %%rcx\n" // if b == 0
      "  je .L_return_a_or_minus_a_%=\n"
      "  movq %0, %%rdx\n"
      "  sarq $63, %%rdx\n"
      "  idivq %%rcx\n" // next_b = a % b
      "  movq %%rcx, %0\n" // a = b
      "  movq %%rdx, %%rcx\n" // b = next_b
      "  jmp .L_for_%=\n"
      ".L_return_a_or_minus_a_%=:\n"
         // next five lines: slightly faster than compare & jump
      "  movq %0, %%rcx\n"
      "  cqto\n"
      "  xorq %%rdx, %%rcx\n"
      "  subq %%rdx, %%rcx\n"
      "  movq %%rcx, %0\n"
        : "=a"(a) /* output */
        : "0"(a), "r"(b) /* input */
        : "cc", "%rcx", "%rdx"); /* clobbered registers */
    return a;
  }
#endif

  // from boost/math/common_factor_rt.hpp, svn trunk rev. 47847
  inline
  unsigned long
  gcd_unsigned_long_binary(
    unsigned long u,
    unsigned long v)
  {
    if ( u && v ) {
      // Shift out common factors of 2
      unsigned shifts = 0;
      while ( !(u & 1u) && !(v & 1u) ) {
        ++shifts;
        u >>= 1;
        v >>= 1;
      }
      // Start with the still-even one, if any
      unsigned long r[] = { u, v };
      unsigned which = static_cast<bool>( u & 1u );
      // Whittle down the values via their differences
      do {
        // Remove factors of two from the even one
        while ( !(r[ which ] & 1u) ) {
          r[ which ] >>= 1;
        }
        // Replace the larger of the two with their difference
        if ( r[!which] > r[which] ) {
          which ^= 1u;
        }
        r[ which ] -= r[ !which ];
      }
      while ( r[which] );
      // Shift-in the common factor of 2 to the residues' GCD
      return r[ !which ] << shifts;
    }
    else {
      // At least one input is zero, return the other
      // (adding since zero is the additive identity)
      // or zero if both are zero.
      return u + v;
    }
  }

  inline
  long
  gcd_long_binary(
    long u,
    long v)
  {
    return static_cast<long>(
      gcd_unsigned_long_binary(
        u < 0 ? -u : u,
        v < 0 ? -v : v));
  }

  inline
  int
  gcd_int(
    int a,
    int b)
  {
#if defined(SCITBX_MATH_GCD_USING_ASM)
    return gcd_int32_asm(a, b);
#else
    return gcd_int_simple(a, b);
#endif
  }

  inline
  long
  gcd_long(
    long a,
    long b)
  {
#if defined(SCITBX_MATH_GCD_USING_ASM)
# if defined(__x86_64__)
    return gcd_int64_asm(a, b);
# else
    return gcd_int32_asm(a, b);
# endif
#else
    return gcd_long_simple(a, b);
#endif
  }

}} // namespace scitbx::math

#endif // GUARD
