#ifndef CCTBX_ARRAY_FAMILY_ALIGNMENT_CALCULATOR_H
#define CCTBX_ARRAY_FAMILY_ALIGNMENT_CALCULATOR_H

#include <cstddef>
#include <boost/config.hpp>

namespace cctbx { namespace af { namespace alignment_calculator {

////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2001-2002 by Andrei Alexandrescu
// Reference:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
// Permission to use, copy, modify, distribute and sell this software for any
//     purpose is hereby granted without fee, provided that the above copyright
//     notice appear in all copies and that both that copyright notice and this
//     permission notice appear in supporting documentation.
// The author or Addison-Welsey Longman make no representations about the
//     suitability of this software for any purpose. It is provided "as is"
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// macros CCTBX_AF_TYPELIST1, TYPELIST_2, ... TYPELIST_50
// Each takes a number of arguments equal to its numeric suffix
// The arguments are type names. CCTBX_AF_TYPELISTNN generates a typelist containing
//     all types passed as arguments, in that order.
// Example: CCTBX_AF_TYPELIST2(char, int) generates a type containing char and int.
////////////////////////////////////////////////////////////////////////////////

#define CCTBX_AF_TYPELIST1(T1) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, ::cctbx::af::alignment_calculator::NullType>

#define CCTBX_AF_TYPELIST2(T1, T2) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST1(T2) >

#define CCTBX_AF_TYPELIST3(T1, T2, T3) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST2(T2, T3) >

#define CCTBX_AF_TYPELIST4(T1, T2, T3, T4) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST3(T2, T3, T4) >

#define CCTBX_AF_TYPELIST5(T1, T2, T3, T4, T5) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST4(T2, T3, T4, T5) >

#define CCTBX_AF_TYPELIST6(T1, T2, T3, T4, T5, T6) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST5(T2, T3, T4, T5, T6) >

#define CCTBX_AF_TYPELIST7(T1, T2, T3, T4, T5, T6, T7) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST6(T2, T3, T4, T5, T6, T7) >

#define CCTBX_AF_TYPELIST8(T1, T2, T3, T4, T5, T6, T7, T8) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST7(T2, T3, T4, T5, T6, T7, T8) >

#define CCTBX_AF_TYPELIST9(T1, T2, T3, T4, T5, T6, T7, T8, T9) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST8(T2, T3, T4, T5, T6, T7, T8, T9) >

#define CCTBX_AF_TYPELIST10(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST9(T2, T3, T4, T5, T6, T7, T8, T9, T10) >

#define CCTBX_AF_TYPELIST11(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST10(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11) >

#define CCTBX_AF_TYPELIST12(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST11(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12) >

#define CCTBX_AF_TYPELIST13(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST12(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13) >

#define CCTBX_AF_TYPELIST14(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST13(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14) >

#define CCTBX_AF_TYPELIST15(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST14(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15) >

#define CCTBX_AF_TYPELIST16(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST15(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16) >

#define CCTBX_AF_TYPELIST17(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST16(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17) >

#define CCTBX_AF_TYPELIST18(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST17(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18) >

#define CCTBX_AF_TYPELIST19(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST18(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19) >

#define CCTBX_AF_TYPELIST20(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST19(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20) >

#define CCTBX_AF_TYPELIST21(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST20(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21) >

#define CCTBX_AF_TYPELIST22(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST21(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22) >

#define CCTBX_AF_TYPELIST23(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST22(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23) >

#define CCTBX_AF_TYPELIST24(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST23(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24) >

#define CCTBX_AF_TYPELIST25(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST24(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25) >

#define CCTBX_AF_TYPELIST26(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST25(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26) >

#define CCTBX_AF_TYPELIST27(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST26(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27) >

#define CCTBX_AF_TYPELIST28(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST27(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28) >

#define CCTBX_AF_TYPELIST29(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST28(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29) >

#define CCTBX_AF_TYPELIST30(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST29(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30) >

#define CCTBX_AF_TYPELIST31(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST30(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31) >

#define CCTBX_AF_TYPELIST32(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST31(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, \
        T32) >

#define CCTBX_AF_TYPELIST33(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, \
        T32, T33) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST32(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, \
        T32, T33) >

#define CCTBX_AF_TYPELIST34(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST33(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34) >

#define CCTBX_AF_TYPELIST35(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST34(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35) >

#define CCTBX_AF_TYPELIST36(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST35(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36) >

#define CCTBX_AF_TYPELIST37(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37) \
        ::cctbx::af::alignment_calculator::Typelist< \
        T1, CCTBX_AF_TYPELIST36(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37) >

////////////////////////////////////////////////////////////////////////////////
// class template Typelist
// The building block of typelists of any length
// Use it through the CCTBX_AF_TYPELISTNN macros
// Defines nested types:
//     Head (first element, a non-typelist type by convention)
//     Tail (second element, can be another typelist)
////////////////////////////////////////////////////////////////////////////////

    class NullType {};

    template <class T, class U>
    struct Typelist
    {
       typedef T Head;
       typedef U Tail;
    };

////////////////////////////////////////////////////////////////////////////////
// class template Select
// Selects one of two types based upon a boolean constant
// Invocation: Select<flag, T, U>::Result
// where:
// flag is a compile-time boolean constant
// T and U are types
// Result evaluates to T if flag is true, and to U otherwise.
////////////////////////////////////////////////////////////////////////////////

    template <bool flag, typename T, typename U>
    struct Select
    {
        typedef T Result;
    };
    template <typename T, typename U>
    struct Select<false, T, U>
    {
        typedef U Result;
    };

////////////////////////////////////////////////////////////////////////////////
// class template MaxSize
// Computes the maximum sizeof for all types in a typelist
// Usage: MaxSize<TList>::result
////////////////////////////////////////////////////////////////////////////////

template <class TList> struct MaxSize;

template <>
struct MaxSize<NullType>
{
    enum { result = 0 };
};

template <class Head, class Tail>
struct MaxSize< Typelist<Head, Tail> >
{
private:
    enum { tailResult = std::size_t(MaxSize<Tail>::result) };
public:
    enum { result = sizeof(Head) > tailResult ?
        sizeof(Head) : std::size_t(tailResult) };
};

////////////////////////////////////////////////////////////////////////////////
// class template Align
// Builds a union that contains each type in a typelist
// Usage: Align<TList> is the very type
////////////////////////////////////////////////////////////////////////////////

template <class U> union Align;

template <> union Align<NullType>
{
};

template <class Head, class Tail>
union Align< Typelist<Head, Tail> >
{
   Head head_;
   Align<Tail> tail_;
};

////////////////////////////////////////////////////////////////////////////////
// class AlignmentCalculatorBase
// Used by AlignmentCalculator
////////////////////////////////////////////////////////////////////////////////

class AlignmentCalculatorBase
{
protected:
    template <class TList, std::size_t size> struct Compute;

    template <std::size_t size>
    struct Compute<NullType, size>
    {
        typedef NullType Result;
    };

    template <class Head, class Tail, std::size_t size>
    struct Compute<Typelist<Head, Tail>, size>
    {
       typedef typename Compute<Tail, size>::Result TailResult;

       typedef typename Select<
               sizeof(Head) <= size,
               Typelist<Head, TailResult>,
               TailResult>::Result
           Result;
    };

    template <typename U> struct Structify
    { U dummy_; };

    class Unknown;
    struct Abstract
    {
#if !(defined(__EDG_VERSION__) && __EDG_VERSION__ <= 245)
#if !(defined(__GNUC__) && ((__GNUC__ < 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 2))))
       virtual void Fun() {}
#endif
#endif
    };

    typedef CCTBX_AF_TYPELIST37(
            char,
            short int,
            int,
            long int,
            float,
            double,
            long double,
            char*,
            short int*,
            int*,
            long int*,
            float*,
            double*,
            long double*,
            void*,
            Unknown (*)(Unknown),
            Unknown* Unknown::*,
            Unknown (Unknown::*)(Unknown),
            Abstract,
            Structify<char>,
            Structify<short int>,
            Structify<int>,
            Structify<long int>,
            Structify<float>,
            Structify<double>,
            Structify<long double>,
            Structify<char*>,
            Structify<short int*>,
            Structify<int*>,
            Structify<long int*>,
            Structify<float*>,
            Structify<double*>,
            Structify<long double*>,
            Structify<void*>,
            Structify<Unknown (*)(Unknown)>,
            Structify<Unknown* Unknown::*>,
            Structify<Unknown (Unknown::*)(Unknown)>
            )
        TypesOfAllAlignments;
};

////////////////////////////////////////////////////////////////////////////////
// class template AlignmentCalculator
// Computes the alignment of all types in a typelist
// Usage: Align<TList> is the very type
////////////////////////////////////////////////////////////////////////////////

template <typename TList>
class AlignmentCalculator : private AlignmentCalculatorBase
{
    enum { maxSize = MaxSize<TList>::result };
    typedef typename Compute<TypesOfAllAlignments, maxSize>::Result
        AlignTypes;
public:
    typedef Align<AlignTypes> Result;
};

}}} // namespace cctbx::af::alignment_calculator

#endif // CCTBX_ARRAY_FAMILY_ALIGNMENT_CALCULATOR_H
