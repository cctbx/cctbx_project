#ifndef NOEXCEPT_FALSE_HPP
#define NOEXCEPT_FALSE_HPP

// Used for destructors that throw exceptions

// If using C++11 or a later version, or VS2015 and later
#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
#define NOEXCEPT_FALSE noexcept(false)
#else
#define NOEXCEPT_FALSE
#endif

#endif // GUARD
