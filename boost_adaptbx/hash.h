#ifndef BOOST_ADAPTBX_HASH_H
#define BOOST_ADAPTBX_HASH_H

namespace boost_adaptbx {

  /// A class to inherit from to provide Python hashing in a wrapper
  template <class WrappedType>
  struct py_hashable
  {
    /// The hash of the object as expected by Python
    /** Boost.Python converts std::size_t into a PyLongObject (an integer of
     arbitrary length) whereas the __hash__ function is expected to return
     a PyIntObject (a C long). This lead to potential Overflow exception
     when the PyLongObject does not fit into a PyIntObject.
     Since Python 2.5, the result of __hash__ is itself hashed into a
     PyIntObject if it is a PyLongObject, thus alleviating the problem.
     To retain compatibility with Python 2.3 and 2.4, we need to do that
     conversion ourselves.
     */
    static long py_hash(WrappedType const &self) {
      std::size_t h = boost::hash<WrappedType>()(self);
      long result;
      if (sizeof(std::size_t) == sizeof(long)) {
        // LP32 and LP64 platforms (Linux, MacOS)
        result = (long)h;
      }
      else {
        /* Prime example: 64-bit Windows which is LLP64 (long is 32 bits)
           This is the algorithm used in Boost.Hash,
           hash/hash.hpp, hash_detail::hash_value_unsigned,
           edited to produce an unsigned long hash instead of std::size_t
         */
        const int ulong_bits = std::numeric_limits<unsigned long>::digits;
        // ceiling(std::numeric_limits<T>::digits / ulong_bits) - 1
        const int length = (std::numeric_limits<std::size_t>::digits - 1)
           / ulong_bits;

        unsigned long seed = 0;

        // Hopefully, this loop can be unrolled.
        for(unsigned int i = length * ulong_bits; i > 0; i -= ulong_bits)
        {
           seed ^= (unsigned long) (h >> i) + (seed<<6) + (seed>>2);
        }
        seed ^= (unsigned long) h + (seed<<6) + (seed>>2);
        result = (long)seed;
      }
      /// Python convention:
      if (result == -1) result = -2;
      return result;
    }
  };

}

#endif // GUARD
