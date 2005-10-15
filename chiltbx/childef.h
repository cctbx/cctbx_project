/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_STDDEF_HPP
#define CHILTBX_STDDEF_HPP

#include<cstddef>

namespace chiltbx {

typedef std::size_t             iterator_t;
typedef const iterator_t        const_iterator_t;
typedef bool                    truth_t;
typedef const truth_t           const_truth_t;

template < const_iterator_t I > struct iterator {
        typedef iterator<I> type;
        static const_iterator_t value = I;
};

template < const_truth_t T > struct truth {
        typedef truth<T> type;
        static const_truth_t value = T;
};

}// end chiltbx namespace

#endif//CHILTBX_STDDEF_HPP
