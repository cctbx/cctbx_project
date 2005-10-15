/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_SELF_HPP
#define CHILTBX_GRAMMAR_SELF_HPP

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename T >
struct self : public state<true,T> {
        typedef self<T> self_type;
        self ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,T,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& zelf ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                if ( itr == end )
                        return result<result_type>(false);
                return zelf.parse<rtype>(itr,end,zelf);
        }
        template < typename rtype, typename iterator >
        result<typename result_type<rtype,T,self_type>::type>
        parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename parser >
struct result_type<rtype,self<T>,parser> {
        typedef typename result_type<rtype,T,parser>::type type;
};

template < typename T > struct has_left_branch< self<T> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T >
self<T> recurse ( T const& t ) {
        return self<T>(t);
}

}// end grammar namespace

}// end chiltbx

#endif//CHILTBX_GRAMMAR_SELF_HPP
