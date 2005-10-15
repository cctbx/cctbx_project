/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_OPTIONAL_HPP
#define CHILTBX_GRAMMAR_OPTIONAL_HPP

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename T >
struct optional : public state<true,T> {
        typedef optional<T> self_type;
        optional ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,T,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                iterator tmp = itr;
                result<result_type> R = this->left().parse<rtype>(tmp,end,self);
                if ( not R.truth )
                        return result<result_type>(true);
                itr = tmp;
                return R;
        }
        template < typename rtype, typename iterator >
        result<typename result_type<rtype,T,self_type>::type>
        parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename parser >
struct result_type<rtype,optional<T>,parser> {
        typedef typename result_type<rtype,T,parser>::type type;
};

template < typename T >
struct has_left_branch< optional<T> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T >
optional<T> q ( T const& t ) {
        return optional<T>(t);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_OPTIONAL_HPP
