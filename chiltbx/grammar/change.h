/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_CHANGE_HPP
#define CHILTBX_GRAMMAR_CHANGE_HPP

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename T, typename qtype >
struct change : public state<true,T> {
        change ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<qtype> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<qtype,T,parser>::type result_type;
                if ( itr == end )
                        return result<qtype>(false);
                result<result_type> Q = this->left().parse<qtype>(itr,end,self);
                if ( Q.truth ) {
                        result<qtype> R(true);
                        R.value = Q.value;
                        return R;
                }
                return result<qtype>(false);
        }
        template < typename rtype, typename iterator >
        result<qtype> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename qtype, typename parser >
struct result_type<rtype,change<T,qtype>,parser> {
        typedef qtype type;
};

template < typename T, typename Q > struct has_left_branch< change<T,Q> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename Q, typename T >
change<T,Q> d ( T const& t ) {
        return change<T,Q>(t);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_CHANGE_HPP
