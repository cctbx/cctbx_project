/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_PASS_THROUGH_H
#define CHILTBX_GRAMMAR_PASS_THROUGH_H

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename T >
struct pass_through : public state<true,T,bool> {
        typedef pass_through<T> self_type;
        pass_through ( T const& t, bool pass=true ) : state<true,T,bool>(t,pass) {}
        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,T,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                if ( itr == end )
                        return result<result_type>(false);
                if ( this->right() )
                        return result<result_type>(true);
                return this->left().parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<typename result_type<rtype,T,self_type>::type>
        parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename parser >
struct result_type<rtype,pass_through<T>,parser> {
        typedef typename result_type<rtype,T,parser>::type type;
};

template < typename T >
struct has_left_branch< pass_through<T> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T >
pass_through<T> skip ( T const& t, bool B=true ) {
        return pass_through<T>(t,B);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_PASS_THROUGH_H
