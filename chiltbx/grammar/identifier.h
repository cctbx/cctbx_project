/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_IDENTIFIER_HPP
#define CHILTBX_GRAMMAR_IDENTIFIER_HPP

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

struct identify_tag;

template < typename T, typename Tag >
struct identifier : public state<true,T> {
        typedef identifier<T,Tag> self_type;
        identifier ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,T,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                if ( itr == end )
                        return result<result_type>(false);
                return this->left().parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<typename result_type<rtype,T,self_type>::type>
        parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename Tag, typename parser >
struct result_type<rtype,identifier<T,Tag>,parser> {
        typedef typename result_type<rtype,T,parser>::type type;
};

template < typename T, typename G > struct has_left_branch< identifier<T,G> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T, typename Tag >
struct has_tag< identifier<T,Tag>, Tag > {
        typedef T type;
        static const_truth_t value = true;
};

template < typename G, typename T >
identifier<T,G> tag ( T const& t ) {
        return identifier<T,G>(t);
}

template < iterator_t I, typename T >
identifier<T,iterator<I> > tag ( T const& t ) {
        return identifier<T,iterator<I> >(t);
}


}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_IDENTIFIER_HPP
