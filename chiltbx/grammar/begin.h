/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_BEGIN_H
#define CHILTBX_GRAMMAR_BEGIN_H

#include<chiltbx/grammar/base.h>
#include<chiltbx/grammar/search.h>

namespace chiltbx {

namespace grammar {

template < typename T, typename Tag >
struct begin : public state<true,T> {
        typedef begin<T,Tag> self_type;
        begin ( T const& t )
        : state<true,T>(t)
        , parser_( grammar::search<T,Tag>::get(this->left()) ) {}
        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,typename search<T,Tag>::type,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,typename search<T,Tag>::type,parser>::type result_type;
                if ( itr == end )
                        return result<result_type>(false);
                return this->parser_.parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<typename result_type<rtype,typename search<T,Tag>::type,self_type>::type>
        parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
private:
        typename grammar::search<T,Tag>::type& parser_;
};

template < typename T, typename G >
struct has_left_branch< begin<T,G> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename rtype, typename T, typename Tag, typename parser >
struct result_type<rtype,begin<T,Tag>,parser> {
        typedef typename result_type<rtype,typename search<T,Tag>::type,parser>::type type;
};

template < typename G, typename T >
begin<T,G> S ( T const& t ) {
        return begin<T,G>(t);
}

template < iterator_t I, typename T >
begin<T,iterator<I> > S ( T const& t ) {
        return begin<T,iterator<I> >(t);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_BEGIN_H
