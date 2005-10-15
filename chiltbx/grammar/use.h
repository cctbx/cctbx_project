/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_USE_H
#define CHILTBX_GRAMMAR_USE_H

#include<chiltbx/grammar/base.h>
#include<chiltbx/grammar/search.h>

namespace chiltbx {

namespace grammar {

template < typename Tag >
struct go_to {
        typedef go_to<Tag> self_type;

        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,typename search<parser,Tag>::type,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,typename search<parser,Tag>::type,parser>::type resultant;
                return search<parser,Tag>::get(self).parse<rtype>(itr,end,self);
        }

        template < typename rtype, typename iterator >
        result<null> parse ( iterator&, iterator const& ) {
                return result<null>(false);
        }
};

template < typename rtype, typename Tag, typename parser >
struct result_type<rtype,go_to<Tag>,parser> {
        typedef typename result_type<rtype,typename search<parser,Tag>::type,parser>::type type;
};

template < typename T >
go_to<T> use () {
        return go_to<T>();
}

template < iterator_t I >
go_to<iterator<I> > use () {
        return go_to<iterator<I> >();
}

}// end grammar namespace

}// end chiltbx

#endif//CHILTBX_GRAMMAR_USE_H
