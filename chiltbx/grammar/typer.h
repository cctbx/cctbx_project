/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_TYPER_HPP
#define CHILTBX_GRAMMAR_TYPER_HPP

#include<chiltbx/grammar/base.h>
#include<chiltbx/grammar/change.h>

namespace chiltbx {

namespace grammar {

template < typename T, typename Tag >
struct typer : public state<true,T> {
        typer ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<Tag> parse ( iterator& itr, iterator const& end, parser& self ) {
                if ( itr == end )
                        return result<Tag>(false);
                if ( this->left().parse<rtype>(itr,end,self).truth )
                        return result<Tag>(true);
                return result<Tag>(false);
        }
        template < typename rtype, typename iterator >
        result<Tag> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename Tag, typename parser >
struct result_type<rtype,typer<T,Tag>,parser> {
        typedef Tag type;
};

template < typename T, typename Tag > struct has_left_branch< typer<T,Tag> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename Tag, typename T >
typer<T,Tag> type ( T const& t ) {
        return typer<T,Tag>(t);
}

template < const_iterator_t I, typename T >
typer<T,iterator<I> > type ( T const& t ) {
        return typer<T,iterator<I> >(t);
}

template < const_iterator_t I, typename RType, typename T >
typer<change<T,RType>,iterator<I> > type ( T const& t ) {
        return typer<change<T,RType>,iterator<I> >( d<RType>(t) );
}

}// end grammar namespace

}// end chiltbx

#endif//CHILTBX_GRAMMAR_TYPER_HPP
