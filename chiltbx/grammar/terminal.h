/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_TERMINAL_HPP
#define CHILTBX_GRAMMAR_TERMINAL_HPP

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename > struct range_tag {};
template < typename T > struct unbox< range_tag<T> > {
        typedef typename unbox<T>::type type;
};
template < typename > struct container_tag {};
template < typename T > struct unbox< container_tag<T> > {
        typedef typename unbox<T>::type type;
};

template < typename T >
struct terminal : public state<true,T> {
        terminal ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<T> parse ( iterator& itr, iterator const& ind, parser& self ) {
                if ( itr == ind )
                        return result<T>(false);
                if ( *itr != this->left() )
                        return result<T>(false);
                result<T> R(true,*itr);
                ++itr;
                return R;
        }
        template < typename rtype, typename iterator >
        result<T> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename T >
struct terminal< range_tag<T> > : public state<true,T,T> {
        terminal ( T const& t1, T const& t2 ) : state<true,T,T>(t1,t2) {}
        template < typename rtype, typename iterator, typename parser >
        result<T> parse ( iterator& itr, iterator const& end, parser& self ) {
                if ( itr == end )
                        return result<T>(false);
                if ( this->left() > *itr
                        || *itr > this->right() )
                        return result<T>(false);
                result<T> R(true,*itr);
                ++itr;
                return R;
        }
        template < typename rtype, typename iterator >
        result<T> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename T >
struct terminal< container_tag<T> > : public state<true,T> {
        terminal ( T const& t ) : state<true,T>(t) {}
        template < typename rtype, typename iterator, typename parser >
        result<T> parse ( iterator& itr, iterator const& end, parser& self ) {
                if ( itr == end )
                        return result<T>(false);
                typename T::const_iterator ctr = this->left().begin();
                typename T::const_iterator cnd = this->left().end();
                while ( itr != end && ctr != cnd )
                        if ( *itr != *ctr )
                                break;
                        else {
                                ++itr;
                                ++ctr;
                        }
                if ( ctr == cnd )
                        return result<T>(true,this->left());
                return result<T>(false);
        }
        template < typename rtype, typename iterator >
        result<T> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename T >
terminal<T> single ( T const& value ) {
        return terminal<T>(value);
}

template < typename T >
terminal< range_tag<T> > range ( T const& lower, T const& upper ) {
        return terminal< range_tag<T> >(lower,upper);
}

template < typename T >
terminal< container_tag<T> > container ( T const& value ) {
        return terminal< container_tag<T> >(value);
}

template < typename rtype, typename T, typename parser >
struct result_type< rtype, terminal<T>, parser > {
        typedef typename unbox<T>::type type;
};

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_TERMINAL_HPP
