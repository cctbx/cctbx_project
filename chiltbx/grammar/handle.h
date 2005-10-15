/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_HANDLE_H
#define CHILTBX_GRAMMAR_HANDLE_H

#include<chiltbx/grammar/base.h>
#include<chiltbx/handle.h>

namespace chiltbx {

namespace grammar {

namespace ns_handle = chiltbx::handle;

template <      typename RType,
                typename Iterator,
                typename T=void
                                        > struct handle;

template < typename RType, typename Iterator >
struct handle<RType,Iterator,ns_handle::interface> {
        virtual ~ handle () {}
        virtual RType parse ( Iterator&, Iterator const& ) = 0;
};

template < typename RType, typename Iterator, typename T >
struct handle : public state<true,T>, public handle<RType,Iterator,void> {
        handle ( T const& t ) : state<true,T>(t) {}
        virtual ~ handle () {}
        template < typename rtype, typename iterator, typename parser >
        result<RType> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                if ( itr == end )
                        return rtype(false);
                result<result_type> R = this->left().parse<RType>(itr,end,self);
                if ( not R.truth )
                        return result<RType>(false);
                return result<RType>(R);
        }
        virtual result<RType> parse ( Iterator& itr, Iterator const& end ) {
                return this->parse<RType>(itr,end,*this);
        }
};

template < typename RType, typename Iterator >
struct handle<RType,Iterator,void> {
public:
        typedef handle<RType,Iterator,ns_handle::interface> h_parser;
        handle () {}
        template < typename T >
        handle ( T const& t ) {
                this->handle_.set( handle<RType,Iterator,T>(t) );
        }
        template < typename rtype, typename parser >
        result<RType> parse ( Iterator& itr, Iterator const& end, parser& self ) {
                if ( itr == end )
                        return rtype(false);
                result<RType> R = this->handle_.get_raw()->parse(itr,end);
                if ( not R )
                        return result<RType>(false);
                return R;
        }
        template < typename rtype >
        result<RType> parse ( Iterator& itr, Iterator const& end ) {
                return this->parse<rtype>(itr,end);
        }
private:
        ns_handle::handle<h_parser,ns_handle::by_value> handle_;
};

template < typename rtype, typename RType, typename Iterator, typename parser >
struct result_type<rtype,handle<RType,Iterator>,parser> {
        typedef RType type;
};

template < typename R, typename I, typename T >
handle<R,I> abstract ( T const& t ) {
        return handle<R,I>(t);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_HANDLE_H
