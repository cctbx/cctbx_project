/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_SEARCH_H
#define CHILTBX_GRAMMAR_SEARCH_H

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

struct left {};
struct right {};
struct here {};
template < typename > struct look {};

template < typename T, typename Tag=void >
struct search : state<true,T> {

        typedef search<T,Tag> self_type;

        typedef typename if_< has_tag<T,Tag>::value,
                typename has_tag<T,Tag>::type,
                typename if_< has_tag<typename has_left_branch<T>::type,Tag>::value,
                        typename has_tag<typename has_left_branch<T>::type,Tag>::type,
                        typename if_< has_tag<typename has_right_branch<T>::type,Tag>::value,
                                typename has_tag<typename has_right_branch<T>::type,Tag>::type,
                                typename if_< not is_null<typename search<typename has_left_branch<T>::type,Tag>::type>::value,
                                        typename search<typename has_left_branch<T>::type,Tag>::type,
                                        typename if_< not is_null<typename search<typename has_right_branch<T>::type,Tag>::type>::value,
                                                typename search<typename has_right_branch<T>::type,Tag>::type,
                                                null
                                        >::type
                                >::type
                        >::type
                >::type
        >::type type;

        search ( T const& value ) : state<true,T>(value), parser_(search<T,Tag>::get(this->left())) {}

        template < typename rtype, typename iterator, typename parser >
        result<typename result_type<rtype,type,parser>::type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,type,parser>::type result_type;
                if ( itr == end )
                        return result<result_type>(false);
                return this->parser_.parse<typename rtype>(itr,end,self);
        }

        template < typename rtype, typename iterator >
        result<typename result_type<rtype,type,self_type>::type>
        parse ( iterator& itr, iterator const& ind ) {
                return this->parse<typename rtype>(itr,ind,*this);
        }

        template < typename U >
        static type& get ( U& u, here const& ) {
                return u.left();
        }

        template < typename U >
        static type& get ( U& u, left const& ) {
                return search<T,Tag>::get(u.left(),here());
        }

        template < typename U >
        static type& get ( U& u, right const& ) {
                return search<T,Tag>::get(u.right(),here());
        }

        template < typename U >
        static type& get ( U& u, look<left> const& ) {
                return search<T,Tag>::get(u.left());
        }

        template < typename U >
        static type& get ( U& u, look<right> const& ) {
                return search<T,Tag>::get(u.right());
        }

        template < typename U >
        static type& get ( U& u, null const& ) {
                static null N;
                return N;
        }

        template < typename U >
        static type& get ( U& u ) {
                typedef
                typename if_< has_tag<U,Tag>::value,
                        here,
                        typename if_< has_tag<typename has_left_branch<U>::type,Tag>::value,
                                left,
                                typename if_< has_tag<typename has_right_branch<U>::type,Tag>::value,
                                        right,
                                        typename if_< not is_null<typename search<typename has_left_branch<U>::type,Tag>::type>::value,
                                                look<left>,
                                                typename if_< not is_null<typename search<typename has_right_branch<U>::type,Tag>::type>::value,
                                                        look<right>,
                                                        null
                                                >::type
                                        >::type
                                >::type
                        >::type
                >::type direction;

                return search<T,Tag>::get(u,direction());
        }

private:
        type &parser_;
};

template < typename Tag > struct search<null,Tag> {
        typedef null type;
        null result;
};

template < typename rtype, typename T, typename Tag, typename parser >
struct result_type<rtype,search<T,Tag>,parser> {
        typedef typename result_type<rtype,typename search<T,Tag>::type,parser>::type type;
};

template < typename Tag > struct search<Tag,void> {
        template < typename T >
        search<T,Tag> operator = ( T const& t ) {
                return search<T,Tag>(t);
        }
};

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_SEARCH_H
