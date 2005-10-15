/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_PRODUCTION_HPP
#define CHILTBX_GRAMMAR_PRODUCTION_HPP

#include<chiltbx/grammar/base.h>
#include<chiltbx/grammar/change.h>

namespace chiltbx {

namespace grammar {

struct law_tag;
struct rule_tag;
template < typename, typename > struct pair_of_tags;

template < typename T, typename Tag=void >
struct production : public state<true,T> {
        typedef production<T,Tag> self_type;
        production ( T const& t ) : state<true,T>(t) {}
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

template < typename Tag > struct production<Tag,rule_tag> {
        template < typename T >
        production<T,Tag> operator = ( T const& t ) {
                return production<T,Tag>(t);
        }
};

template < typename Tag > struct production<Tag,law_tag> {
        template < typename T >
        production<change<T,Tag>,Tag> operator = ( T const& t ) {
                return production<change<T,Tag>,Tag>( d<Tag>(t) );
        }
};

template < typename Tag, typename RType >
struct production<pair_of_tags<Tag,RType>,law_tag> {
        template < typename T >
        production<change<T,RType>,Tag> operator = ( T const& t ) {
                return production<change<T,RType>,Tag>( d<RType>(t) );
        }
};

template < typename T, typename G > struct has_left_branch< production<T,G> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T, typename Tag >
struct has_tag< production<T,Tag>, Tag > {
        typedef T type;
        static const_truth_t value = true;
};

template < typename Tag >
production<Tag,rule_tag> rule () {
        return production<Tag,rule_tag>();
};

template < const_iterator_t I >
production<iterator<I>,rule_tag> rule () {
        return production<iterator<I>,rule_tag>();
};


template < typename Tag >
production<Tag,law_tag> law () {
        return production<Tag,law_tag>();
};

template < typename Tag, typename RType >
production<pair_of_tags<Tag,RType>,law_tag> law () {
        return production<pair_of_tags<Tag,RType>,law_tag>();
};

template < const_iterator_t I, typename RType >
production<pair_of_tags<iterator<I>,RType>,law_tag> law () {
        return production<pair_of_tags<iterator<I>,RType>,law_tag>();
};


}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_PRODUCTION_HPP
