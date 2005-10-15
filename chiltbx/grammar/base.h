/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_BASE_HPP
#define CHILTBX_GRAMMAR_BASE_HPP

#include<chiltbx/childef.h>
#include<string>
#include<sstream>
#include<fstream>

namespace chiltbx {

namespace grammar {

struct grammar_error {
        virtual const char* what () const {
                return "Grammar error.";
        }
};

// character based templatized type definitions
template < typename Char > struct                                               char_type {
        typedef std::basic_string<Char>                                         string;
        typedef std::basic_stringstream<Char>                           sstream;
        typedef std::basic_ifstream<Char>                                       ifstream;
        typedef std::basic_ofstream<Char>                                       ofstream;
};

// the result-type of a parser
template < typename T > struct                                                  result {
        result ( truth_t worked=false ) : truth(worked) {}
        result ( truth_t worked, T const& val )
        : value(val), truth(worked) {}
        T value;
        truth_t truth;
};

/* definition of null & is_null partialization template */
struct                                                                                                  null {
        typedef null self_type;
        template < typename rtype, typename iterator, typename parser >
        result<rtype> parse ( iterator&, iterator const&, parser& ) {
                return result<rtype>(false);
        }
};

template < typename T > struct                                                  is_null {
        typedef T type;
        static const_truth_t value = false;
};
template <> struct                                                                              is_null<null> {
        typedef null type;
        static const_truth_t value = true;
};

/* default branch-values */
template < typename T > struct                                                  has_left_branch {
        typedef null type;
        static const_truth_t value = false;
};
template < typename T > struct                                                  has_right_branch {
        typedef null type;
        static const_truth_t value = false;
};

/* compile time conditional */
template < const_truth_t, typename, typename > struct   if_;
template < typename L, typename R > struct                              if_<true,L,R> {
        typedef L type;
        static const_truth_t value = true;
};
template < typename L, typename R > struct                              if_<false,L,R> {
        typedef R type;
        static const_truth_t value = false;
};

/* resultant type from a parser */
template < typename T, typename specializer, typename parser >
struct                                                                                                  result_type {
        typedef T type;
};

/* whether or not a parser is "tagged" */
template < typename T, typename U > struct                              has_tag {
        typedef null type;
        static const_truth_t value = false;
};

/*
the state-types for use with parsers to make them contractually safe for searching
*/
template < truth_t, typename T=void, typename U=void >
struct                                                                                                  state {
        state ( T const& t, U const& u )
        : left_(t), right_(u) {}
        T const& left () const { return this->left_; }
        T& left () { return this->left_; }
        U const& right () const { return this->right_; }
        U& right () { return this->right_; }
private:
        T left_;
        U right_;
};

template < truth_t Truth, typename T > struct                   state<Truth,T,void> {
        state ( T const& t ) : left_(t) {}
        T const& left () const { return this->left_; }
        T& left () { return this->left_; }
private:
        T left_;
};

template < typename T, typename U > struct                              state<false,T,U> {
        typedef null pass_through;
        T left () const { return T(); }
        U right () const { return U(); }
};

template < typename T > struct                                                  state<false,T,void> {
        T left () const { return T(); }
};

template <> struct                                                                              state<false,void,void> {
};

// for unboxing "boxed" values
template < typename T > struct                                                  unbox {
        typedef T type;
};

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_BASE_HPP
