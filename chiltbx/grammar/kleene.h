/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_KLEENE_H
#define CHILTBX_GRAMMAR_KLEENE_H

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

struct infinite_loop : public grammar_error {
        virtual const char* what () const {
                return "Kleene-star detected an infinite loop.";
        }
};

template<typename> struct bounded;
template < typename T > struct unbox< bounded<T> > {
        typedef typename unbox<T>::type type;
};

template < typename T >
struct kleene : public state<true,T,iterator_t> {
        kleene ( T const& value, iterator_t const& itr=iterator_t(-1) )
        : state<true,T,iterator_t>(value,itr) {}
        template < typename rtype, typename iterator, typename parser >
        result<rtype> parse ( iterator& itr, iterator const& ind, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                if ( itr == ind )
                        return result<rtype>(false);
                result<rtype> R(true);
                iterator_t repetitions = 0;
                iterator tmp = itr;
                while ( tmp != ind ) {
                        iterator loc = tmp;
                        result<result_type> Q = this->left().parse<rtype>(tmp,ind,self);
                        if ( not Q.truth )
                                break;
                        if ( tmp == loc )
                                throw infinite_loop();
                        ++repetitions;
                        R.value *= Q.value;
                }
                if ( iterator_t(-1) != this->right() )
                        if ( repetitions < this->right() )
                                return result<rtype>(false);
                itr = tmp;
                return R;
        }
        template < typename rtype, typename iterator >
        result<rtype> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename T >
struct kleene< bounded<T> > : public state<true,T,iterator_t> {
        kleene ( T const& value, iterator_t const& itr, bool exact=false )
        : state<true,T,iterator_t>(value,itr), exactly_(exact) {}
        template < typename rtype, typename iterator, typename parser >
        result<rtype> parse ( iterator& itr, iterator const& ind, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type result_type;
                if ( itr == ind )
                        return result<rtype>(false);
                result<rtype> R(true);
                iterator tmp = itr;
                iterator_t repetitions = 0;
                for ( iterator_t i=0; i<this->right(); ++i ) {
                        iterator loc = tmp;
                        result<result_type> Q = this->left().parse<rtype>(tmp,ind,self);
                        if ( loc == tmp )
                                break;
                        if ( not Q.truth )
                                break;
                        ++repetitions;
                        R.value *= Q.value;
                }
                if ( this->exactly_ )
                        if ( this->right() != repetitions )
                                return result<rtype>(false);
                itr = tmp;
                return R;
        }
        template < typename rtype, typename iterator >
        result<rtype> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
        bool exactly_;
};

template < typename T > struct has_left_branch< kleene<T> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T >
kleene<T> operator * ( T const& value ) {
        return kleene<T>(value);
}

template < typename T >
kleene<T> operator + ( T const& value ) {
        return kleene<T>(value,1);
}

template < typename T >
kleene<T> operator >= ( T const& value, iterator_t const& itr ) {
        return kleene<T>(value,itr);
}

template < typename T >
kleene< bounded<T> > operator <= ( T const& value, iterator_t const& itr ) {
        return kleene< bounded<T> >(value,itr);
}

template < typename T >
kleene< bounded<T> > operator *= ( T const& value, iterator_t const& itr ) {
        return kleene< bounded<T> >(value,itr,true);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_KLEENE_H
