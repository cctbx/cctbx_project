/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_EITHER_HPP
#define CHILTBX_GRAMMAR_EITHER_HPP

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename T, typename U, typename Disam >
struct either : public state<true,T,U>, Disam {
        either ( T const& t, U const& u ) : state<true,T,U>(t,u) {}
        template < typename rtype, typename iterator, typename parser >
        result<rtype> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type left_type;
                typedef typename result_type<rtype,U,parser>::type right_type;
                if ( itr == end )
                        return result<rtype>(false);
                iterator ltr = itr;
                iterator rtr = itr;
                result<left_type> L = this->left().parse<rtype>(ltr,end,self);
                result<right_type> R = this->right().parse<rtype>(rtr,end,self);
                result<rtype> Q(true);
                if ( L.truth and not R.truth ) {
                        itr = ltr;
                        Q.value |= L.value;
                        return Q;
                } else if ( R.truth and not L.truth ) {
                        itr = rtr;
                        Q.value |= R.value;
                        return Q;
                } else if ( not L.truth and not R.truth )
                        return result<rtype>(false);
                return this->template disambiguate<rtype>(ltr,L,rtr,R,itr,end);
                return result<rtype>(false);
        }
        template < typename rtype, typename iterator >
        result<rtype> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename T, typename U, typename D >
struct has_left_branch< either<T,U,D> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T, typename U, typename D >
struct has_right_branch< either<T,U,D> > {
        typedef typename unbox<U>::type type;
        static const_truth_t value = true;
};

struct longest {
        template < typename type, typename iterator, typename ltype, typename rtype >
        result<type> disambiguate ( iterator& ltr, ltype& L,
                                                iterator& rtr, rtype& R,
                                                iterator& itr, iterator const& end ) const {
                std::size_t len = end - ltr;
                std::size_t ren = end - rtr;
                result<type> Q(true);
                if ( len < ren ) {
                        itr = ltr;
                        Q.value |= L.value;
                        return Q;
                } else if ( ren < len ) {
                        itr = rtr;
                        Q.value |= R.value;
                        return Q;
                }
                return result<type>(false);
        }
};

struct barf {
        template < typename type, typename iterator, typename ltype, typename rtype >
        result<type> disambiguate ( iterator&, ltype&,
                                                iterator&, rtype&,
                                                iterator&, iterator const& ) {
                throw grammar_error();
        }
};

struct first {
        template < typename type, typename iterator, typename ltype, typename rtype >
        result<type> disambiguate ( iterator& ltr, ltype& L,
                                                iterator&, rtype&,
                                                iterator& itr, iterator const& ) const {
                itr = ltr;
                result<type> R(true);
                R.value |= L.value;
                return R;
        }
};

struct exclusive {
        template < typename type, typename iterator, typename ltype, typename rtype >
        result<type> disambiguate ( iterator& ltr, ltype&,
                                                iterator& rtr, rtype&,
                                                iterator& itr, iterator const& ) const {
                throw result<type>(false);
        }
};

template < typename T, typename U >
either<T,U,barf> vomit ( T const& t, U const& u ) {
        return either<T,U,barf>(t,u);
}

template < typename T, typename U >
either<T,U,longest> operator | ( T const& t, U const& u ) {
        return either<T,U,longest>(t,u);
}

template < typename T, typename U >
either<T,U,first> operator , ( T const& t, U const& u ) {
        return either<T,U,first>(t,u);
}

template < typename T, typename U >
either<T,U,exclusive> operator ^ ( T const& t, U const& u ) {
        return either<T,U,exclusive>(t,u);
}

}//end grammar

}//end chiltbx

#endif//CHILTBX_GRAMMAR_EITHER_HPP
