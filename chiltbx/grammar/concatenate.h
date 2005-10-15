/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_GRAMMAR_CONCATENATE_H
#define CHILTBX_GRAMMAR_CONCATENATE_H

#include<chiltbx/grammar/base.h>

namespace chiltbx {

namespace grammar {

template < typename T, typename U >
struct concatenate : public state<true,T,U> {
        concatenate ( T const& left, U const& right ) : state<true,T,U>(left,right) {}
        template < typename rtype, typename iterator, typename parser >
        result<rtype> parse ( iterator& itr, iterator const& ind, parser& self ) {
                typedef typename result_type<rtype,T,parser>::type left_type;
                typedef typename result_type<rtype,U,parser>::type right_type;
                iterator tmp = itr;
                result<left_type> lresult = this->left().parse<rtype>(tmp,ind,self);
                if ( not lresult.truth )
                        return result<rtype>(false);
                result<right_type> rresult = this->right().parse<rtype>(tmp,ind,self);
                if ( not rresult.truth )
                        return result<rtype>(false);
                result<rtype> R(true);
                itr = tmp;
                R.value >>= std::make_pair(lresult.value,rresult.value);
                return R;
        }
        template < typename rtype, typename iterator >
        result<rtype> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename T, typename U > struct has_left_branch<concatenate<T,U> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T, typename U > struct has_right_branch<concatenate<T,U> > {
        typedef typename unbox<T>::type type;
        static const_truth_t value = true;
};

template < typename T, typename U >
concatenate<T,U> operator >> ( T const& left, U const& right ) {
        return concatenate<T,U>(left,right);
}

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_CONCATENATE_H
