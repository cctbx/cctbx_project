/*
  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.
*/

#ifndef CHILTBX_GRAMMAR_STANDARD_PARSERS_H
#define CHILTBX_GRAMMAR_STANDARD_PARSERS_H

#include<chiltbx/grammar.h>

namespace chiltbx {

namespace grammar {

struct lower {};
struct upper {};
template < const_iterator_t Radix >
struct numeral {
        typedef numeral<Radix> type;
        static const_iterator_t value = Radix;
};
template < typename Char, Char Value >
struct character {
        typedef character<Char,Value> type;
        static const Char value = Value;
};

template < typename Case, typename Char=char > struct alpha;

template < typename rtype, typename T, typename Char, typename parser >
struct result_type<rtype,alpha<T,Char>,parser> {
        typedef Char type;
};

template < typename Char >
struct alpha<lower,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                return range('a','z').parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char >
struct alpha<upper,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                return range('A','Z').parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char >
struct alpha<void,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                iterator loc = itr;
                result<Char> r = alpha<upper,Char>().parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = alpha<lower,Char>().parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                return result<Char>(false);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char >
struct alpha<numeral<2>,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                return range('0','1').parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char >
struct alpha<numeral<8>,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                return range('0','7').parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char >
struct alpha<numeral<10>,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                return range('0','9').parse<rtype>(itr,end,self);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char >
struct alpha<numeral<16>,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                iterator loc = itr;
                result<Char> r = range('0','9').parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = range('a','f').parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = range('A','F').parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                return result<Char>(false);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char, Char Value >
struct alpha<character<Char,Value>,Char> {
        template < typename rtype, typename iterator, typename parser >
        result<Char>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                iterator tmp = itr;
                result<Char> R = single(Value).parse<rtype>(tmp,end,self);
                if ( not R.truth )
                        return result<Char>(false);
                itr = tmp;
                return result<Char>(true,Value);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename Char, Char Value, typename parser >
struct result_type<rtype,alpha<character<Char,Value>,Char>,parser> {
        typedef Char type;
};

template < typename Char, Char Value >
alpha<character<Char,Value>,Char> alpha_value () {
        return alpha<character<Char,Value>,Char>();
}

template < typename Char > struct alpha_values {
        static const Char under_score_                                          = '_';
        static const Char dot_                                                          = '.';
        static const Char equal_                                                        = '=';
        static const Char paren_                                                        = '(';
        static const Char thesis_                                                       = ')';
        static const Char bra_                                                          = '{';
        static const Char ket_                                                          = '}';
        static const Char car_                                                          = '<';
        static const Char ret_                                                          = '>';
        static const Char space_                                                        = ' ';
        static const Char tab_                                                          = '\t';
        static const Char eol_                                                          = '\n';
        typedef alpha<character<Char,under_score_>,Char>        under_score;
        typedef alpha<character<Char,dot_>,Char>                        dot;
        typedef alpha<character<Char,equal_>,Char>                      equal;
        typedef alpha<character<Char,paren_>,Char>                      paren;
        typedef alpha<character<Char,thesis_>,Char>                     thesis;
        typedef alpha<character<Char,bra_>,Char>                        bra;
        typedef alpha<character<Char,ket_>,Char>                        ket;
        typedef alpha<character<Char,space_>,Char>                      space;
        typedef alpha<character<Char,tab_>,Char>                        tab;
        typedef alpha<character<Char,eol_>,Char>                        eol;
};

template < typename Char >
struct white_space {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                iterator loc = itr;
                result<Char> r = typename alpha_values<Char>::space().parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = typename alpha_values<Char>::tab().parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = typename alpha_values<Char>::eol().parse<rtype>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                return result<Char>(false);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,white_space<Char>,parser> {
        typedef Char type;
};

template < typename Char >
struct c_alpha {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                iterator loc = itr;
                result<Char> r = alpha<void,Char>().parse<Char>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = typename alpha_values<Char>::under_score().parse<Char>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                return result<Char>(false);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<type>(itr,ind,*this);
        }
};

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,c_alpha<Char>,parser > {
        typedef Char type;
};

template < typename Char >
struct c_alphanumeric {
        template < typename rtype, typename iterator, typename parser >
        result<Char> parse ( iterator& itr, iterator const& end, parser& self ) {
                iterator loc = itr;
                result<Char> r = c_alpha<Char>().parse<Char>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                r = alpha<numeral<10>,Char>().parse<Char>(itr,end,self);
                if ( r.truth )
                        return r;
                itr = loc;
                return result<Char>(false);
        }
        template < typename rtype, typename iterator >
        result<Char> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<type>(itr,ind,*this);
        }
};

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,c_alphanumeric<Char>,parser> {
        typedef Char type;
};

template < typename String >
struct string {
        typedef typename String::value_type Char;
        typedef string<String> self_type;
        void operator |= ( Char const& chr ) {
                this->value_ += chr;
        }
        void operator |= ( String const& str ) {
                this->value_ += str;
        }
        void operator *= ( string<String> const& str ) {
                this->value_ += str.value_;
        }
        void operator *= ( Char const& chr ) {
                this->value_ += chr;
        }
        void operator >>= ( std::pair<Char,self_type> const& concat ) {
                this->value_ += concat.first;
                this->value_ += concat.second.value_;
        }
        inline std::size_t size () const { return this->value_.size(); }
        char const& operator [] ( std::size_t const& idx ) const {
                return this->value_[idx];
        }
        String value_;
};

template < typename Char >
struct c_word {
        typedef std::basic_string<Char> string_type;
        template < typename rtype, typename iterator, typename parser >
        result<string_type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef string<string_type> stracc;
                iterator tmp = itr;
                result<stracc> r = ( c_alpha<Char>() >>
                        q( *c_alphanumeric<Char>()) ).parse<stracc>(tmp,end,self);
                if ( not r.truth )
                        return result<string_type>(false);
                return result<string_type>(true,r.value.value_);
        }
        template < typename rtype, typename iterator >
        result<string_type> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,c_word<Char>,parser> {
        typedef std::basic_string<Char> type;
};

template < typename Parser >
struct eat : public state<true,Parser> {
        eat ( Parser const& P ) : state<true,Parser>(P) {}
        eat () : state<true,Parser>(Parser()) {}
        template < typename rtype, typename iterator, typename parser >
        result<null> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef typename result_type<rtype,Parser,parser>::type result_type;
                if ( end == itr )
                        return result<null>(false);
                iterator tmp = itr;
                result<result_type> r(false);
                do {
                        r = this->left().parse<rtype>(tmp,end,self);
                        if ( r.truth )
                                itr = tmp;
                } while ( r.truth && end != tmp );
                return result<null>(true);
        }
        template < typename rtype, typename iterator >
        result<null> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename T, typename parser >
struct result_type<rtype,eat<T>,parser> {
        typedef null type;
};

template < typename Char >
struct to_eol {
        typedef std::basic_string<Char> string_type;
        template < typename rtype, typename iterator, typename parser >
        result<string_type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                if ( end == itr )
                        return result<string_type>(false);
                iterator tmp = itr;
                result<Char> r(true);
                typename alpha_values<Char>::eol eol;
                string_type strt;
                bool found_eol = false;
                do {
                        r = eol.parse<Char>(tmp,end,self);
                        if ( not r.truth )
                                strt += *tmp++;
                        else
                                found_eol = true;
                } while ( not r.truth && end != tmp );
                if ( not found_eol )
                        return result<string_type>(false);
                itr = tmp;
                return result<string_type>(true,strt);
        }
        template < typename rtype, typename iterator >
        result<string_type> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,to_eol<Char>,parser> {
        typedef std::basic_string<Char> type;
};

template < typename Char >
int char_to_int ( Char const& C ) {
        if ( '0' <= C && C <= '9' )
                return C-'0';
        if ( 'a' <= C && C <= 'f' )
                return C-'a'+10;
        if ( 'A' <= C && C <= 'F' )
                return C-'A'+10;
}

template < const_iterator_t Radix, typename Char, typename Integer >
struct radix_accumulator {
        radix_accumulator () : value_(0) {}
        void operator *= ( Char const& chr ) {
                this->value_ *= Radix;
                this->value_ += char_to_int(chr);
        }
        void operator *= ( radix_accumulator const& bacc ) {
                this->value_ *= Radix;
                this->value_ += bacc.value_;
        }
        Integer value_;
};

template < const_iterator_t Radix=10,
                        typename Char=char,
                        typename Integer=unsigned long >
                        struct natural;

template < typename rtype,
                        typename Char,
                        typename Integer,
                        Integer Radix,
                        typename parser >
struct result_type<rtype,natural<Radix,Char,Integer>,parser> {
        typedef Integer type;
};

template < typename Char, typename Integer >
struct natural<2,Char,Integer> {
        template < typename rtype, typename iterator, typename parser >
        result<Integer> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef radix_accumulator<2,Char,Integer> binacc;
                if ( end == itr )
                        return result<Integer>(false,0);
                iterator tmp = itr;
                result<binacc> R = ( +alpha<numeral<2>,Char>() ).parse<binacc>(tmp,end);
                if ( not R.truth )
                        return result<Integer>(false,0);
                itr = tmp;
                return result<Integer>(true,R.value.value_);
        }
        template < typename rtype, typename iterator >
        result<Integer> parse ( iterator& itr, iterator const& end ) {
                return this->parse<rtype>(itr,end,*this);
        }
};

template < typename Char, typename Integer >
struct natural<8,Char,Integer> {
        template < typename rtype, typename iterator, typename parser >
        result<Integer> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef radix_accumulator<8,Char,Integer> octacc;
                if ( end == itr )
                        return result<Integer>(false,0);
                iterator tmp = itr;
                result<octacc> R = ( +alpha<numeral<8>,Char>() ).parse<octacc>(tmp,end);
                if ( not R.truth )
                        return result<Integer>(false,0);
                return result<Integer>(true,R.value.value_);
        }
        template < typename rtype, typename iterator >
        result<Integer> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char, typename Integer >
struct natural<10,Char,Integer> {
        template < typename rtype, typename iterator, typename parser >
        result<Integer> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef radix_accumulator<10,Char,Integer> decacc;
                if ( end == itr )
                        return result<Integer>(false,0);
                iterator tmp = itr;
                result<decacc> R = ( +alpha<numeral<10>,Char>() ).parse<decacc>(tmp,end);
                if ( not R.truth )
                        return result<Integer>(false,0);
                return result<Integer>(true,R.value.value_);
        }
        template < typename rtype, typename iterator >
        result<Integer> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

template < typename Char, typename Integer >
struct natural<16,Char,Integer> {
        template < typename rtype, typename iterator, typename parser >
        result<Integer> parse ( iterator& itr, iterator const& end, parser& self ) {
                typedef radix_accumulator<16,Char,Integer> decacc;
                if ( end == itr )
                        return result<Integer>(false,0);
                iterator tmp = itr;
                result<decacc> R = ( +alpha<numeral<16>,Char>() ).parse<decacc>(tmp,end);
                if ( not R.truth )
                        return result<Integer>(false,0);
                return result<Integer>(true,R.value.value_);
        }
        template < typename rtype, typename iterator >
        result<Integer> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
};

}// end grammar

}// end chiltbx

#endif//CHILTBX_GRAMMAR_STANDARD_PARSERS_H
