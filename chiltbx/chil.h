/*
  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.
*/

#ifndef CHILTBX_CHIL_H
#define CHILTBX_CHIL_H

#include<chiltbx/chil/hil.h>
#include<chiltbx/grammar.h>
#include<chiltbx/grammar/standard_parsers.h>
#include<vector>
#include<map>

namespace chiltbx {

namespace chil {

namespace G = chiltbx::grammar;

struct dot {};
struct equal {};

template < typename Char >
struct attributes {
        typedef std::basic_string<Char>                         string_type;
        typedef std::map<string_type,string_type>       attribute_map;
        typedef G::container_tag<string_type>           c_str_type;
        typedef G::terminal<c_str_type>                         strparser;
        typedef std::vector<strparser>                          strparser_vector;
        typedef std::vector<string_type>                        str_vector;
        attributes ( str_vector const& strv = str_vector() ) {
                str_vector std = attributes<Char>::standard_attributes();
                for ( std::size_t i=0; i<std.size(); ++i )
                        this->attributes_.push_back( strparser(std[i]) );
                for ( std::size_t i=0; i<strv.size(); ++i )
                        this->attributes_.push_back( strparser(strv[i]) );
        }
        static str_vector standard_attributes () {
                str_vector result;
                result.push_back( "help" );
                result.push_back( "caption" );
                result.push_back( "short_caption" );
                result.push_back( "optional" );
                result.push_back( "type" );
                result.push_back( "multiple" );
                result.push_back( "input_size" );
                result.push_back( "expert_level" );
                return result;
        }
        template < typename rtype, typename iterator, typename parser >
        G::result<string_type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                if ( end == itr )
                        return G::result<attribute_map>(false);
                for ( std::size_t i=0; i<this->attributes_.size(); ++i ) {
                        iterator tmp = itr;
                        G::result<string_type> r = this->attributes_.parse<rtype>(tmp,end,self);
                        if ( not r.truth )
                                continue;
                        itr = tmp;
                        return r;
                }
                return G::result<string_type>(false);
        }
        template < typename rtype, typename iterator >
        attribute_map parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
        strparser_vector attributes_;
};

template < typename Char >
struct dot_word {
        typedef std::basic_string<Char> string_type;
        void operator *= ( dot_word const& dw ) {
                this->value_ += ":" + dw.value_;
        }
        void operator >>= ( std::pair<dot,string_type> const& word ) {
                this->value_ = word.second;
        }
        string_type value_;
};

template < typename Char >
struct include {
        include ( bool expand=false ) : expand_(expand) {}
        typedef include<Char> self_type;
        typedef std::basic_string<Char> string_type;
        template < typename rtype, typename iterator, typename parser >
        grammar::result<self_type>
        parse ( iterator& itr, iterator const& end, parser& self ) {
                using namespace grammar;
                if ( end == itr )
                        return result<self_type>(false);
                iterator tmp = itr;
                result<string_type> r = container(string_type("include")).template parse<rtype>(tmp,end,self);
                if ( not r.truth )
                        return result<self_type>(false);
                eat<white_space<Char> >().template parse<rtype>(tmp,end,self);
                r = to_eol<Char>().template parse<rtype>(tmp,end,self);
                if ( not r.truth )
                        return result<self_type>(false);
                this->value_ = r.value;
                itr = tmp;
                return result<self_type>(true,*this);
        }
        template < typename rtype, typename iterator >
        grammar::result<self_type> parse ( iterator& itr, iterator const& ind ) {
                return this->parse<rtype>(itr,ind,*this);
        }
        bool expand_;
        string_type value_;
};

template < typename Char >
struct nspace {
        typedef std::basic_string<Char> string_type;
        typedef dot_word<Char>            dot_word_type;
        void operator >>= ( std::pair<string_type,dot_word_type> const& strdw ) {
                std::cout << strdw.first << ":" << strdw.second.value_ << std::endl;
        }
};

template < typename Char >
struct chill {
        void operator |= ( include<Char> const& inc ) {
                if ( inc.expand_ )
                        std::cout << "Expanding: " << inc.value_ << std::endl;
                else
                        std::cout << "Not Expanding: " << inc.value_ << std::endl;
        }
        void operator |= ( nspace<Char> const& ns ) {
                std::cout << "namespace" << std::endl;
        }
        void operator |= ( Char const& ns ) {
                std::cout << "namespace" << std::endl;
        }
        void operator *= ( chill<Char> const& chil ) {
                std::cout << "*= chill" << std::endl;
        }
        void operator >>= ( std::pair<grammar::null,chill<Char> > const& chil ) {
                std::cout << ">>= ws chill" << std::endl;
        }
        void operator >>= ( std::pair<chill<Char>,grammar::null> const& chil ) {
                std::cout << ">>= chill ws" << std::endl;
        }
};

template < typename Char=char >
class parser {
public:

        typedef chiltbx::chil::hil<Char>                hil_type;
        typedef chiltbx::chil::chill<Char>              chil;
        typedef chiltbx::chil::nspace<Char>     nspace;
        typedef chiltbx::chil::dot_word<Char>   dot_word;

        parser ( bool expand_includes=true, bool cool=false )
        : expand_includes_(expand_includes), cool_(cool) {}

        bool expand_includes () const {
                return this->expand_includes_;
        }
        bool& expand_includes () {
                return this->expand_includes_;
        }

        bool cool () const {
                return this->cool_;
        }
        bool& cool () {
                return this->cool_;
        }

        void chill () {
                this->cool_ = true;
        }
        void square () {
                this->cool_ = false;
        }

        template < typename iterator >
        G::result<hil_type> parse ( iterator& itr, iterator const& end ) {
                if ( end == itr )
                        return G::result<hil_type>(false);

                typedef attributes<Char>                attribute;

                namespace ozc = chiltbx::chil;
                using namespace G;

                typedef alpha_values<Char>              alpha_values;

                typename alpha_values::dot              dot;
                typename alpha_values::equal    equal;
                c_word<Char>                                    c_word;
                eat<white_space<Char> >                 eat_ws;
                chiltbx::chil::include<Char>    include(this->expand_includes_);
                attributes<Char>                                attrs;

                S<chil> ((

                        ( law<chil>()                   = +(q(eat_ws) >> (include|use<nspace>()) >> q(eat_ws)) ),
                        ( law<nspace>()                 = c_word >> use<dot_word>() ),
                        ( law<dot_word>()       = *(type<ozc::dot>(dot)>>c_word) ),
                        ( law<attribute>()      = *(q(eat_ws) >> attrs ))//,
                        //( law<assignment>()   = q(eat_ws) >> type<ozc::equal>(equal) >> q(eat_ws) >> any_but_eol_p() ),

                )).parse<chil>(itr,end);

                return G::result<hil_type>(false);
        }

private:
        bool expand_includes_;
        bool cool_;
};

}// end chil

namespace grammar {

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,chil::include<Char>,parser> {
        typedef chil::include<Char> type;
};

template < typename rtype, typename Char, typename parser >
struct result_type<rtype,chiltbx::chil::attributes<Char>,parser> {
        typedef std::basic_string<Char> type;
};

}// end grammar

}// end chiltbx

#endif//CHILTBX_CHIL_H
