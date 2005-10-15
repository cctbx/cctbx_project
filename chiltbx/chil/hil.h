/*
  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.
*/

#ifndef CHILTBX_CHIL_HIL_H
#define CHILTBX_CHIL_HIL_H

#include<chiltbx/handle.h>
#include<chiltbx/any.h>
#include<string>
#include<map>
#include<vector>

namespace chiltbx {

namespace chil {

namespace ns_any = chiltbx::any;

template < typename Char >
struct hil_error {
        typedef std::basic_string<Char> string;
        hil_error ( string const& str="" )
        : what_("Error chiltbx::chi::hil: "+str) {}
        char const* what () const {
                return this->what_.c_str();
        }
        string what_;
};

// heirarchical interface language object
template < typename Char >
struct hil {
        typedef hil_error<Char>                                                 hil_error_type;
        typedef hil<Char>                                                               hil_type;
        typedef std::basic_string<Char>                                 string;
        typedef std::map<string,hil_type>                               hil_map;
        typedef typename hil_map::iterator                              hil_iterator;
        typedef typename hil_map::const_iterator                const_hil_iterator;
        typedef std::map<string,string>                                 attribute_map;
        typedef typename attribute_map::iterator                attribute_iterator;
        typedef typename attribute_map::const_iterator  const_attribute_iterator;
        typedef std::vector<string>                                             string_vector;
        string name;
        string string_value;
        ns_any::any<> any_value;
        hil_map children;
        attribute_map attributes;
        string_vector includes;
        string_vector parse_command_line_arguments ( string_vector const& args ) {
                string_vector not_chill_arguments;
                for ( std::size_t i=0; i<args.size(); ++i )
                        try {
                                std::pair<string_vector,string> svs = this->split(args[i]);
                                hil &hil_match = this->match(svs.first);
                                hil_match.string_value = svs.second;
                        } catch ( hil_error_type const& ) {
                                not_chill_arguments.push_back( args[i] );
                        }
                return not_chill_arguments;
        }
        hil& match ( string_vector const& sv ) {
                hil *response = this->find_match(sv);
                if ( 0 == response )
                        throw hil_error_type("No viable match found.");
                return *response;
        }
        hil* find_match ( string_vector const& sv, std::size_t pos=std::size_t(-1) ) {
                std::vector<hil*> responses;
                if ( std::size_t(-1) == pos ) {
                        hil *self = this->find_match(sv,0);
                        if ( 0 != self )
                                responses.push_back( self );
                } else {
                        if ( string::npos == this->name.find(sv[pos]) )
                                return 0;
                        if ( sv.size() == (pos+1) )
                                return this;
                }
                for ( hil_iterator htr=this->children.begin();
                        htr!=this->children.end(); ++htr ) {
                        hil *response = htr->second.find_match(sv,(pos!=std::size_t(-1)?pos+1:pos));
                        if ( 0 != response )
                                responses.push_back( response );
                }
                if ( 1==responses.size() )
                        return responses[0];
                return 0;
        }
        std::pair<string_vector,string> split ( string const& S ) {
                std::size_t equals = S.find_first_of("=");
                if ( string::npos == equals )
                        throw hil_error_type("Invalid command-line argument.");
                string dots = S.substr(0,equals);
                std::size_t notws = dots.find_first_of("\t \n");
                if ( string::npos != notws )
                        dots = dots.substr(0,notws);
                while ( dots[0] == '-' )
                        dots.erase(dots.begin());
                string value = S.substr(equals+1,S.size());
                notws = value.find_first_not_of("\t \n");
                if ( string::npos != notws )
                        value = value.substr(notws,value.size());
                string_vector parts;
                std::size_t start = 0;
                std::size_t pos = dots.find(".",start);
                while ( string::npos != pos ) {
                        if ( pos-start>0 && start<dots.size() )
                                parts.push_back( dots.substr(start,pos-start) );
                        start = pos+1;
                        pos = dots.find(".",start);
                }
                if ( pos-start>0 && start<dots.size() )
                        parts.push_back( dots.substr(start,pos-start) );
                return std::make_pair(parts,value);
        }
        template < typename T >
        T& get () {
                return this->any_value.template get<T>();
        }
        hil& operator [] ( string const& str ) {
                hil_iterator htr = this->children.find(str);
                if ( this->children.end() == htr )
                        return *this;
                return htr->second;
        }
        hil const& operator [] ( string const& str ) const {
                const_hil_iterator htr = this->children.find(str);
                if ( this->children.end() == htr )
                        return *this;
                return htr->second;
        }
        void hil_union ( hil const& H ) {
                if ( not H.name.empty() )
                        this->name = H.name;
                if ( not H.string_value.empty() )
                        this->string_value = H.string_value;
                this->attribute_union(H.attributes);
                this->children_union(H.children);
                this->include_union(H.includes);
        }
        void include_union ( string_vector const& paths ) {
                this->includes.insert(this->includes.end(),paths.begin(),paths.end());
        }
        void children_union ( hil_map const& hils ) {
                for ( const_hil_iterator htr=hils.begin();
                        htr!=hils.end(); ++htr )
                        this->children[htr->first].hil_union(htr->second);
        }
        void attribute_union ( attribute_map const& attrs ) {
                this->attributes.insert(attrs.begin(),attrs.end());
        }
        void clear () {
                this->name = "";
                this->string_value = "";
                this->children.clear();
        }
        string print ( string const& indent="" ) const {
                string more = "  ";
                string result = indent + this->name +
                                                (this->string_value.empty()?"":" = " + this->string_value) +
                                                (this->children.empty()?"\n":" {\n");
                for ( std::size_t i=0; i<this->includes.size(); ++i )
                        result += more + indent + "include " + this->includes[i] + "\n";
                for ( const_attribute_iterator itr=this->attributes.begin();
                        itr!=this->attributes.end(); ++itr )
                        result += more + indent + "." + itr->first + " = " + itr->second + "\n";
                if ( not this->children.empty() ) {
                        for ( const_hil_iterator htr=this->children.begin();
                                htr!=this->children.end(); ++htr )
                                result += htr->second.print(more+indent);
                        result += indent + "}\n";
                }
                return result;
        }
};

}// end chil

}// end chiltbx

#endif//CHILTBX_CHIL_HIL_H
