/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_FILE_ITERATOR_H
#define CHILTBX_FILE_ITERATOR_H

#include<fstream>
#include<string>

namespace chiltbx {

namespace hil {

template < typename Char, typename PosType=std::size_t, PosType BufferSize=512 >
class input_file_iterator {
public:
        typedef std::basic_ifstream<Char>       ifstream;
        typedef std::basic_string<Char>         string;
        enum iterator_type {
                end = 0,
                beg = 1,
        };
        typedef input_file_iterator self;
        input_file_iterator ( char const* fname, iterator_type type=beg ) {
                this->stream_ = 0;
                this->set(fname,type);
        }
        input_file_iterator ( self const& S ) {
                this->stream_ = 0;
                this->clone(S);
        }
        self& operator = ( self const& S ) {
                this->clone(S);
                return *this;
        }
        friend std::size_t operator - ( self const& L, self const& R ) {
                return L.stream_->tellg() - R.stream_->tellg();
        }
        friend bool operator == ( self const& L, self const& R ) {
                if ( L.file_name_ != R.file_name_ )
                        return false;
                return L.stream_->tellg() == R.stream_->tellg();
        }
        friend bool operator != ( self const& L, self const& R ) {
                return not (L==R);
        }
        self& operator ++ () {
                if ( this->is_good() )
                        this->current_ = this->stream_->get();
                return *this;
        }
        self& operator ++ (int) {
                if ( this->is_good() )
                        return ++*this;
                return *this;
        }
        self& operator -- () {
                if ( not this->is_good() )
                        return *this;
                if ( this->stream_->tellg() > 0 ) {
                        this->stream_->seekg(-2,std::ios::cur);
                        this->current_ = this->stream_->get();
                }
                return *this;
        }
        self& operator -- (int) {
                if ( this->is_good() )
                        return --*this;
                return *this;
        }
        char operator * () {
                return this->current_;
        }
        virtual ~ input_file_iterator () {
                if ( 0 != this->stream_ )
                        delete this->stream_;
                this->stream_ = 0;
        }
        bool is_good () const {
                return 0 != this->stream_ && bool(*this->stream_);
        }
protected:
        void set ( const char *fname, iterator_type type ) {
                this->clear();
                this->file_name_ = string(fname);
                this->stream_ = new ifstream(fname);
                if ( this->is_good() )
                        switch ( type ) {
                        case beg:
                                 this->current_ = char(this->stream_->get()); break;
                        case end:
                                 this->stream_->seekg(0,std::ios::end);
                                 this->current_ = 0;
                                 break;
                        default:
                                this->clear();
                        }
                else
                        this->current_ = 0;
        }
        void clone ( self const& S ) {
                this->clear();
                this->file_name_ = S.file_name_;
                this->stream_ = new ifstream(this->file_name_.c_str());
                if ( not this->is_good() ) {
                        this->clear();
                        return;
                }
                if ( S.is_good() && this->is_good() )
                        this->stream_->seekg(S.stream_->tellg());
                else if ( std::size_t(-1) == S.stream_->tellg() )
                        this->stream_->seekg(std::ios::end);
                this->current_ = S.current_;
        }
        void clear () {
                if ( 0 != this->stream_ )
                        delete this->stream_;
                this->stream_ = 0;
                this->file_name_ = "";
                this->current_ = 0;
        }
private:
        char current_;
        string file_name_;
        ifstream *stream_;
};

template < typename Char, typename PosType >
input_file_iterator<Char,PosType> begin ( char const* fname ) {
        return input_file_iterator<Char,PosType>(fname,input_file_iterator<Char,PosType>::beg);
}

template < typename Char, typename PosType >
input_file_iterator<Char,PosType> end ( char const* fname ) {
        return input_file_iterator<Char,PosType>(fname,input_file_iterator<Char,PosType>::end);
}

}//end hil namespace

}//end chiltbx

#endif//CHILTBX_FILE_ITERATOR_H
