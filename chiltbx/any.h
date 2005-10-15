/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_ANY_H
#define CHILTBX_ANY_H

#include<chiltbx/handle.h>
#include<chiltbx/runtimetype.h>
#include<chiltbx/str.h>

namespace chiltbx {

namespace any {

namespace ns_handle = chiltbx::handle;

typedef ns_string ns_string;

// convert this to use the cctbx::error mechanism
struct any_error {
        any_error ( ns_string const& str="" )
        : what_("Error in chiltbx::any::any: "+str) {}
        string const& what () const {
                return this->what_;
        }
        string what_;
};

// the type store keeps both the data
// value and the type of the data that is
// being stored for recovery by the Any class
// the interface is used so that the type
// information can be extracted without knowing
// the type
template < typename T > struct type_store;
template <> struct type_store<void> {
        virtual string type () const = 0;
};
template < typename T >
struct type_store : public type_store<void> {
        type_store ( T const& t ) : value(t) {}
        virtual string type () const {
                return runtimetype<T>::name(value);
        }
        T value;
};

// a simple function like "std::make_pair" to create
// a specialized store
template < typename T >
type_store<T> make_store ( T const& t ) {
        return type_store<T>(t);
};

// the any class, it is specialized over the Storage
// policy that will be used as default by the back-end
// handle which stores the any-data.
template < typename Storage=ns_handle::by_value, typename Synch=ns_handle::unsafe >
class any {
public:
        any () {}
        any ( any const& a ) {
                this->handle_ = a.handle_;
        }
        template < typename StorageT, typename SynchT >
        any ( any<StorageT,SynchT> const& a ) {
                this->handle_ = a.handle_;
        }
        template < typename T >
        any ( T const& t ) {
                this->set(t);
        }
        template < typename T, typename Box >
        any ( T const& t, Box const& B ) {
                this->set(t,B);
        }
        any& operator = ( any const& a ) {
                this->handle_ = a.handle_;
                return *this;
        }
        template < typename StorageT, typename SynchT >
        any& operator = ( any<StorageT,SynchT> const& a ) {
                this->handle_ = a.handle_;
                return *this;
        }
        ~ any () {}
        // assignment by value either tries to replace the value
        // in the any with a new type value, or attempts to modify
        // the value stored in the any
        template < typename T >
        any& operator = ( T const& t ) {
                if ( this->type() == runtimetype<T>::name(t) )
                        this->get<T>() = t;
                else
                        this->set(t);
                return *this;
        }
        // sets the type and creates a new back-end
        template < typename T >
        void set ( T const& t ) {
                this->handle_.set<Storage>(make_store(t));
        }
        template < typename StorageP, typename T >
        void set ( T const& t ) {
                this->handle_.set<StorageP>(make_store(t));
        }
        // get's the back-end but does not allow non-const side-effects
        // guaranteed COW et al. safe
        template < typename T >
        T const& const_get () const {
                if ( not this->handle_.empty() ) {
                        type_store<T> const& ts = *reinterpret_cast<type_store<T> const*>(this->handle_.const_get_raw());
                        if ( runtimetype<T>().name() == ts.type() )
                                return ts.value;
                        throw any_error(runtimetype<T>().name()+" does not match stored type "+ts.type());
                }
                throw any_error("There is no handle.");
        }
        template < typename T >
        T& get () {
                if ( not this->handle_.empty() ) {
                        type_store<T>& ts = *reinterpret_cast<type_store<T>*>(this->handle_.get_raw());
                        if ( runtimetype<T>().name() == ts.type() )
                                return ts.value;
                        throw any_error(runtimetype<T>().name()+" does not match stored type "+ts.type());
                }
                throw any_error("There is no handle.");
        }
        template < typename T >
        T const& get () const {
                return this->const_get<T>();
        }
        ns_string type () const {
                if ( not this->handle_.empty() )
                        return this->handle_.const_get_raw()->type();
                return ns_string("");
        }
        bool empty () const {
                return this->handle_.empty();
        }
private:
        ns_handle::handle<type_store<void>,Storage,Synch> handle_;
};

}// end Any namespace

// runtimetype specialization for the any class
template < typename Storage, typename Synch >
struct runtimetype<any::any<Storage,Synch> > {
        static string name () {
                return string("any::any");
        }
        static string name ( any::any<Storage,Synch> const& value ) {
                return value.type();
        }
};

}// end chiltbx

#endif//CHILTBX_ANY_H
