/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_HANDLE_H
#define CHILTBX_HANDLE_H

#include<chiltbx/childef.h>
#include<chiltbx/str.h>

namespace chiltbx {

namespace handle {

struct handle_error {
    handle_error ( string const& str=string() ) : what_("Error in chiltbx::handle::handle: "+str) {}
        const string what () const {
                return this->what_;
        }
        string what_;
};

// defines an arbitrary type for the default storage policy
#ifndef CHILTBX_GENERIC_HANDLE_DEFAULT_STORAGE_POLICY_TYPE
    #define CHILTBX_GENERIC_HANDLE_DEFAULT_STORAGE_POLICY_TYPE chiltbx::handle::by_value
#endif

//
// ---- Policies ----
//
//      interface: defines classes which templatize to be interfaces
//
//      by_value: holds the value, T
//      on_heap: holds the value, T* = new T
//      counted: reference counted
//      copy_on_write: derives from counted, copies when "written" to
//
//      unsafe: the standard synchronization method
//              NOTE: only synchronizes reference-counting related mechanisms!!!
//              My personal experience is limited, this may need to be rewritten!!!
//              IMPORTANT UPDATE: Ralph! Please help the synchronization!
//
// ---- type_* ----
//
//      type_interface: carries a type and an interface to that type (a "base class")
//
//      type_store: carries the back-end storage type and the front-end storage type
//                      allowing differently-typed back-ends to be friended to their front-end
//                      basically ... cheap "friend injection"
//
struct interface;
struct by_value;
struct on_heap;
struct counted;
struct copy_on_write;
struct unsafe {};
template < typename T, typename Interface > struct type_interface;
template < typename BackStore, typename FrontStore > struct type_store;

// a simple class for counting and storage
template < typename T >
struct counter {
        counter () : count(0) {}
        counter ( T const& t ) : count(1), value(t) {}
        std::size_t count;
        T value;
};

// the handle type with defaults, ole!
template < typename Interface=void,
           typename Storage=CHILTBX_GENERIC_HANDLE_DEFAULT_STORAGE_POLICY_TYPE,
           typename Synchronization=unsafe > class handle;

// the interface to the back-ends
template < typename Interface, typename Synch >
class handle<type_interface<interface,Interface>,void,Synch> {
public:
        typedef handle<type_interface<interface,Interface>,void,Synch> h_interface;
        typedef Interface* interface_type;
        typedef Interface const* const_interface_type;

        virtual ~ handle () {}
        virtual h_interface * clone () const = 0;
        virtual interface_type get () = 0;
        virtual const_interface_type const_get () const = 0;
};

//
// A set of back-ends specialized over their storage mechanisms
//

//
// by_value
//
template < typename T, typename Interface, typename Storage, typename Synch >
class handle<type_interface<T,Interface>,type_store<by_value,Storage>,Synch>
: protected handle<type_interface<interface,Interface>,void,Synch> {
protected:
        typedef handle<type_interface<interface,Interface>,void,Synch> h_interface;
        typedef typename h_interface::interface_type interface_type;
        typedef typename h_interface::const_interface_type const_interface_type;
        friend class handle<Interface,Storage>;

        handle () {}
        handle ( T const& t ) : value_(t) {}
        virtual ~ handle () {}
        virtual h_interface * clone () const {
                return new handle(this->value_);
        }
        virtual interface_type get () {
                return &this->value_;
        }
        virtual const_interface_type const_get () const {
                return &this->value_;
        }
        T value_;
};

//
// on_heap
// this is here for completeness' sake; I'm too lazy to figure out how to get it
// to work without breaking the API.
//
template < typename T, typename Interface, typename Storage, typename Synch >
class handle<type_interface<T,Interface>,type_store<on_heap,Storage>,Synch>
: protected handle<type_interface<interface,Interface>,void,Synch> {
protected:
        typedef handle<type_interface<interface,Interface>,void,Synch> h_interface;
        typedef typename h_interface::interface_type interface_type;
        typedef typename h_interface::const_interface_type const_interface_type;
        friend class handle<Interface,Storage>;

        handle () : value_(0) {}
        handle ( T const& t ) {
                this->value_ = new T[1];
                this->value_[0] = t;
        }
        virtual ~ handle () {
                this->clear();
        }
        virtual h_interface * clone () const {
                return new handle(*this->value_);
        }
        virtual interface_type get () {
                return this->value_;
        }
        virtual const_interface_type const_get () const {
                return this->value_;
        }
        void clear () {
                if ( 0 != this->value_ )
                        delete [] this->value_;
                this->value_ = 0;
        }
        T *value_;
};

//
// reference counted back-end
// this needs to be made synchronous-safe
//
template < typename T, typename Interface, typename Storage, typename Synch >
class handle<type_interface<T,Interface>,type_store<counted,Storage>,Synch>
: protected handle<type_interface<interface,Interface>,void,Synch> {
protected:
        typedef handle<type_interface<interface,Interface>,void,Synch> h_interface;
        typedef typename h_interface::interface_type interface_type;
        typedef typename h_interface::const_interface_type const_interface_type;
        friend class handle<Interface,Storage>;

        handle () : reference_(0) {}
        handle ( T const& t ) {
                this->reference_ = new counter<T>(t);
        }
        handle ( counter<T> *ref ) {
                Synch lock;
                this->reference_ = ref;
                ++this->reference_->count;
        }
        virtual ~ handle () {
                this->clear();
        }
        virtual h_interface * clone () const {
                return new handle(this->reference_);
        }
        virtual interface_type get () {
                return &this->reference_->value;
        }
        virtual const_interface_type const_get () const {
                return &this->reference_->value;
        }
        void clear () {
                Synch lock;
                if ( 0 != this->reference_ )
                        --this->reference_->count;
                if ( 0 == this->reference_->count )
                        delete this->reference_;
                this->reference_ = 0;
        }
        counter<T> *reference_;
};

//
// copy-on-write (really simple); inherits from counted, also needs to be made
// synchronous-safe.
//
template < typename T, typename Interface, typename Storage, typename Synch >
class handle<type_interface<T,Interface>,type_store<copy_on_write,Storage>,Synch>
: protected handle<type_interface<T,Interface>,type_store<counted,Storage>,Synch> {
protected:
        typedef handle<type_interface<T,Interface>,type_store<counted,Storage>,Synch> base;
        typedef typename base::h_interface h_interface;
        typedef typename base::interface_type interface_type;
        friend class handle<Interface,Storage>;

        typedef handle<type_interface<T,Interface>,type_store<copy_on_write,Storage>,Synch> self_type;

        handle () : base() {}
        handle ( T const& t ) : base(t) {}
        handle ( counter<T> *ref ) : base(ref) {}
        virtual h_interface * clone () const {
                return new self_type(this->reference_);
        }
        virtual interface_type get () {
                Synch lock;
                if ( 1 == this->reference_->count )
                        return &this->reference_->value;
                counter<T> *local_reference = new counter<T>(this->reference_->value);
                --this->reference_->count;
                this->reference_ = local_reference;
                return &this->reference_->value;
        }
};

// the end-user's handle class
//
//      Interface: the base type of the pointer which is returned
//              e.g.
//                      int -> int*
//                      my_base_class -> my_base_class*
//                      double* -> double**
//
//      Storage: the storage policy
//
//      Synch: the Synchronization (threading) policy
//
template < typename Interface, typename Storage, typename Synch >
class handle {
private:
        typedef handle<type_interface<interface,Interface>,void,Synch> h_interface;
public:
        typedef typename h_interface::interface_type interface_type;
        typedef typename h_interface::const_interface_type const_interface_type;
        typedef handle<Interface,Storage> self_type;

        handle () {
                this->handle_ = 0;
        }
        handle ( handle const& H ) {
                this->handle_ = 0;
                this->copy(H);
        }
        handle& operator = ( handle const& H ) {
                this->copy(H);
                return *this;
        }
        template < typename StorageP >
        handle ( handle<Interface,StorageP> const& H ) {
                this->handle_ = 0;
                this->copy(H);
        }
        template < typename T >
        handle ( T const& t ) {
                this->handle_ = 0;
                this->set(t);
        }
        virtual ~ handle () {
                this->clear();
        }
        template < typename T, typename StorageP >
        handle ( T const& t, StorageP ) {
                this->handle_ = 0;
                this->set<StorageP>(t);
        }
        template < typename StorageP >
        handle& operator = ( handle<Interface,StorageP> const& H ) {
                this->clear();
                this->copy(H);
                return *this;
        }
        template < typename T >
        handle& operator = ( T const& t ) {
                this->set(t);
                return *this;
        }
        template < typename T >
        void set ( T const& t ) {
                this->set<Storage>(t);
        }
        template < typename StorageP, typename T >
        void set ( T const& t ) {
                this->clear();
                this->handle_ = new handle<type_interface<T,Interface>,type_store<StorageP,Storage> >(t);
        }
        template < typename T >
        T& get () {
                if ( !this->empty() )
                        return *reinterpret_cast<T*>(this->handle_->get());
                throw handle_error("Handle is empty.");
        }
        template < typename T >
        T const& get () const {
                return this->const_get<T>();
        }
        // this is here to /guarantee/ a constant access
        template < typename T >
        T const& const_get () const {
                if ( !this->empty() )
                        return *reinterpret_cast<T const*>(this->handle_->const_get());
                throw handle_error("Handle is empty.");
        }
        interface_type get_raw () {
                if ( !this->empty() )
                        return this->handle_->get();
                return 0;
        }
        const_interface_type get_raw () const {
                if ( !this->empty() )
                        return this->handle_->const_get();
                return 0;
        }
        // this is here to /guarantee/ a constant access
        const_interface_type const_get_raw () const {
                if ( !this->empty() )
                        return this->handle_->const_get();
                return 0;
        }
        // pointer-like interface
        interface_type operator -> () {
                return this->get_raw();
        }
        const_interface_type operator -> () const {
                return this->const_get_raw();
        }
        void clear () {
                if ( !this->empty() )
                        delete this->handle_;
                this->handle_ = 0;
        }
        template < typename StorageP >
        void copy ( handle<Interface,StorageP> const& H ) {
                this->clear();
                if ( !H.empty() )
                        this->handle_ = H.handle_->clone();
        }
        bool empty () const {
                return 0 == this->handle_;
        }
private:
        h_interface * handle_;
};

}// end handle namespace

}// end chiltbx namespace

#endif//CHILTBX_HANDLE_H
