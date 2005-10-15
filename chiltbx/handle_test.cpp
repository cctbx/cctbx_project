/*
  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.
*/

#include<iostream>
#include<chiltbx/handle.h>

struct base {
        virtual int foo () const = 0;
        virtual int foo () = 0;
};

struct derived : public base {
        virtual int foo () const { return 5001; }
        virtual int foo () { return 10001; }
};

typedef chiltbx::handle::handle<base> base_handle;

struct bar {
        bar () {
                this->bh_ = derived();
        }
        void constant () const {
                std::cout << this->bh_->foo() << std::endl;
        }
        void non_constant () {
                std::cout << this->bh_->foo() << std::endl;
        }
        base_handle bh_;
};

int main ( int argc, char *argv[] ) {

        bar B;
        B.constant();
        B.non_constant();

}
