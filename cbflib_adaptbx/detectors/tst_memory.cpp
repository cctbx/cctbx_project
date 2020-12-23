#include <cstdio>
#include <string>
#include <exception>
#include <iostream>
#include <include/cbf.h>
#include <include/cbf_simple.h>

struct Error : public std::exception {
  std::string s;
  Error(std::string s):s(s){}
  virtual const char* what() const throw() {return s.c_str();}
  virtual ~Error() throw() {}
};

int main() {
  std::string file("adscconverted_orig.cbf");

  for (int cc=0; cc<20000; ++cc) {
    printf("Iteration %8d\n",cc);
    cbf_handle cbf_h;
    FILE* private_file = std::fopen(file.c_str(),"rb");
    if (!private_file) throw Error("cbf file BAD_OPEN");

    cbf_failnez ( cbf_make_handle (&cbf_h) )
    cbf_failnez ( cbf_read_file (cbf_h, private_file, MSG_DIGEST))
      //file handle must be left open & is closed by the cbf library.

    cbf_detector detector1;
    cbf_failnez ( cbf_construct_detector(cbf_h,&detector1,0) )
    cbf_failnez ( cbf_free_detector(detector1) )
    cbf_failnez ( cbf_free_handle (cbf_h))
    //fclose(private_file);
  }
}
