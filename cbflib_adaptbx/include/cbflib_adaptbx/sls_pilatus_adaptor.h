#ifndef SLS_PILATUS_AD_H
#define SLS_PILATUS_AD_H
#include <cbflib_adaptbx/cbf_adaptor.h>
#include <cbflib_adaptbx/cbf_byte_offset_optimized.h>
#include "cbf_binary.h"
#include "cbf_compress.h"

namespace iotbx {
  namespace detectors {

class wrapper_of_byte_decompression {
  cbf_handle* cbf_h;
  size_t elsize,nelem;
  int elsign;
  cbf_file *file;

  //variables required by cbf_get_bintext
  void *file_text;
  int id_text, checked_digest_text, bits_text, sign_text, realarray_text;
  unsigned long start_text, size_text;
  unsigned long dimover_text, dimfast_text, dimmid_text, dimslow_text, padding_text;
  char digest_text [25];
  char byteorder_text [14];
  unsigned int compression_text;

 public:
  wrapper_of_byte_decompression(cbf_handle* handle,const std::size_t &nelem):
    cbf_h(handle),
    elsize(sizeof(int)),
    elsign(1),
    nelem(nelem){
    SCITBX_ASSERT(elsize==4);
    SCITBX_ASSERT(cbf_h != NULL);
    SCITBX_ASSERT( CHAR_BIT == 8 );
    int data_bits = elsize * CHAR_BIT;
    SCITBX_ASSERT( data_bits == 32 );
    int numints = (data_bits + CHAR_BIT*sizeof (int)-1)/(CHAR_BIT*sizeof (int));
    SCITBX_ASSERT( numints == 1 );
  }

  void set_file_position(){
    /* adapted from the cbf_get_binary() function */
    cbf_node *column = (*cbf_h)->node;
    unsigned int row = (*cbf_h)->row;

    int eltype_file, elsigned_file, elunsigned_file,
                   minelem_file, maxelem_file;

    std::size_t nelem_file;

    /* Check the digest (this will also decode it if necessary) */
    cbf_failnez (cbf_check_digest (column, row))

    /* Is it an encoded binary section? */
    SCITBX_ASSERT(!cbf_is_mimebinary (column, row));

    /* Parse the value */
    get_bintext(column, row);
    //show_bintext();

    /* Position the file at the start of the binary section */
    cbf_failnez (cbf_set_fileposition (file, start_text, SEEK_SET))

    /* Get the parameters and position the file */
    cbf_failnez (cbf_decompress_parameters (&eltype_file, NULL,
                                          &elsigned_file, &elunsigned_file,
                                          &nelem_file,
                                          &minelem_file, &maxelem_file,
                                          compression_text,
                                          file))
  }

  void get_bintext(cbf_node* & column, unsigned int & row){
    SCITBX_ASSERT(cbf_is_binary (column, row));
    const char *text;
    /* Get the value */
    cbf_get_columnrow (&text, column, row);

    sscanf (text + 1, " %x %p %lx %lx %d %24s %x %d %d %14s %lu %lu %lu %lu %lu %u",
                      (unsigned int *)&id_text,
                      &file_text,
                      (unsigned long *)&start_text,
                      (unsigned long *)&size_text,
                      &checked_digest_text,
                       digest_text,
                      (unsigned int *)&bits_text,
                      &sign_text,
                      &realarray_text,
                      byteorder_text,
                      (unsigned long *)&dimover_text,
                      (unsigned long *)&dimfast_text,
                      (unsigned long *)&dimmid_text,
                      (unsigned long *)&dimslow_text,
                      (unsigned long *)&padding_text,
                      &compression_text);
    file = (cbf_file *)file_text;

  }

  void show_bintext(){
    SCITBX_EXAMINE(id_text);
    SCITBX_EXAMINE(checked_digest_text);
    SCITBX_EXAMINE(bits_text);
    SCITBX_EXAMINE(sign_text);
    SCITBX_EXAMINE(realarray_text);
    SCITBX_EXAMINE(start_text);
    SCITBX_EXAMINE(size_text);
    SCITBX_EXAMINE(dimover_text);
    SCITBX_EXAMINE(dimfast_text);
    SCITBX_EXAMINE(dimmid_text);
    SCITBX_EXAMINE(dimslow_text);
    SCITBX_EXAMINE(padding_text);
    SCITBX_EXAMINE(digest_text);
    SCITBX_EXAMINE(byteorder_text);
    SCITBX_EXAMINE(compression_text);
  }

  void decompress_byte_offset_optimized(void *value){
    size_t nelem_read;
    /* Read the binary data */
    SCITBX_ASSERT(!file->temporary);
    cbf_decompress_byte_offset_optimized(
             value,
             elsize,
             elsign,  //always 1
             nelem,
             &nelem_read,
             elsize * CHAR_BIT,
             1,
             file,
             0
    );
    SCITBX_ASSERT( nelem == nelem_read );
  }

};

class MiniCBFAdaptor: public CBFAdaptor {
 public:
  inline MiniCBFAdaptor(const std::string& filename):
    CBFAdaptor(filename){/* Create the cbf */}

  inline scitbx::af::flex_int read_data(const int& slow, const int& fast){
    private_file = std::fopen(filename.c_str(),"rb");
    if (!private_file) throw Error("minicbf file BAD_OPEN");
    cbf_failnez (cbf_read_widefile (cbf_h, private_file, MSG_DIGEST))

    /* The Python wrapper takes care of the SLS header read
    char* header_info = NULL;
    if ( !cbf_find_tag(cbf_h,"_array_data.header_contents")) {
      cbf_failnez(cbf_get_value(cbf_h,(const char * *)&header_info))
      cbf_parse_sls_header(cbf_h, header_info, 0);}  */

    /* Find the binary data */
    cbf_failnez (cbf_find_tag (cbf_h, "_array_data.data"))
    cbf_failnez (cbf_rewind_row (cbf_h))

    unsigned int compression;
    int binary_id,elsigned,elunsigned,minelement,maxelement;
    size_t elsize,elements,dim1,dim2,dim3,padding;
    char *byteorder ="little_endian";

    cbf_failnez (cbf_get_integerarrayparameters_wdims (
     cbf_h, &compression, &binary_id, &elsize, &elsigned, &elunsigned,
     &elements, &minelement, &maxelement,(const char **) &byteorder,
     &dim1, &dim2, &dim3, &padding))

    SCITBX_ASSERT(elements == slow*fast);
    SCITBX_ASSERT(elsize == sizeof(int));
    SCITBX_ASSERT(elsigned==1);

    //C++ weirdness
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(slow,fast)),scitbx::af::init_functor_null<int>());
    int* begin = z.begin();
    std::size_t sz = z.size();

    /* Read the binary data */
    cbf_failnez (cbf_get_integerarray (cbf_h, //cbf handle
                 &id,                         //ptr to binary section identifier
                 begin,                       //array ptr
                 sizeof (int),                //element size
                 1,                           //flag of signed data type
                 sz,                          //elements requested
                 &nelem_read                  //elements read
                 ))
    SCITBX_ASSERT(sz==nelem_read);

    return z;
  }

  inline scitbx::af::flex_int optimized_read_data(const int& slow, const int& fast){
    private_file = std::fopen(filename.c_str(),"rb");
    if (!private_file) throw Error("minicbf file BAD_OPEN");
    cbf_failnez (cbf_read_widefile (cbf_h, private_file, MSG_DIGEST))

    /* Find the binary data */
    cbf_failnez (cbf_find_tag (cbf_h, "_array_data.data"))
    cbf_failnez (cbf_rewind_row (cbf_h))

    unsigned int compression;
    int binary_id,elsigned,elunsigned,minelement,maxelement;
    size_t elsize,elements,dim1,dim2,dim3,padding;
    char *byteorder ="little_endian";

    cbf_failnez (cbf_get_integerarrayparameters_wdims (
     cbf_h, &compression, &binary_id, &elsize, &elsigned, &elunsigned,
     &elements, &minelement, &maxelement,(const char **) &byteorder,
     &dim1, &dim2, &dim3, &padding))

    SCITBX_ASSERT(elements == slow*fast);
    SCITBX_ASSERT(elsize == sizeof(int));
    SCITBX_ASSERT(elsigned==1);

    //C++ weirdness
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(slow,fast)),scitbx::af::init_functor_null<int>());
    int* begin = z.begin();
    std::size_t sz = z.size();

    wrapper_of_byte_decompression wrap_dee(&cbf_h,sz);
    wrap_dee.set_file_position();
    wrap_dee.decompress_byte_offset_optimized((void *)begin);

    return z;
  }

  inline scitbx::af::flex_int optimized_read_data(){
    private_file = std::fopen(filename.c_str(),"rb");
    if (!private_file) throw Error("minicbf file BAD_OPEN");
    cbf_failnez (cbf_read_widefile (cbf_h, private_file, MSG_DIGEST))

    /* Find the binary data */
    cbf_failnez (cbf_find_tag (cbf_h, "_array_data.data"))
    cbf_failnez (cbf_rewind_row (cbf_h))

    unsigned int compression;
    int binary_id,elsigned,elunsigned,minelement,maxelement;
    size_t elsize,elements,fast,slow,dim3,padding;
    char *byteorder ="little_endian";

    cbf_failnez (cbf_get_integerarrayparameters_wdims (
     cbf_h, &compression, &binary_id, &elsize, &elsigned, &elunsigned,
     &elements, &minelement, &maxelement,(const char **) &byteorder,
     &fast, &slow, &dim3, &padding))

    SCITBX_ASSERT(elements == fast*slow);
    SCITBX_ASSERT(elsize == sizeof(int));
    SCITBX_ASSERT(elsigned==1);

    //C++ weirdness
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(slow,fast)),scitbx::af::init_functor_null<int>());
    int* begin = z.begin();
    std::size_t sz = z.size();

    wrapper_of_byte_decompression wrap_dee(&cbf_h,sz);
    wrap_dee.set_file_position();
    wrap_dee.decompress_byte_offset_optimized((void *)begin);

    return z;
  }

};

  }//namespace detectors
}//namespace iotbx

#endif
