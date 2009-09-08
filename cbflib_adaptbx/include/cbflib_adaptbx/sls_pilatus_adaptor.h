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
    int *id;
    long start;
    int eltype_file, elsigned_file, elunsigned_file,
                   minelem_file, maxelem_file, bits, sign;
    unsigned int compression;
    std::size_t nelem_file;
    std::size_t text_dimover;
    int *realarray;
    const char **byteorder;
    size_t *dimfast;
    size_t *dimslow;
    size_t *padding;
    size_t *dimmid;
    size_t *dimover;

    /* Check the digest (this will also decode it if necessary) */
    cbf_failnez (cbf_check_digest (column, row))

    /* Is it an encoded binary section? */
    SCITBX_ASSERT(!cbf_is_mimebinary (column, row));

    /* Parse the value */
    cbf_failnez (cbf_get_bintext (column, row, NULL,
                                id, &file, &start, NULL,
                                 NULL, NULL, &bits, &sign, realarray,
                                 byteorder, &text_dimover, dimfast, dimmid,
                                 dimslow, padding,
                                 &compression))

    if (dimover) {*dimover = text_dimover;}

    /* Position the file at the start of the binary section */
    cbf_failnez (cbf_set_fileposition (file, start, SEEK_SET))

    /* Get the parameters and position the file */
    cbf_failnez (cbf_decompress_parameters (&eltype_file, NULL,
                                          &elsigned_file, &elunsigned_file,
                                          &nelem_file,
                                          &minelem_file, &maxelem_file,
                                          compression,
                                          file))
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
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(slow,fast)));
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
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(slow,fast)));
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
