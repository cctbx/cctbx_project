#ifndef SLS_PILATUS_AD_H
#define SLS_PILATUS_AD_H
#include <cbflib_adaptbx/cbf_adaptor.h>

namespace iotbx {
  namespace detectors {

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

};

  }//namespace detectors
}//namespace iotbx

#endif
