#ifndef GENERAL_CBF_WRITE_H
#define GENERAL_CBF_WRITE_H
#include <cbflib_adaptbx/cbf_adaptor.h>
#include "cbf_binary.h"
#include "cbf_compress.h"

namespace iotbx {
  namespace detectors {

class CBFWriteAdaptor: public CBFAdaptor {
 public:
  inline CBFWriteAdaptor(const std::string& filename):
    CBFAdaptor(filename){/* Create the cbf */}

  inline void write_data(scitbx::af::flex_int data){
    private_file = std::fopen(filename.c_str(),"wb");
    if (!private_file) throw Error("minicbf file BAD_OPEN");

    //C++ weirdness
    int* begin = data.begin();

    /* Make a new data block */

    cbf_failnez (cbf_new_datablock (cbf_h, "image_1"))

    /* Make the _array_data category */

    cbf_failnez (cbf_new_category     (cbf_h, "array_data"))
    cbf_failnez (cbf_new_column       (cbf_h, "array_id"))
    cbf_failnez (cbf_set_value        (cbf_h, "image_1"))
    cbf_failnez (cbf_new_column       (cbf_h, "binary_id"))
    cbf_failnez (cbf_set_integervalue (cbf_h, 1))
    cbf_failnez (cbf_new_column       (cbf_h, "data"))

    /* Save the binary data */

    cbf_failnez (cbf_set_integerarray_wdims (cbf_h,
                         CBF_BYTE_OFFSET,
                         1,
                         begin,
                         sizeof (int),
                         1,
                         data.size(),
                         "little_endian",
                         data.accessor().focus()[1],  //the fast size
                         data.accessor().focus()[0],  //the slow size
                         0,
                         0))


    /* Write the new file */
    cbf_failnez (cbf_write_file (cbf_h, private_file, 1, CBF, MSG_DIGEST | MIME_HEADERS, 0))
  }
};

  }//namespace detectors
}//namespace iotbx

#endif
