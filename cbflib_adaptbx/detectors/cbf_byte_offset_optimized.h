#ifndef CBF_BYTE_OFFSET_OPTIMIZED_H
#define CBF_BYTE_OFFSET_OPTIMIZED_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <include/cbf_file.h>

  /* Decompress an array with the byte-offset algorithm */

int cbf_decompress_byte_offset_optimized
                               (void         *destination,
                                size_t        elsize,
                                int           elsign,
                                size_t        nelem,
                                size_t       *nelem_read,
                                int           bits,
                                int           sign,
                                cbf_file     *file,
                                int           realarray
                                );

#ifdef __cplusplus
}
#endif
#endif /* CBF_BYTE_OFFSET_OPTIMIZED_H */
