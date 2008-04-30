/* If you change the include guards, please be sure to also rename the
   functions below. Otherwise your project will clash with the original
   iotbx declarations and definitions.
 */
#ifndef IOTBX_PDB_HYBRID_36_C_H
#define IOTBX_PDB_HYBRID_36_C_H

#ifdef __cplusplus
extern "C" {
#endif

#define HY36_WIDTH_4_MIN -999
#define HY36_WIDTH_4_MAX 2436111 /* 10000 + 2*26*36*36*36 - 1 */
#define HY36_WIDTH_5_MIN -9999
#define HY36_WIDTH_5_MAX 87440031 /* 100000 + 2*26*36*36*36*36 - 1 */

const char*
hy36encode(unsigned width, int value, char* result);

const char*
hy36decode(unsigned width, const char* s, unsigned s_size, int* result);

#ifdef __cplusplus
}
#endif
#endif /* IOTBX_PDB_HYBRID_36_C_H */
