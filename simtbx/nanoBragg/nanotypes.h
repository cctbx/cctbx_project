/*
 * nanotypes.h
 *
 *  Created on: Jan 2, 2018
 *      Author: giles
 */

#ifndef NANOTYPES_H_
#define NANOTYPES_H_

typedef enum { SAMPLE, BEAM } pivot;
typedef enum { UNKNOWN, SQUARE, ROUND, GAUSS, GAUSS_ARGCHK, TOPHAT, FIBER } shapetype;
// GAUSS_ARGCHK provides a lightweight backdoor for efficient implementation of
// the GAUSS shapetype in CUDA. Precaluates the argument of the exponential.
typedef enum { CUSTOM, ADXV, MOSFLM, XDS, DIALS, DENZO } convention;

#endif /* NANOTYPES_H_ */
