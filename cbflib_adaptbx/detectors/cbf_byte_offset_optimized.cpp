/**********************************************************************
 * An optimized byte-offset decompression                             *
 *                                                                    *
 * Based on CBFlib Version 0.8.0 20 July 2008                         *
 *                                                                    *
 * Assumptions:                                                       *
 *  elsize:  the size of the integer element is 4 bytes               *
 *  elsign:  is always 1.  Means integer element type is signed?      *
 *  data_bits: the number of bits in the integer element is 32        *
 *  CHAR_BIT: the number of bits in a byte is 8                       *
 *  data_sign:  is always 1.  Signed integer type for the data?       *
 *  realarray: is 0; the type of the data is always integer, not real *
 *                                                                    *
 * (C) Copyright 2009 Nicholas K. Sauter                              *
 *                                                                    *
 **********************************************************************/
/**********************************************************************
 * cbf_byte_offset -- byte-offset compression                         *
 *                                                                    *
 * Version 0.8.0 20 July 2008                                         *
 *                                                                    *
 *                          Paul Ellis and                            *
 *         Herbert J. Bernstein (yaya@bernstein-plus-sons.com)        *
 *                                                                    *
 * (C) Copyright 2006, 2007 Herbert J. Bernstein                      *
 *                                                                    *
 **********************************************************************/

/************************* LGPL NOTICES *******************************
 *                                                                    *
 * This library is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU Lesser General Public         *
 * License as published by the Free Software Foundation; either       *
 * version 2.1 of the License, or (at your option) any later version. *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
 * Lesser General Public License for more details.                    *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License along with this library; if not, write to the Free         *
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
 * MA  02110-1301  USA                                                *
 *                                                                    *
 **********************************************************************/

/**********************************************************************
 *                                                                    *
 *                    Stanford University Notices                     *
 *  for the CBFlib software package that incorporates SLAC software   *
 *                 on which copyright is disclaimed                   *
 *                                                                    *
 * This software                                                      *
 * -------------                                                      *
 * The term ‘this software’, as used in these Notices, refers to      *
 * those portions of the software package CBFlib that were created by *
 * employees of the Stanford Linear Accelerator Center, Stanford      *
 * University.                                                        *
 *                                                                    *
 * Stanford disclaimer of copyright                                   *
 * --------------------------------                                   *
 * Stanford University, owner of the copyright, hereby disclaims its  *
 * copyright and all other rights in this software.  Hence, anyone    *
 * may freely use it for any purpose without restriction.             *
 *                                                                    *
 * Acknowledgement of sponsorship                                     *
 * ------------------------------                                     *
 * This software was produced by the Stanford Linear Accelerator      *
 * Center, Stanford University, under Contract DE-AC03-76SFO0515 with *
 * the Department of Energy.                                          *
 *                                                                    *
 * Government disclaimer of liability                                 *
 * ----------------------------------                                 *
 * Neither the United States nor the United States Department of      *
 * Energy, nor any of their employees, makes any warranty, express or *
 * implied, or assumes any legal liability or responsibility for the  *
 * accuracy, completeness, or usefulness of any data, apparatus,      *
 * product, or process disclosed, or represents that its use would    *
 * not infringe privately owned rights.                               *
 *                                                                    *
 * Stanford disclaimer of liability                                   *
 * --------------------------------                                   *
 * Stanford University makes no representations or warranties,        *
 * express or implied, nor assumes any liability for the use of this  *
 * software.                                                          *
 *                                                                    *
 * Maintenance of notices                                             *
 * ----------------------                                             *
 * In the interest of clarity regarding the origin and status of this *
 * software, this and all the preceding Stanford University notices   *
 * are to remain affixed to any copy or derivative of this software   *
 * made or distributed by the recipient and are to be affixed to any  *
 * copy of software made or distributed by the recipient that         *
 * contains a copy or derivative of this software.                    *
 *                                                                    *
 * Based on SLAC Software Notices, Set 4                              *
 * OTT.002a, 2004 FEB 03                                              *
 **********************************************************************/
/**********************************************************************
 *                               NOTICE                               *
 * Creative endeavors depend on the lively exchange of ideas. There   *
 * are laws and customs which establish rights and responsibilities   *
 * for authors and the users of what authors create.  This notice     *
 * is not intended to prevent you from using the software and         *
 * documents in this package, but to ensure that there are no         *
 * misunderstandings about terms and conditions of such use.          *
 *                                                                    *
 * Please read the following notice carefully.  If you do not         *
 * understand any portion of this notice, please seek appropriate     *
 * professional legal advice before making use of the software and    *
 * documents included in this software package.  In addition to       *
 * whatever other steps you may be obliged to take to respect the     *
 * intellectual property rights of the various parties involved, if   *
 * you do make use of the software and documents in this package,     *
 * please give credit where credit is due by citing this package,     *
 * its authors and the URL or other source from which you obtained    *
 * it, or equivalent primary references in the literature with the    *
 * same authors.                                                      *
 *                                                                    *
 * Some of the software and documents included within this software   *
 * package are the intellectual property of various parties, and      *
 * placement in this package does not in any way imply that any       *
 * such rights have in any way been waived or diminished.             *
 *                                                                    *
 * With respect to any software or documents for which a copyright    *
 * exists, ALL RIGHTS ARE RESERVED TO THE OWNERS OF SUCH COPYRIGHT.   *
 *                                                                    *
 * Even though the authors of the various documents and software      *
 * found here have made a good faith effort to ensure that the        *
 * documents are correct and that the software performs according     *
 * to its documentation, and we would greatly appreciate hearing of   *
 * any problems you may encounter, the programs and documents any     *
 * files created by the programs are provided **AS IS** without any   *
 * warranty as to correctness, merchantability or fitness for any     *
 * particular or general use.                                         *
 *                                                                    *
 * THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF    *
 * PROGRAMS OR DOCUMENTS OR ANY FILE OR FILES CREATED BY USE OF THE   *
 * PROGRAMS OR DOCUMENTS LIES SOLELY WITH THE USERS OF THE PROGRAMS   *
 * OR DOCUMENTS OR FILE OR FILES AND NOT WITH AUTHORS OF THE          *
 * PROGRAMS OR DOCUMENTS.                                             *
 **********************************************************************/

/**********************************************************************
 *                                                                    *
 *                           The IUCr Policy                          *
 *      for the Protection and the Promotion of the STAR File and     *
 *     CIF Standards for Exchanging and Archiving Electronic Data     *
 *                                                                    *
 * Overview                                                           *
 *                                                                    *
 * The Crystallographic Information File (CIF)[1] is a standard for   *
 * information interchange promulgated by the International Union of  *
 * Crystallography (IUCr). CIF (Hall, Allen & Brown, 1991) is the     *
 * recommended method for submitting publications to Acta             *
 * Crystallographica Section C and reports of crystal structure       *
 * determinations to other sections of Acta Crystallographica         *
 * and many other journals. The syntax of a CIF is a subset of the    *
 * more general STAR File[2] format. The CIF and STAR File approaches *
 * are used increasingly in the structural sciences for data exchange *
 * and archiving, and are having a significant influence on these     *
 * activities in other fields.                                        *
 *                                                                    *
 * Statement of intent                                                *
 *                                                                    *
 * The IUCr's interest in the STAR File is as a general data          *
 * interchange standard for science, and its interest in the CIF,     *
 * a conformant derivative of the STAR File, is as a concise data     *
 * exchange and archival standard for crystallography and structural  *
 * science.                                                           *
 *                                                                    *
 * Protection of the standards                                        *
 *                                                                    *
 * To protect the STAR File and the CIF as standards for              *
 * interchanging and archiving electronic data, the IUCr, on behalf   *
 * of the scientific community,                                       *
 *                                                                    *
 * * holds the copyrights on the standards themselves,                *
 *                                                                    *
 * * owns the associated trademarks and service marks, and            *
 *                                                                    *
 * * holds a patent on the STAR File.                                 *
 *                                                                    *
 * These intellectual property rights relate solely to the            *
 * interchange formats, not to the data contained therein, nor to     *
 * the software used in the generation, access or manipulation of     *
 * the data.                                                          *
 *                                                                    *
 * Promotion of the standards                                         *
 *                                                                    *
 * The sole requirement that the IUCr, in its protective role,        *
 * imposes on software purporting to process STAR File or CIF data    *
 * is that the following conditions be met prior to sale or           *
 * distribution.                                                      *
 *                                                                    *
 * * Software claiming to read files written to either the STAR       *
 * File or the CIF standard must be able to extract the pertinent     *
 * data from a file conformant to the STAR File syntax, or the CIF    *
 * syntax, respectively.                                              *
 *                                                                    *
 * * Software claiming to write files in either the STAR File, or     *
 * the CIF, standard must produce files that are conformant to the    *
 * STAR File syntax, or the CIF syntax, respectively.                 *
 *                                                                    *
 * * Software claiming to read definitions from a specific data       *
 * dictionary approved by the IUCr must be able to extract any        *
 * pertinent definition which is conformant to the dictionary         *
 * definition language (DDL)[3] associated with that dictionary.      *
 *                                                                    *
 * The IUCr, through its Committee on CIF Standards, will assist      *
 * any developer to verify that software meets these conformance      *
 * conditions.                                                        *
 *                                                                    *
 * Glossary of terms                                                  *
 *                                                                    *
 * [1] CIF:  is a data file conformant to the file syntax defined     *
 * at http://www.iucr.org/iucr-top/cif/spec/index.html                *
 *                                                                    *
 * [2] STAR File:  is a data file conformant to the file syntax       *
 * defined at http://www.iucr.org/iucr-top/cif/spec/star/index.html   *
 *                                                                    *
 * [3] DDL:  is a language used in a data dictionary to define data   *
 * items in terms of "attributes". Dictionaries currently approved    *
 * by the IUCr, and the DDL versions used to construct these          *
 * dictionaries, are listed at                                        *
 * http://www.iucr.org/iucr-top/cif/spec/ddl/index.html               *
 *                                                                    *
 * Last modified: 30 September 2000                                   *
 *                                                                    *
 * IUCr Policy Copyright (C) 2000 International Union of              *
 * Crystallography                                                    *
 **********************************************************************/

#ifdef __cplusplus

extern "C" {

#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "cbf.h"
#include "cbf_file.h"
#include "cbf_byte_offset.h"
#include <cbflib_adaptbx/cbf_byte_offset_optimized.h>

  /* Decompress an array with the byte-offset algorithm */

int cbf_decompress_byte_offset_optimized (void         *destination,
                                size_t        elsize,
                                int           elsign,
                                size_t        nelem,
                                size_t       *nelem_read,
                                int           data_bits,
                                int           data_sign,
                                cbf_file     *file,
                                int           realarray)
{
  unsigned int element[4], prevelement[4], sign, unsign, limit;

  unsigned int data_unsign;

  unsigned char *unsigned_char_data;

  int errorcode, overflow, numints, iint, carry;

  int delta[4];

  char * border;

  char * rformat;

  size_t numread;

    /* prepare the errorcode */

  errorcode = 0;

    /* Is the element size valid? */

  if (elsize != sizeof (int) &&
      elsize != 2* sizeof (int) &&
      elsize != 4* sizeof (int) &&
      elsize != sizeof (short) &&
      elsize != sizeof (char))

    return CBF_ARGUMENT;

    /* check for compatible real format */

  if ( realarray ) {

    cbf_failnez (cbf_get_local_real_format(&rformat) )

    if ( strncmp(rformat,"ieee",4) ) return CBF_ARGUMENT;

  }

    /* Check the stored element size */

  if (data_bits < 1 || data_bits > 64)

    return CBF_ARGUMENT;

  numints = (data_bits + CHAR_BIT*sizeof (int) -1)/(CHAR_BIT*sizeof (int));


    /* Initialise the pointer */

  unsigned_char_data = (unsigned char *) destination;


    /* Maximum limits */

  sign = 1 << ((elsize-(numints-1)*sizeof(int))* CHAR_BIT - 1);

  if (elsize == sizeof (int) || elsize == numints*sizeof(int))

    limit = ~0;

  else

    if (numints == 1 ) {

      limit = ~(-(1 << (elsize * CHAR_BIT)));

    } else {

      limit = ~(-(1 << ((elsize-(numints-1)*sizeof(int)) * CHAR_BIT)));

    }


    /* Offsets to make the value unsigned */

  if (data_sign)

    data_unsign = sign;

  else

    data_unsign = 0;

  if (elsign)

    unsign = sign;

  else

    unsign = 0;

    /* Get the local byte order */

  if (realarray) {

    cbf_get_local_real_byte_order(&border);

  } else {

    cbf_get_local_integer_byte_order(&border);

  }


    /* Set up the previous element for increments */

  prevelement[0] = prevelement[1] = prevelement[2] = prevelement[3] = 0;

  prevelement[numints-1] = data_unsign;


    /* Read the elements */

  overflow = 0;

  numread = 0;

  while (numread < nelem)
  {

    for (iint=0; iint < numints; iint++){

      element[iint] = prevelement[iint];

      delta[iint] = 0;

    }

    carry = 0;

    cbf_failnez(cbf_get_bits(file,delta,8))

    if ((delta[0]&0xFF) == 0x80) {


      cbf_failnez(cbf_get_bits(file,delta,16))


      if ( (delta[0]& 0xFFFF) == 0x8000)  {

        cbf_failnez(cbf_get_bits(file,delta,32))

        if ( (sizeof(int)==2 && delta[0] == 0 && delta[1] == 0x8000)
           || (sizeof(int)> 3 && (delta[0]&0xFFFFFFFF)==0x80000000) )  {

          cbf_failnez(cbf_get_bits(file,delta,64))

        } else {

          if (sizeof(int) == 2) {

            if (delta[1] & 0x8000)  {

              for (iint = 2; iint < numints; iint++) delta[iint] = ~0;

            }

          } else  {

            if (delta[0] & 0x80000000) {

              delta[0] |= ~0xFFFFFFFF;

              for (iint = 1; iint < numints; iint++) {

                delta[iint] = ~0;

              }

            }

          }

        }


      }  else  {

        if (delta[0] & 0x80000) {

          delta[0] |= ~0xFFFF;

          for (iint = 1; iint < numints; iint++) {

            delta[iint] = ~0;

          }
        }

      }


    } else {

      if (delta[0]&0x80)
      {
        delta[0] |= ~0xFF;

        for (iint = 1; iint < numints; iint++) {

          delta[iint] = ~0;

        }
      }

    }


    if (numints > 1) {

      for (iint = 0; iint < numints; iint++) element[iint] = prevelement[iint];

      cbf_failnez(cbf_mpint_add_acc(element,numints, (unsigned int *)delta,numints))

    } else {

      element[0] = prevelement[0] + delta[0];

      element[0] &= limit;

    }

    for (iint = 0; iint < numints; iint++)   {

      prevelement[iint] = element[iint];

    }


      /* Make the element signed? */

    element[numints-1] -= unsign;


      /* Save the element */

    if (numints > 1) {

      if (border[0] == 'b') {

        for (iint = numints; iint; iint--) {

            *((unsigned int *) unsigned_char_data) = element[iint-1];

            unsigned_char_data += sizeof (int);

        }

      } else {

        for (iint = 0; iint < numints; iint++) {

            *((unsigned int *) unsigned_char_data) = element[iint];

            unsigned_char_data += sizeof (int);
        }
      }

    } else {

      if (elsize == sizeof (int))

        *((unsigned int *) unsigned_char_data) = element[0];

      else

        if (elsize == sizeof (short))

          *((unsigned short *) unsigned_char_data) = element[0];

        else

          *unsigned_char_data = element[0];

      unsigned_char_data += elsize;

    }

    numread++;
  }

    /* Number read */

  if (nelem_read)

    *nelem_read = numread;


    /* Success */

  return overflow;
}


#ifdef __cplusplus

}

#endif
