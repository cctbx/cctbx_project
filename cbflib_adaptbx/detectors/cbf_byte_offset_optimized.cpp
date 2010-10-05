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
#include <cbflib_adaptbx/detectors/cbf_byte_offset_optimized.h>

/* Changes made in the byte-offset algorithm

Baseline performance for one cbf decompression                 0.26 sec
Eliminated all code that tests numints, since it is always==1  0.23 sec
Remove dependency on data_bits                                 0.23 sec
Remove dependency on realarray -- it is always 0               0.23 sec
Remove cbf_failnez macro resolution overhead                   0.23 sec
Eliminate function call overhead to cbf_get_bits(*,*,8)        0.18 sec
When bitcount==8 also omit variable count==0; while loop once  0.14 sec
Remove dependency on elsize--always 4 bytes                    0.13 sec
Eliminate cases where sizeof(int) != 4                         0.13 sec
file->characters_used == 0 for bitcount==8                     0.13 sec
double buffering--substitute fread for getc                    0.09 sec

*/

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

  int errorcode, overflow, iint, carry;

  int delta[4];

  char * border;

  char * rformat;

  size_t numread;

  /* home-brewed buffered I/O */
      unsigned char dbuf[128];
      static int dbuf_size = 128;
      size_t ifile_still_buffered;
      long cbf_file_position;

  /* ported from cbf_get_bits */
  int get_bits_bitcode, get_bits_m, get_bits_maxbits, get_bits_bitcount;
  get_bits_maxbits = sizeof (int) * CHAR_BIT; /* is always 32 */

  ifile_still_buffered = 0;


    /* prepare the errorcode */

  errorcode = 0;

    /* Initialise the pointer */

  unsigned_char_data = (unsigned char *) destination;


    /* Maximum limits */

  sign = 1 << (elsize - 1);

  limit = ~0;


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

  cbf_get_local_integer_byte_order(&border);


    /* Set up the previous element for increments */

  prevelement[0] = prevelement[1] = prevelement[2] = prevelement[3] = 0;

  prevelement[0] = data_unsign;


    /* Read the elements */

  overflow = 0;

  numread = 0;

  while (numread < nelem)
  {

    {

      element[0] = prevelement[0];

      delta[0] = 0;

    }

    carry = 0;




  get_bits_bitcount = 8;
        /* Read the bits into an int */

  get_bits_bitcode = file->bits [1] & 0x0ff;

  /*file->bits [1] = getc (file->stream);*/
  if (!ifile_still_buffered) {
    cbf_file_position = ftell(file->stream);
    ifile_still_buffered = fread(dbuf, 1, dbuf_size, file->stream);
  }

  file->bits[1] = (int)(dbuf[dbuf_size-ifile_still_buffered]);
  ifile_still_buffered--;
  cbf_file_position++;


    if (file->bits [1] == EOF)

      return CBF_FILEREAD;

    get_bits_bitcode |= (file->bits [1] << 0) & -(1 << 0);


  file->bits [1] =(file->bits [1] >> (8 + get_bits_bitcount));

  file->bits [0] = 0;


    /* Sign-extend */

  get_bits_m = 1 << (get_bits_bitcount - 1);

  if (get_bits_bitcode & get_bits_m)

    *delta = get_bits_bitcode | -get_bits_m;

  else

    *delta = get_bits_bitcode & ~-get_bits_m;






    if ((delta[0]&0xFF) == 0x80) {

      ifile_still_buffered = 0;
      fseek(file->stream, cbf_file_position, SEEK_SET);
      cbf_get_bits(file,delta,16);
      cbf_file_position = ftell(file->stream);

      if ( (delta[0]& 0xFFFF) == 0x8000)  {

        ifile_still_buffered = 0;
        fseek(file->stream, cbf_file_position, SEEK_SET);
        cbf_get_bits(file,delta,32);
        cbf_file_position = ftell(file->stream);

        if ( (delta[0]&0xFFFFFFFF)==0x80000000 )  {

          ifile_still_buffered = 0;
          fseek(file->stream, cbf_file_position, SEEK_SET);
          cbf_get_bits(file,delta,64);
          cbf_file_position = ftell(file->stream);

        } else {

          if (delta[0] & 0x80000000) {

              delta[0] |= ~0xFFFFFFFF;

          }

        }


      }  else  {

        if (delta[0] & 0x80000) {

          delta[0] |= ~0xFFFF;

        }

      }


    } else {

      if (delta[0]&0x80)
      {
        delta[0] |= ~0xFF;

      }

    }



    element[0] = prevelement[0] + delta[0];

    element[0] &= limit;



    {
      prevelement[0] = element[0];

    }


      /* Make the element signed? */

    element[0] -= unsign;


      /* Save the element */
    *((unsigned int *) unsigned_char_data) = element[0];
    unsigned_char_data += elsize;

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
