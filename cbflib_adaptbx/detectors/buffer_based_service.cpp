#include <cbflib_adaptbx/buffer_based_service.h>

namespace ide = iotbx::detectors;

//Code contributed by Graeme Winter, Diamond Light Source:
typedef union {
  char b[2];
  short s;
} u_s;

typedef union {
  char b[4];
  short i;
} u_i;

// functions for byte swapping

void byte_swap_short(char * b)
{
  char c;
  c = b[0];
  b[0] = b[1];
  b[1] = c;
  return;
}

void byte_swap_int(char * b)
{
  char c;
  c = b[0];
  b[0] = b[3];
  b[3] = c;
  c = b[1];
  b[1] = b[2];
  b[2] = c;
  return;
}

// helper function: is this machine little endian? CBF files are

bool little_endian()
{
  int i = 0x1;
  char b = ((u_i *) &i)[0].b[0];
  if (b == 0)
    {
      return false;
    }
  else
    {
      return true;
    }
}

void ide::buffer_uncompress(char* packed, std::size_t packed_sz, int* values)
{
  int current = 0;
  unsigned int j = 0;
  short s;
  char c;
  int i;
  bool le = little_endian();

  while (j < packed_sz)
    {
      c = packed[j];
      j += 1;

      if (c != -128)
        {
          current += c;
          *values=current;
          values++;
          continue;
        }

      ((u_s *) & s)[0].b[0] = packed[j];
      ((u_s *) & s)[0].b[1] = packed[j + 1];
      j += 2;

      if (!le)
        {
          byte_swap_short((char *) &s);
        }

      if (s != -32768)
        {
          current += s;
          *values=current;
          values++;
          continue;
        }

      ((u_i *) & i)[0].b[0] = packed[j];
      ((u_i *) & i)[0].b[1] = packed[j + 1];
      ((u_i *) & i)[0].b[2] = packed[j + 2];
      ((u_i *) & i)[0].b[3] = packed[j + 3];
      j += 4;

      if (!le)
        {
          byte_swap_int((char *) &i);
        }

      current += i;
      *values=current;
      values++;
    }

}
