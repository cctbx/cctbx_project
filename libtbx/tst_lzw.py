from __future__ import absolute_import, division, print_function
def run(args):
  assert len(args) == 0
  from libtbx import lzw
  enc = lzw.ByteEncoder(12)
  bigstr = b"gabba gabba yo gabba gabba gabba yo gabba gabba gabba" \
    b" yo gabba gabba gabba yo"
  encoding = enc.encodetobytes(bigstr)
  encoded = b"".join([b for b in encoding])
  assert(len(bigstr) > len(encoded))
  dec = lzw.ByteDecoder()
  decoding = dec.decodefrombytes(encoded)
  decoded = b"".join(decoding)
  assert decoded == bigstr
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
