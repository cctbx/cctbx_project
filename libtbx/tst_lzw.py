def run(args):
  assert len(args) == 0
  from libtbx import lzw
  enc = lzw.ByteEncoder(12)
  bigstr = "gabba gabba yo gabba gabba gabba yo gabba gabba gabba" \
    " yo gabba gabba gabba yo"
  encoding = enc.encodetobytes(bigstr)
  encoded = "".join([b for b in encoding])
  dec = lzw.ByteDecoder()
  decoding = dec.decodefrombytes(encoded)
  decoded = "".join(decoding)
  assert decoded == bigstr
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
