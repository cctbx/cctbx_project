from libtbx import lzw
enc = lzw.ByteEncoder(12)
bigstr = b"gabba gabba yo gabba gabba gabba yo gabba gabba gabba yo gabba gabba gabba yo"
encoding = enc.encodetobytes(bigstr)
encoded = b"".join( b for b in encoding )
#print encoded
#    '3\\x98LF#\\x08\\x82\\x05\\x04\\x83\\x1eM\\xf0x\\x1c\\x16\\x1b\\t\\x88C\\xe1q(4"\\x1f\\x17\\x85C#1X\\xec.\\x00'

dec = lzw.ByteDecoder()
decoding = dec.decodefrombytes(encoded)
decoded = b"".join(decoding)
assert decoded == bigstr
print "OK"
