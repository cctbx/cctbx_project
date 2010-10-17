set_lower = set("""\
abs
acos
aimag
aint
alog
alog10
amax0
amax1
amin0
amin1
amod
anint
asin
atan
atan2
cabs
cdabs
char
cmplx
conjg
cos
cosh
dabs
dacos
dasin
datan2
dble
dcmplx
dconjg
dcos
dexp
dim
dimag
dlog
dlog10
dmax1
dmin1
dprod
dsign
dsin
dsqrt
dtan
exp
float
iabs
iargc
iand
ichar
idnint
index
int
ishft
len
len_trim
lge
lgt
lle
llt
lnblnk
log
log10
max
max0
max1
maxloc
min
min0
min1
mod
nint
real
sign
sin
sinh
sngl
sqrt
tan
tanh
transfer
""".splitlines())

extra_set_lower = set("""\
cpu_time
date
getarg
getenv
system
time
""".splitlines())

io_set_lower = set("""\
flush
""".splitlines())

is_modified_info_by_name = {
  "cpu_time": (0,),
  "date": (0,),
  "getarg": (1,),
  "getenv": (1,),
  "time": (0,)
}
