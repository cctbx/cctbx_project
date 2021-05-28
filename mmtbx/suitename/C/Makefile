# suitename

ifeq ($(MAKECMDGOALS),debug)
	CFLAGS = -g -isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 -arch ppc
else
	CFLAGS = -isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 -arch ppc
endif

LIBS = -lm -Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk -arch i386 -arch ppc

# ---------------------------------------------------------------------
PROG_FLGS   =  
OBJS = suitename.o suitenscrt.o suiteninit.o suiteninpt.o \
       suitenout.o suitenutil.o

# ---------------------------------------------------------------------
HEADERS = suitename.h suitenscrt.h suitenutil.h suiteninit.h suiteninpt.h \
          suitenout.h

suitename: $(OBJS)
	cc -o suitename $(CFLAGS) $(OBJS) $(LIBS)

debug:     $(OBJS)
	cc -o suitename $(CFLAGS) $(OBJS) $(LIBS)

clean:
	rm -f *.o
# ------------------------------------------------------------------------
# Dependencies  (presume .o<-.c by standard cc compiler)

suitename.o:  $(HEADERS)
suiteninit.o: $(HEADERS)
suitenscrt.o: $(HEADERS)
suiteninpt.o: $(HEADERS)
suitenout.o:  $(HEADERS)
suitenutil.o: $(HEADERS)

