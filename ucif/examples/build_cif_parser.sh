gcc -c -I /../../antlr3/include/ \
  "../../../antlr3/src/antlr3baserecognizer.c"\
  "../../../antlr3/src/antlr3basetree.c"\
  "../../../antlr3/src/antlr3basetreeadaptor.c"\
  "../../../antlr3/src/antlr3bitset.c"\
  "../../../antlr3/src/antlr3collections.c"\
  "../../../antlr3/src/antlr3commontoken.c"\
  "../../../antlr3/src/antlr3commontree.c"\
  "../../../antlr3/src/antlr3commontreeadaptor.c"\
  "../../../antlr3/src/antlr3commontreenodestream.c"\
  "../../../antlr3/src/antlr3convertutf.c"\
  "../../../antlr3/src/antlr3cyclicdfa.c"\
  "../../../antlr3/src/antlr3debughandlers.c"\
  "../../../antlr3/src/antlr3encodings.c"\
  "../../../antlr3/src/antlr3exception.c"\
  "../../../antlr3/src/antlr3filestream.c"\
  "../../../antlr3/src/antlr3inputstream.c"\
  "../../../antlr3/src/antlr3intstream.c"\
  "../../../antlr3/src/antlr3lexer.c"\
  "../../../antlr3/src/antlr3parser.c"\
  "../../../antlr3/src/antlr3rewritestreams.c"\
  "../../../antlr3/src/antlr3string.c"\
  "../../../antlr3/src/antlr3stringstream.c"\
  "../../../antlr3/src/antlr3tokenstream.c"\
  "../../../antlr3/src/antlr3treeparser.c"\
  "../../../antlr3/src/antlr3ucs2inputstream.c"

ar -r libantlr3.a            \
antlr3baserecognizer.o       \
antlr3commontree.o           \
antlr3encodings.o            \
antlr3parser.o               \
antlr3ucs2inputstream.o      \
antlr3basetree.o             \
antlr3commontreeadaptor.o    \
antlr3exception.o            \
antlr3rewritestreams.o       \
antlr3basetreeadaptor.o      \
antlr3commontreenodestream.o \
antlr3filestream.o           \
antlr3string.o               \
antlr3bitset.o               \
antlr3convertutf.o           \
antlr3inputstream.o          \
antlr3stringstream.o         \
antlr3collections.o          \
antlr3cyclicdfa.o            \
antlr3intstream.o            \
antlr3tokenstream.o          \
antlr3commontoken.o          \
antlr3debughandlers.o        \
antlr3lexer.o                \
antlr3treeparser.o

g++ -o cif_parser -I ../../../antlr3/include/ -I ../../ main.cpp \
../cifLexer.cpp ../cifParser.cpp ../cifWalker.cpp libantlr3.a
