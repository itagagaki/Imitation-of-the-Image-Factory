PROGRAM  = isc
OBJS     = iscom.o subscan.o symbol.o express.o gencsg.o error.o \
           matrix.o vector.o output.o
EXTLIBS  = -lm
IDIR     = ../Hdrs
INCLUDES = -I$(IDIR)
OPTIMIZE = -O
CFLAGS   = $(OPTIMIZE) $(INCLUDES)
LD       = $(CC)
LDFLAGS  = $(OPTIMIZE)
RM       = -rm -f


$(PROGRAM) : $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(EXTLIBS)

clean :
	$(RM) $(PROGRAM) $(OBJS)

iscom.o   : $(IDIR)/config.h $(IDIR)/form.h $(IDIR)/math2.h $(IDIR)/vector.h \
            $(IDIR)/object.h $(IDIR)/light.h $(IDIR)/intcode.h \
            matrix.h tree.h symbol.h token.h subscan.h default.h \
            error.h output.h

subscan.o : $(IDIR)/config.h $(IDIR)/vector.h $(IDIR)/object.h \
            symbol.h token.h subscan.h error.h

symbol.o  : $(IDIR)/config.h $(IDIR)/vector.h $(IDIR)/object.h \
            symbol.h subscan.h

error.o   : $(IDIR)/config.h error.h

express.o : $(IDIR)/config.h $(IDIR)/vector.h $(IDIR)/object.h \
            tree.h symbol.h default.h

gencsg.o  : $(IDIR)/config.h $(IDIR)/vector.h $(IDIR)/object.h \
            matrix.h tree.h symbol.h default.h error.h

matrix.o  : $(IDIR)/config.h $(IDIR)/vector.h $(IDIR)/object.h \
            matrix.h error.h

vector.o  : $(IDIR)/config.h $(IDIR)/vector.h
