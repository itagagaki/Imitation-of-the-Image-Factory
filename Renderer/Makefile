PROGRAM   = is
OBJS      = is.o scan.o input.o
EXTLIBS   = -lm
IDIR      = ../Hdrs
INCLUDES  = -I$(IDIR)
OPTIMIZE  = -O3
DEBUGFLAG = -full_report
CFLAGS    = $(OPTIMIZE) $(DEBUGFLAG) $(INCLUDES)
LD        = $(CC)
LDFLAGS   = $(OPTIMIZE)
RM        = rm -f

$(PROGRAM) : $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(EXTLIBS)

clean :
	$(RM) $(PROGRAM) $(OBJS)

is.o      : $(IDIR)/config.h $(IDIR)/form.h $(IDIR)/vector.h \
            $(IDIR)/object.h $(IDIR)/light.h $(IDIR)/intcode.h \
            input.h

scan.o    : $(IDIR)/config.h $(IDIR)/form.h $(IDIR)/math2.h \
            $(IDIR)/vector.h $(IDIR)/object.h $(IDIR)/light.h \
            toleran.h

scan2.o   : $(IDIR)/config.h $(IDIR)/form.h $(IDIR)/vector.h \
            $(IDIR)/object.h
