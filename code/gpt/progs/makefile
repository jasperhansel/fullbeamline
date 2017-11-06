# UNIX makefile for PROGS

include ../makeinc.gcc

TARGETS = gdfa$(EXEEXT)
all: $(TARGETS)

include makefile.pro
PROGOBJSO = $(PROGOBJS:.obj=.o)

CFLAGS += -I../gdfa -I.

gdfa.c: proglist
	../bin/makeprogs ../gdfa/progkrnl proglist
	@echo Please rerun make.

gdfa: gdfa.o $(PROGOBJSO) ../gdfa/progs.a ../gdfa/gdfa.a ../kernel/gdf.a ../kernel/utils.a
	$(CC) $(LDFLAGS) $(TARGET_ARCH) -o $@ gdfa.o $(PROGOBJSO) ../gdfa/progs.a ../gdfa/gdfa.a ../kernel/gdf.a ../kernel/utils.a
	$(CP) gdfa ../bin

clean:
	$(RM) gdfa $(PROGOBJSO) gdfa.c
