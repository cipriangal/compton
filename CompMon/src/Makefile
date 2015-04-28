# -------------------------------------------------------------------------
#   Sources for all modules
# -------------------------------------------------------------------------

SRC	= mainanal_CompMon.cc fadcTriggered.cc fadcdata.cc vmeauxdata.cc comptonStatus.cc comptonParams.cc fadcAccums.cc THaEpics.C

OBJ     = $(SRC:.cc=.o) THaCodaFile.o THaCodaData.o

NAME    = compmon

# -------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -------------------------------------------------------------------------

BINDIR		= ../

LIBS		= $(shell $(ROOTSYS)/bin/root-config --libs)

GLIBS		= $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11

EVIO		= /home/compton/evio/Linux-x86_64/lib
EVIO_LIB	= $(EVIO)/libevio.so

ALL_LIBS	= $(EVIO_LIB) $(GLIBS) $(ROOTLIBS)

CC		= g++ -D_FILE_OFFSET_BITS=64

ifdef ONLINE
	ALL_LIBS += $(ONLIBS)
endif

# -------------------------------------------------------------------------
#   Compiler & Linker flags
# -------------------------------------------------------------------------
ROOT_FLAGS      = $(shell $(ROOTSYS)/bin/root-config --cflags)

CXXFLAGS	= -g -Wall $(ROOT_FLAGS) $(LIBS) $(INCLUDES) -DSTANDALONE

# =========================================================================
#   Targets & Rules
# =========================================================================
all: object

object: Makefile $(OBJ)
	$(CC)  $(CXXFLAGS) -o $(BINDIR)/$(NAME) $(OBJ)  $(ALL_LIBS)

clean:
	rm -fv *.o
