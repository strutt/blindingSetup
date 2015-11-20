#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include ../Makefile.arch

#Site Specific  Flags
SYSINCLUDES	= -I/usr/local/include
SYSLIBS         = -L/usr/local/include
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}


ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR) -lBensAnitaTools -lAnitaCorrelator -lAnitaEvent 
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
ANITA_UTIL_CALIB_DIR=$(ANITA_UTIL_INSTALL_DIR)/share/anitaCalib
else
ANITA_UTIL_LIB_DIR=/usr/local/lib
ANITA_UTIL_INC_DIR=/usr/local/include
ANITA_UTIL_CALIB_DIR=/usr/local/share/anitaCalib
ifdef EVENT_READER_DIR
LD_ANITA_UTIL=-L$(EVENT_READER_DIR)  -lAnitaEvent
INC_ANITA_UTIL=-I$(EVENT_READER_DIR)
endif
endif

#Toggles the FFT functions on and off
USE_FFT_TOOLS=1

ifdef USE_FFT_TOOLS
FFTLIBS = -L/usr/local/lib -lRootFftwWrapper -lfftw3
FFTFLAG = -DUSE_FFT_TOOLS
else
FFTLIBS =
FFTFLAG =
endif

#Generic and Site Specific Flags
CXXFLAGS     = -g -O2 -fPIC $(ROOTCFLAGS) $(FFTFLAG) $(SYSINCLUDES) $(INC_ANITA_UTIL) #-std=c++11 
LDFLAGS      = -g


LIBS          = $(ROOTLIBS) -lMathMore -lMinuit -lGpad $(FFTLIBS) $(SYSLIBS) $(LD_ANITA_UTIL) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Toggles google performance profile functionality on and off
#USE_GPERFTOOLS=1

ifdef USE_GPERFTOOLS
LDFLAGS	+= -Wl,-no_pie -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
LIBS += -lprofiler -ltcmalloc
endif

# I use git svn and am too infrequent a commiter by nature...
# Use this flag to prompt a local commit after every compile
# (You probably don't actually want to commit after every compile, but the prompt is nice)
#FORCE_GIT=1 


# For those who like really bloody pedantic compiler warnings... like me
HARDCORE_MODE=1
ifdef HARDCORE_MODE
CXXFLAGS += -Wall -Wextra -Wshadow -Werror #-Wpedantic
endif


OBJS = 


BINARIES = makeTreesOfWaisPulsesWithSwappedPolarizations overwriteSoftwareTriggeredEventsWithSwappedWaisPulses

#Now the bits we're actually compiling
all: $(OBJS) $(BINARIES) commit

.PHONY: commit clean

$(BINARIES): %: %.$(SRCSUF) #$(ROOT_LIBRARY) 
	@echo "<**Compiling**> "
	@echo $<
	-$(LD) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $< -o $@
ifdef FORCE_GIT
	-@if test $$? == 0; then git add $<; fi
endif

%.$(OBJSUF) : %.$(SRCSUF) %.h
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o $@
ifdef FORCE_GIT	
	@if test $$? == 0; then git add $^; fi
endif

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(OBJS)
	@rm -f $(BINARIES)

commit: 
ifdef FORCE_GIT
	-@git add Makefile
	-@git commit
endif

check-syntax:
	gcc -o nul -S ${CHK_SOURCES}
