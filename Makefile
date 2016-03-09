########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

########################################################################
EXECFILE    := EXKalRTPC
LIBNAME     := EXKalRTPC
VERSION     := 1.0.1

LIBFILE     := lib$(LIBNAME).so
USERDICT    := $(LIBNAME)Dict
########################################################################
SRCDIR      := src
INCDIR      := include
OBJDIR      := obj.$(ARCH)

########################################################################
MODELDIR    := #geomlib:kallib:kaltracklib:utils
MODELLIST   := $(subst :, ,$(MODELDIR))
INCDIRS     := $(INCDIR):$(MODELDIR)

#OTHERINC just define the header file, will not be used for cint
#OTHERLIBS is for those third-party libs
ifeq ("$(KALMANROOT)","")
KALMANROOT := /media/DISK500G/work/KalmanFilter
endif
OTHERINC    := -I$(KALMANROOT)/include -I$(CLHEP_INCLUDE_DIR)
OTHERLIBS   := -L$(KALMANROOT)/lib -lS4Goem -lS4Kalman -lS4KalTrack -lS4Utils
OTHERLIBS   += -L$(CLHEP_LIB_DIR) -lCLHEP

########################################################################
# Compiler
AR          := ar
CC          := gcc
CXX         := g++
FF          := gfortran
LD          := g++

########################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
INCDIRS     := $(patsubst %,-I%,$(subst :, ,$(INCDIRS)))
INCDIRS     += $(OTHERINC)
CFLAGS      := -Wall -fPIC -O3 -g $(MODE)
CXXFLAGS    := -Wall -fPIC -O3 -g $(MODE)
FFLAGS      := -Wall -fPIC -O3 -g $(MODE)
ifeq ($(MYOS),Darwin)
#in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
               -fno-range-check
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fno-automatic -fbounds-check -funroll-all-loops \
               -fdollar-ok -ffixed-line-length-none -fno-range-check
endif
GPPFLAGS    := -MM
LDFLAGS     := -O3 -g $(MODE)

########################################################################
# Generate obj file list
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)
# header files
HEADERS     := $(foreach n,$(subst :, ,$(INCDIR)),$(wildcard $(n)/*.hh))
HEADERS     += $(foreach n,$(subst :, ,$(INCDIR)),$(wildcard $(n)/*.h))
LINKDEF     := $(foreach n,$(subst :, ,$(INCDIR)),$(wildcard $(n)/*LinkDef.h*))
#remove LinkDef.h from header because cint requires its position be at the end
#and only cint use it. No other source will use it
HEADERS     := $(filter-out $(LINKDEF), $(HEADERS))
HEADERS     := $(filter-out include/BField_Helm.hh, $(HEADERS))
# add .o to all the source files
OBJS        := $(addsuffix .o, $(basename $(SOURCES)))
OBJS        := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

########################################################################
# Libs
SYSLIBS     := -lstdc++ -lgfortran

ifdef LIBCONFIG
INCDIRS     += -I$(LIBCONFIG)/include
SYSLIBS     += -L$(LIBCONFIG)/lib -lconfig
else
#$(error $$LIBCONFIG environment variable not defined)
endif

CFLAGS      += $(INCDIRS)
CXXFLAGS    += $(INCDIRS)
FFLAGS      += $(INCDIRS)

########################################################################
# ROOT configure
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTGLIBS   := $(shell root-config --glibs) -lMinuit

CXXFLAGS    += $(ROOTCFLAGS)
LIBS        := $(SYSLIBS) $(ROOTLIBS)
GLIBS       := $(SYSLIBS) $(ROOTGLIBS)

########################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.PHONY: all clean test
VPATH       := $(SRCDIR)

########################################################################
all: lib  exe

########################################################################
# Make the $(TARGET).d file and include it.
$(OBJDIR)/%.d: %.c
	@echo Making dependency for file $< ......
	@set -e; \
	$(CC) $(GPPFLAGS) $(CFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.C
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cc
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.CC
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cpp
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cxx
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.f
#	@echo Making dependency for file $< ......
#	@set -e; \
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.F
#	@echo Making dependency for file $< ......
#	@set -e; \
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

ifneq ($(DEPS),)
-include $(DEPS)
endif

########################################################################
exe: dir $(OBJS) $(OBJDIR)/Main.o $(OBJDIR)/$(USERDICT).o
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJDIR)/Main.o $(OBJS)\
           $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@echo "Linking $(EXECFILE) ...... done!"

$(OBJDIR)/Main.o: Main.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
lib: dir $(OBJS) $(OBJDIR)/$(USERDICT).o
	@for model in $(MODELLIST); do \
		if [ -d $$model ]; then \
			make -s -C $$model; \
		fi; \
	done;
	@$(LD) -shared $(LDFLAGS) -o $(LIBFILE).$(VERSION) \
           $(OBJS) $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@ln -sf $(LIBFILE).$(VERSION) $(LIBFILE)
	@echo "Linking $(LIBFILE) ...... done!"


$(USERDICT).cxx: $(HEADERS) include/LinkDef.h
		@echo "Generating dictionary $(USERDICT).cxx ......"
		@$(ROOTSYS)/bin/rootcint -f $@ -c $(INCDIRS) $^

$(OBJDIR)/$(USERDICT).o: $(USERDICT).cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
$(OBJDIR)/%.o: %.c
	@echo Compiling $< ......
	@$(CC) -c $< -o $@  $(CFLAGS)

$(OBJDIR)/%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.CC
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

$(OBJDIR)/%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

dir:
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

########################################################################
clean: dir
	@rm -f $(OBJDIR)/*
	@rm -f $(USERDICT).cxx $(USERDICT).h
	@rm -f $(EXECFILE) $(LIBFILE) $(LIBFILE).$(VERSION)
	@rm -f *~ *# */*~ */*#

distclean: clean
	@for model in $(MODELLIST); do \
		if [[ -d $$model ]]; then \
			make distclean -s -C $$model; \
		fi; \
	done;

test:
	@echo ==================================
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo ==================================
	@echo \\LIBNAME\:$(LIBNAME) \\VERSION\:$(VERSION)
	@echo ==================================
	@echo \\USERDICT\:$(USERDICT) \\LINKDEF\:$(LINKDEF)
	@echo ==================================
	@echo \\CFLAGS\:$(CFLAGS)
	@echo ==================================	
	@echo \\CXXFLAGS\:$(CXXFLAGS)   
	@echo ==================================     
	@echo \\FFLAGS\:$(FFLAGS)
	@echo ==================================
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo ==================================
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo ==================================
	@echo \\fsources\: $(FSOURCES)	
	@echo ==================================
	@echo \\sources\: $(SOURCES)
	@echo ==================================
	@echo \\headers\: $(HEADERS)
	@echo ==================================
	@echo \\objs\: $(OBJS)	
	@echo ==================================
	@echo \\dependencies: \$(DEPS)
	@echo ==================================

help: test

env: test
