# Collie Makefile V00-04-13

# This is the configurable region
# ROOT is located through root-config; make sure it is on your PATH
# (e.g. by sourcing your ROOT installation's bin/thisroot.sh).

ROOT_CONFIG := $(shell command -v root-config 2>/dev/null)
ROOT_REQUIRED_GOALS := $(if $(MAKECMDGOALS),$(filter-out clean tidy,$(MAKECMDGOALS)),all)

ifeq ($(strip $(ROOT_CONFIG)),)
ifneq ($(strip $(ROOT_REQUIRED_GOALS)),)
$(error root-config not found on PATH; source your ROOT environment first)
endif
else
ROOTINC := $(shell $(ROOT_CONFIG) --incdir)
ROOTLIB := $(shell $(ROOT_CONFIG) --libdir)
ROOTBIN := $(shell $(ROOT_CONFIG) --bindir)
ROOTLIBS := $(shell $(ROOT_CONFIG) --libs)
ROOTCXXFLAGS := $(shell $(ROOT_CONFIG) --cflags)
endif

# Collie directories
IO_INC := io/include
IO_SRC := io/src
LIM_INC := limit/include
LIM_SRC := limit/src
COLLIELIB := lib
COLLIEEXMP := examples

# Directories to search for header files
SEARCHDIRS := -I. -I$(ROOTINC) -I$(IO_INC) -I$(LIM_INC)

# Note: -ffast-math is deliberately not used; it breaks IEEE semantics
# (e.g. NaN checks) that the statistics code relies on.
DEPENDFLAGS := -g -O3 -Wall -std=c++17 $(SEARCHDIRS) -fPIC -pthread
TOUCHHEADERS := $(MYCODEDIR)/*.h

# C
CC := gcc
CFLAGS := $(DEPENDFLAGS)

# C++
CXX := g++
CXXFLAGS := $(DEPENDFLAGS) $(ROOTCXXFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
%.o: %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# C preprocessor
CPPFLAGS :=

# Linker
LOADLIBS := -lm -ldl -pthread
LIMLIBS := -lMinuit2

EXAMPLELIBS := -lCollieLimit -lCollieIO $(ROOTLIBS)
LDFLAGS := -L$(ROOTLIB) -L$(COLLIELIB)
COLLIE_LINK_DEPS := $(COLLIELIB)/libCollieIO.so $(COLLIELIB)/libCollieLimit.so

# Canned recipe for linking example and driver executables
LINKEXE = $(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(EXAMPLELIBS) $(LOADLIBS)

LIMITSRCS := \
  $(LIM_SRC)/SigBkgdDist.cc \
  $(LIM_SRC)/CLpoint.cc \
  $(LIM_SRC)/CLcompute.cc \
  $(LIM_SRC)/FitTest.cc \
  $(LIM_SRC)/Loader.cc \
  $(LIM_SRC)/CLfit2.cc \
  $(LIM_SRC)/CLfitSS.cc \
  $(LIM_SRC)/CLfast.cc \
  $(LIM_SRC)/CLsyst.cc \
  $(LIM_SRC)/CLfit.cc \
  $(LIM_SRC)/CollieLoader.cc \
  $(LIM_SRC)/FileSetLoader.cc \
  $(LIM_SRC)/ProfileLH.cc \
  $(LIM_SRC)/ProfileLH_2D.cc \
  $(LIM_SRC)/CrossSectionLimit_FastApprox.cc \
  $(LIM_SRC)/CrossSectionLimit.cc \
  $(LIM_SRC)/ThreeSigmaEvidence.cc \
  $(LIM_SRC)/ExclusionLimit.cc \
  $(LIM_SRC)/CrossSectionCalc.cc \
  $(LIM_SRC)/timeBasedSeed.cc \
  $(LIM_SRC)/CrossSectionCalc2D.cc \
  $(LIM_SRC)/BinnedModel.cc \
  $(LIM_SRC)/ProfileLikelihoodEngine.cc \
  $(LIM_SRC)/RunContext.cc \
  $(LIM_SRC)/SigBkgdDistAdapter.cc
LIMITOBJS := $(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.cc,%.o,$(LIMITSRCS))))

LIMITHEADERS := $(LIM_INC)/CLpoint.hh

IOSRCS := \
  $(IO_SRC)/CollieDistribution.cc \
  $(IO_SRC)/CollieEventList.cc \
  $(IO_SRC)/CollieChannel.cc \
  $(IO_SRC)/CollieIterator.cc \
  $(IO_SRC)/CollieHistogram.cc \
  $(IO_SRC)/CollieQuickKEYS.cc \
  $(IO_SRC)/CollieIOFile.cc \
  $(IO_SRC)/CollieXsec.cc \
  $(IO_SRC)/CollieMasspoint.cc \
  $(IO_SRC)/CDF_DZero_Distribution.cc \
  $(IO_SRC)/CDF_DZero_IOpoint.cc \
  $(IO_SRC)/CDF_DZero_IOfile.cc \
  $(IO_SRC)/collieio_dict.cxx
IOOBJS := $(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.cc,%.o,$(IOSRCS))))

IOHEADERS := \
  $(IO_INC)/CollieDistribution.hh \
  $(IO_INC)/CollieEventList.hh \
  $(IO_INC)/CollieChannel.hh \
  $(IO_INC)/CollieIterator.hh \
  $(IO_INC)/CollieHistogram.hh \
  $(IO_INC)/CollieQuickKEYS.hh \
  $(IO_INC)/CollieIOFile.hh \
  $(IO_INC)/CollieXsec.hh \
  $(IO_INC)/CollieMasspoint.hh \
  $(IO_INC)/CDF_DZero_Distribution.hh \
  $(IO_INC)/CDF_DZero_IOpoint.hh \
  $(IO_INC)/CDF_DZero_IOfile.hh \
  $(IO_INC)/InfoForEfficiencyCalculation.hh

ALLSRCS := $(LIMITSRCS) $(IOSRCS)
ALLOBJS := $(patsubst %.cpp,%.o,$(patsubst %.cxx,%.o,$(patsubst %.cc,%.o,$(ALLSRCS))))

EXAMPLE_TARGETS := \
  examples/collieIOmaker.exe \
  examples/collieLimitCalc.exe \
  examples/collieXsecCalc.exe \
  examples/resultsFileCombiner.exe \
  examples/collieIOcondenser.exe \
  examples/collieFastApprox.exe \
  examples/collieXsecFitLoader.exe \
  examples/IO.exe \
  examples/Sigmas.exe \
  examples/CLs.exe \
  examples/IO_contour.exe \
  examples/CLsb.exe

GENERATED_DICTIONARIES := \
  $(IO_SRC)/collieio_dict.cxx \
  $(IO_SRC)/collieio_dict.h \
  $(IO_SRC)/collieio_dict.o \
  $(LIM_SRC)/collielimit_dict.cxx \
  $(LIM_SRC)/collielimit_dict.h \
  $(LIM_SRC)/collielimit_dict.o

TIDY_FILES := \
  core \
  $(ALLOBJS) \
  $(COLLIEEXMP)/*.o \
  $(COLLIEEXMP)/*.exe \
  tools/*.o \
  tools/*.exe \
  tools/collie \
  $(GENERATED_DICTIONARIES)

all: $(COLLIELIB)/libCollieIO.so $(COLLIELIB)/libCollieLimit.so tools/collie $(EXAMPLE_TARGETS)

$(IO_SRC)/collieio_dict.cxx: Makefile $(IOHEADERS)
	$(ROOTBIN)/rootcint -f $(IO_SRC)/collieio_dict.cxx -c $(SEARCHDIRS) $(IOHEADERS)

$(LIM_SRC)/collielimit_dict.cxx: Makefile $(LIMITHEADERS)
	$(ROOTBIN)/rootcint -f $(LIM_SRC)/collielimit_dict.cxx -c $(SEARCHDIRS) $(LIMITHEADERS)

$(COLLIELIB)/libCollieIO.so: $(IOOBJS) $(IO_SRC)/collieio_dict.o
	$(CXX) -shared -Wl,-soname,libCollieIO.so -o $@ $^ $(LOADLIBS) $(LDFLAGS)

$(COLLIELIB)/libCollieLimit.so: \
    $(LIMITOBJS) $(LIM_SRC)/collielimit_dict.o | $(COLLIELIB)/libCollieIO.so
	$(CXX) -shared -Wl,-soname,libCollieLimit.so -o $@ $^ $(LIMLIBS) -lCollieIO $(ROOTLIBS) $(LOADLIBS) $(LDFLAGS)

tools/collie: tools/collie.o tools/CollieCli.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/resultsFileCombiner.exe: examples/resultsFileCombiner.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/collieIOcondenser.exe: examples/collieIOcondenser.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/collieXsecFitLoader.exe: examples/xsecFitLoader.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/collieIOmaker.exe: examples/collieIOexample.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/collieLimitCalc.exe: examples/exampleLimitCalculation.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/collieFastApprox.exe: examples/exampleFastApproxLimitCalculation.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/collieXsecCalc.exe: examples/exampleXsecMeasurement.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/IO.exe: examples/IO.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/IO_contour.exe: examples/IO_contour.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/Sigmas.exe: examples/Sigmas.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/CLsb.exe: examples/CLsb.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

examples/CLs.exe: examples/CLs.o | $(COLLIE_LINK_DEPS)
	$(LINKEXE)

tidy::
	@$(RM) $(TIDY_FILES)

clean:: tidy
	@$(RM) $(COLLIELIB)/libCollieIO.so $(COLLIELIB)/libCollieLimit.so make.depend make.depend.bak

depend: make.depend

make.depend: Makefile
	touch make.depend
	makedepend -K -fmake.depend -- $(CXXFLAGS) -- $(ALLSRCS) 2>/dev/null

# include make.depend
# DO NOT DELETE
