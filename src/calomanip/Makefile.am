AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  CaloManip.h \
  CaloTree.h

lib_LTLIBRARIES = \
   libcalomanip.la

libcalomanip_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib 

libcalomanip_la_SOURCES = \
  CaloManip.cc \
  CaloTree.cc
  
libcalomanip_la_LIBADD = \
  -lcalo_io \
  -lcalotrigger_io \
  -lphhepmc_io \
  -lfun4all \
  -lphool \
  -lSubsysReco


BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals 

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libcalomanip.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
