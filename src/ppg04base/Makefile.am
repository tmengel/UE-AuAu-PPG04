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
  CaloSpy.h \
  CaloTowerManip.h \
  PPG04AnaWriter.h

lib_LTLIBRARIES = \
   libppg04base.la

libppg04base_la_SOURCES = \
  CaloSpy.cc \
  CaloTowerManip.cc \
  PPG04AnaWriter.cc
  
libppg04base_la_LIBADD = \
  -lcalo_io \
  -lcalotrigger_io \
  -lcentrality_io \
  -lunderlyingevent_io \
  -lunderlyingevent \
  -ljetbase \
  -ljetbackground_io \
  -ljetbackground \
  -lfun4all \
  -lglobalvertex_io \
  -lphhepmc_io \
  -lphool \
  -lSubsysReco


BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals 

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libppg04base.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
