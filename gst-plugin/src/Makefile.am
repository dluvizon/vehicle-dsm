## run autogen.sh to create Makefile.in from this file
noinst_LTLIBRARIES = libvehicledsm.la

.PHONY: backgroundsub common grouping gtruth seg svm thog klt

# common subdir
COMMON_SRC = \
	common/xml/xml.cpp \
	common/image.cpp \
	common/vector.cpp \
	common/shapes.cpp \
	common/straux.cpp \
	common/fsaux.cpp \
	common/utils.cpp \
	common/clanczos.cpp \
	common/lanczos.cpp \
	common/font.cpp

COMMON_H = \
	common/xml/pugiconfig.hpp \
	common/xml/pugixml.hpp \
	common/image.h \
	common/vector.h \
	common/shapes.h \
	common/straux.h \
	common/fsaux.h \
	common/dsm-common.h \
	common/term-colors.h \
	common/debug.h \
	common/index.h \
	common/utils.h \
	common/clanczos.h \
	common/lanczos.cpp \
	common/font.h \
	common/video-defs.h

KLT_SRC = \
	klt/klt.cpp \
	klt/convolve.cpp \
	klt/error.cpp \
	klt/klt_util.cpp \
	klt/pnmio.cpp \
	klt/pyramid.cpp \
	klt/storeFeatures.cpp \
	klt/trackFeatures.cpp \
	klt/writeFeatures.cpp \
	klt/selectGoodFeatures.cpp

KLT_H = \
	klt/klt.h \
	klt/base.h \
	klt/convolve.h \
	klt/error.h \
	klt/klt_util.h \
	klt/pnmio.h \
	klt/pyramid.h

SIFT_SRC = \
	sift/key.cpp \
	sift/sfunc.cpp

SIFT_H = \
	sift/key.h

# sources used to compile this plug-in
libvehicledsm_la_SOURCES = \
	vehicle-speed.cpp \
	backgroundsub/backgroundsub.cpp \
	snoopertext/snoopertext.cpp \
	snoopertext/detection-minetto.cpp \
	snoopertext/detection-zheng.cpp \
	$(COMMON_SRC) \
	$(KLT_SRC) \
	$(SIFT_SRC) \
	grouping/grouping.cpp \
	gtruth/vehicles.cpp \
	seg/toggle.cpp \
	svm/svm.cpp \
	speed/speed.cpp \
	thog/thog.cpp

# compiler and linker flags used to compile this plugin, set in configure.ac
libvehicledsm_la_CPPFLAGS = $(CPPFLAGS) \
	-I${top_srcdir}/src \
	-I${top_srcdir}/src/common

libvehicledsm_la_LDFLAGS = -ljpeg

# headers we need but don't want installed
noinst_HEADERS = \
	vehicle-speed.h \
	backgroundsub/backgroundsub.h \
	snoopertext/snoopertext.h \
	snoopertext/detection-minetto.h \
	snoopertext/detection-zheng.h \
	lanes.h \
	$(COMMON_H) \
	$(KLT_H) \
	$(SIFT_H) \
	grouping/grouping.h \
	gtruth/vehicles.h \
	seg/toggle.h \
	svm/svm.h \
	speed/speed.h \
	thog/thog.h

