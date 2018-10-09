bin_PROGRAMS = train svm-train

train_SOURCES = train.cpp png.cpp
train_CPPFLAGS = \
	-I${top_srcdir}/src \
	-I${top_srcdir}/src/common
train_LDADD = ${top_srcdir}/src/libvehicledsm.la

svm_train_SOURCES = svm-train.cpp
svm_train_CPPFLAGS = -I${top_srcdir}/src
svm_train_LDADD = ${top_srcdir}/src/libvehicledsm.la

noinst_HEADERS = png.h
