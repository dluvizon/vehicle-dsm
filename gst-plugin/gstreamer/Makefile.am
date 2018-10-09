## run autogen.sh to create Makefile.in from this file
plugin_LTLIBRARIES = libgstplugin.la

# sources used to compile this plug-in
libgstplugin_la_SOURCES = gst-vehicle-dsm.c

# compiler and linker flags used to compile this plugin, set in configure.ac
libgstplugin_la_CFLAGS = $(GST_CFLAGS) -I${top_srcdir}/src/common

libgstplugin_la_LIBADD = \
		$(GST_LIBS) \
		$(top_srcdir)/src/libvehicledsm.la

libgstplugin_la_LDFLAGS = $(GST_PLUGIN_LDFLAGS)
libgstplugin_la_LIBTOOLFLAGS = --tag=disable-static

# headers we need but don't want installed
noinst_HEADERS = gst-vehicle-dsm.h
