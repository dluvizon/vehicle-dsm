dnl required version of autoconf
AC_PREREQ([2.53])
AC_INIT([vdsm-package], [1.02], [diogo@luvizon.com])

# Release versioning
GENERIC_MAJOR_VERSION=1
GENERIC_MINOR_VERSION=2
GENERIC_MICRO_VERSION=0
PACKAGE_YEAR=2015
PACKAGE_DATE="03/2015"

dnl required versions of gstreamer and plugins-base
GST_REQUIRED=1.0.0
GSTPB_REQUIRED=1.0.0

AC_CONFIG_MACRO_DIR([m4])
AC_REVISION([$Revision$])

AC_CONFIG_AUX_DIR(config)

AC_DEFINE_UNQUOTED(PACKAGE_NAME,["${PACKAGE_NAME}"],[Name of package])
AC_DEFINE_UNQUOTED(PACKAGE_VERSION,["${PACKAGE_VERSION}"],[Version number])
AC_DEFINE_UNQUOTED(PACKAGE_YEAR,"$PACKAGE_YEAR",[Official year for this release])
AC_DEFINE_UNQUOTED(PACKAGE_DATE,"$PACKAGE_DATE",[Official date of release])

AC_SUBST(PACKAGE_NAME)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(PACKAGE_YEAR)
AC_SUBST(PACKAGE_DATE)

# API version (often = GENERIC_MAJOR_VERSION.GENERIC_MINOR_VERSION)
GENERIC_API_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
GENERIC_LIBRARY_VERSION=$GENERIC_MAJOR_VERSION:$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_API_VERSION)
AC_SUBST(GENERIC_MAJOR_VERSION)

AC_SUBST(GENERIC_LIBRARY_VERSION)
PACKAGE=$GENERIC_LIBRARY_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)

GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_VERSION)

# ----------------------------------------
# Automake configuration
# ----------------------------------------

AM_INIT_AUTOMAKE([subdir-objects])
AM_PATH_PYTHON([2.5])
AC_CONFIG_HEADERS(config.h:config/config.h.in)
AC_CONFIG_SRCDIR([gstreamer/gst-vehicle-dsm.c])
AM_MAINTAINER_MODE

AC_DISABLE_STATIC
AM_PROG_LIBTOOL
AC_LIBTOOL_CXX

dnl enable mainainer mode by default
AM_MAINTAINER_MODE([enable])

dnl check for tools (compiler etc.)
AC_LANG(C)
AC_LANG(C++)
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_MKDIR_P

AC_SUBST(DEPS_CFLAGS)
AC_SUBST(DEPS_CXXFLAGS)
AC_SUBST(DEPS_LIBS)

dnl required version of libtool
LT_PREREQ([2.2.6])
LT_INIT

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strstr])

dnl give error and exit if we don't have pkgconfig
AC_CHECK_PROG(HAVE_PKGCONFIG, pkg-config, [ ], [
  AC_MSG_ERROR([You need to have pkg-config installed!])
])

dnl Check for the required version of GStreamer core (and gst-plugins-base)
dnl This will export GST_CFLAGS and GST_LIBS variables for use in Makefile.am
dnl
dnl If you need libraries from gst-plugins-base here, also add:
dnl for libgstaudio-1.0: gstreamer-audio-1.0 >= $GST_REQUIRED
dnl for libgstvideo-1.0: gstreamer-video-1.0 >= $GST_REQUIRED
dnl for libgsttag-1.0: gstreamer-tag-1.0 >= $GST_REQUIRED
dnl for libgstpbutils-1.0: gstreamer-pbutils-1.0 >= $GST_REQUIRED
dnl for libgstfft-1.0: gstreamer-fft-1.0 >= $GST_REQUIRED
dnl for libgstinterfaces-1.0: gstreamer-interfaces-1.0 >= $GST_REQUIRED
dnl for libgstrtp-1.0: gstreamer-rtp-1.0 >= $GST_REQUIRED
dnl for libgstrtsp-1.0: gstreamer-rtsp-1.0 >= $GST_REQUIRED
dnl etc.
PKG_CHECK_MODULES(GST, [
  gstreamer-1.0 >= $GST_REQUIRED
  gstreamer-base-1.0 >= $GST_REQUIRED
  gstreamer-controller-1.0 >= $GST_REQUIRED
], [
  AC_SUBST(GST_CFLAGS)
  AC_SUBST(GST_LIBS)
], [
  AC_MSG_ERROR([
      You need to install or upgrade the GStreamer development
      packages on your system. On debian-based systems these are
      libgstreamer1.0-dev and libgstreamer-plugins-base1.0-dev.
      on RPM-based systems gstreamer1.0-devel, libgstreamer1.0-devel
      or similar. The minimum version required is $GST_REQUIRED.
  ])
])

dnl check if compiler understands -Wall (if yes, add -Wall to GST_CFLAGS)
AC_MSG_CHECKING([to see if compiler understands -Wall])
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -Wall -O2"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([ ], [ ])], [
  GST_CFLAGS="$GST_CFLAGS -Wall"
  AC_MSG_RESULT([yes])
], [
  AC_MSG_RESULT([no])
])

dnl set the plugindir where plugins should be installed (for src/Makefile.am)
if test "x${prefix}" = "x$HOME"; then
  plugindir="$HOME/.gstreamer-1.0/plugins"
else
  plugindir="\$(libdir)/gstreamer-1.0"
fi
AC_SUBST(plugindir)

dnl set proper LDFLAGS for plugins
GST_PLUGIN_LDFLAGS='-module -avoid-version -export-symbols-regex [_]*\(gst_\|Gst\|GST_\).* -lstdc++'
AC_SUBST(GST_PLUGIN_LDFLAGS)

AC_CONFIG_FILES([Makefile
		 src/Makefile
		 tools/Makefile
		 tools/train/Makefile
		 tools/perspective/Makefile
		 gstreamer/Makefile])

AC_OUTPUT

# Final menssage
echo ""
echo "Configuration is done."

# ----------------------------------------
# CONFIG Template
# ----------------------------------------

# Fence added in configuration file
AH_TOP([
#ifndef __VDSM_CONFIG_H
#define __VDSM_CONFIG_H

/* config.h: begin */
/* This file is generated automatically by configure.ac. Do not edit it! */
])

# Stuff added at bottom of file
AH_BOTTOM([

/* Miscellaneous defines */
#define AUTOCONF 1

/* Define the maximum number of lanes to be monitored. */
#define GTRUTH_MAX_LANES 3

/* By default, enable debug in verbose level */
#define VDSM_DLEVEL 1

/**
 * Enable or disable the debug for each module here:
 */
#define BACKGROUNDSUB_DLEVEL 2
#define GROUNDTRUTH_DLEVEL 1
#define SNOOPERTEXT_DLEVEL 2

/**
 * Enable or disable the generation of log files.
 */
#define VDSM_EN_LOG 1

/* config.h: end */
#endif /* __VDSM_CONFIG_H */
])
