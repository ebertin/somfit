# Program Makefile for SOMFit
# Copyright (C) 1996-2009 Emmanuel Bertin.
SUBDIRS		= fits

bin_PROGRAMS	= somfit
somfit_SOURCES	= main.c makeit.c prefs.c som.c vignet.c \
		  define.h globals.h key.h preflist.h prefs.h som.h \
		  threads.h types.h vignet.h

somfit_LDADD	= $(top_builddir)/src/fits/libfits.a
DATE=`date +"%Y-%m-%d"`
