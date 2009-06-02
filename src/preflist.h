 /*
 				preflist.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMFit
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Keywords for the configuration file.
*
*	Last modify:	15/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include "key.h"

#ifndef _PREFS_H_
#include "prefs.h"
#endif

#ifdef	USE_THREADS
#define	THREADS_PREFMAX THREADS_NMAX
#else
#define	THREADS_PREFMAX 65535
#endif

int idummy;

pkeystruct key[] =
 {
  {"AMP_MAPPING", P_BOOL, &prefs.amp_flag},
  {"CHECKIMAGE_NAME", P_STRING, prefs.check_name},
  {"CHECKIMAGE_TYPE", P_KEY, &prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "PROTOTYPES", "HITS", "RESIDUALS",""}},
  {"FILTER", P_BOOL, &prefs.filter_flag},
  {"FILTER_NAME", P_STRING, prefs.filter_name},
  {"FILTER_THRESH", P_FLOAT, &prefs.filter_thresh, 0,0, 0.0, 1e6},
  {"KERNEL_WIDTH", P_FLOATLIST, prefs.somkernw, 0,0, 0.0, 1e6, {""},
     2,2, &prefs.nsomkernw},
  {"LEARN_TYPE", P_KEY, &prefs.learn_type, 0,0, 0.0,0.0,
   {"NEW","UPDATE", "READ_ONLY", ""}},
  {"LEARNING_RATE", P_FLOATLIST, prefs.learnrate, 0,0, 1e-6,1e6, {""},
     2,2, &prefs.nlearnrate},
  {"MINFLUX", P_FLOAT, &prefs.minflux, 0,0, 1e-6,1e15},
  {"NN_NAME", P_STRING, prefs.nnw_name},
  {"NN_TYPE", P_KEY, &prefs.nnw_type, 0,0, 0.0,0.0,
   {"SOM","BP_RETINA", "BP_FILTER", ""}},
  {"NTHREADS", P_INT, &prefs.nthreads, 0, THREADS_PREFMAX},
  {"PASSES", P_INT, &prefs.niter, 1,1000000000},
  {"PSF_FWHM", P_FLOAT, &prefs.fwhm, 0,0, 1e-1,1e3},
  {"RETINA_SIZE", P_INTLIST, prefs.retisize, 1,1024, 0.0,0.0, {""},
     1,2, &prefs.nretisize},
  {"SOM_SIZE", P_INTLIST, prefs.somsize, 1,1024, 0.0,0.0, {""},
     1,4, &prefs.nsomsize},
  {"TESTIMAGE_NAME", P_STRING, prefs.test_name},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","LOG", "FULL",""}},
  {"XY_MAPPING", P_BOOL, &prefs.spatial_mapflag},
  {"XY_STIFFNESS", P_FLOATLIST, prefs.spatial_stiffness, 0,0, 0.0,1e6, {""},
     1,2, &prefs.nspatial_stiffness},
  {""}
 };

char			keylist[sizeof(key)/sizeof(pkeystruct)][32];
const char		notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
" ",
"#--------------------------------- learning -----------------------------------",
" ",
"NN_NAME         default.som     # Name of the file containing the NN weights",
"*NN_TYPE         SOM             # Type of Neural Network",
"*LEARN_TYPE      NEW             # NEW, UPDATE or READ_ONLY",
"LEARNING_RATE   0.1,10000       # Learning rate and decay",
"KERNEL_WIDTH    3,20000       # Kernel width and decay",
"PASSES          50000           # Number of passes through all sources",
" ",
"#------------------------------- architecture ---------------------------------",
" ",
"RETINA_SIZE     21,21            # Size of sensitive area",
"SOM_SIZE        16,12             # Size of the Self Organizing Map",
" ",
"#---------------------------- degrees of freedom ------------------------------",
" ",
"AMP_MAPPING     N               # Amplitude is constrained?",
"XY_MAPPING      N               # Add X,Y mapping constraint?",
"XY_STIFFNESS    0.2             # X,Y stiffness (in image-size units)",
" ",
"#-------------------------------- selection -----------------------------------",
" ",
"MINFLUX         1000.0         # Minimum flux (ADU) for a source to be used",
"PSF_FWHM        3.0             # Typical PSF FWHM",
"FILTER          N               # Apply filtering with another SOM (Y/N)?",
"FILTER_NAME     filter.som      # Filtering SOM filename",
"FILTER_THRESH   3               # Filtering threshold",
" ",
"#------------------------------- check-images ---------------------------------",
" ",
"CHECKIMAGE_TYPE PROTOTYPES      # Type of check-image: NONE, PROTOTYPES, HITS",
"                                # or RESIDUALS",
"CHECKIMAGE_NAME check.fits      # File name of the check-image",
"*TESTIMAGE_NAME  test.fits       # Name of a test-image",
" ",
"#------------------------------ Miscellaneous ---------------------------------",
" ",
"VERBOSE_TYPE    NORMAL          # QUIET,NORMAL, LOG or FULL",
#ifdef USE_THREADS
"NTHREADS        0               # Number of simultaneous threads for",
"                                # the SMP version of " BANNER,
"                                # 0 = automatic",
#else
"NTHREADS        1               # 1 single thread",
#endif
""
 };

