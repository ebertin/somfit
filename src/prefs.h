 /*
 				prefs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMFit
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for prefs.c.
*
*	Last modify:	14/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHARL	16384	/* max. nb of chars in a string list */
#define		MAXLIST		256	/* max. nb of list members */

/*--------------------------------- typedefs --------------------------------*/
typedef struct
  {
  char		*(file_name[MAXFILE]);		/* Filename(s) of input images */
  int		nfile;				/* Number of input images */
  char		nnw_name[MAXCHAR];		/* NNW filename */
  enum {SOM, BP_RETINA, BP_FILTER} nnw_type;	/* Type of NNW */
  enum {NEW, UPDATE, READONLY}	learn_type;	/* Type of learning */
  int		retisize[2], nretisize;		/* Retina size */
  double	learnrate[2];			/* Learning rate and decay */
  int		nlearnrate;			/* Number of learnparams */
  double	somkernw[2];			/* Kernel width and decay */
  int		nsomkernw;			/* Number of learnparams */
  int		niter;				/* Number of iterations */
  int		somsize[6], nsomsize;		/* SOM size */
  double	minflux;			/* Minimum flux for patterns */
  double	fwhm;				/* Typical PSF FWHM */
  int		amp_flag;			/* Amplitude free ? */
  int		spatial_mapflag;		/* Use x,y to constrain ? */
  double	spatial_stiffness[2];		/* Std Error for x,y mapping */
  int		nspatial_stiffness;		/* Number of stiff params */
  char		filter_name[MAXCHAR];		/* Filter SOM filename */
  double	filter_thresh;			/* Prefiltering threshold */
  int		filter_flag;			/* Prefiltering flag */
  char		test_name[MAXCHAR];		/* Test-image filename */
/* Check-images */
  enum {SOM_NONE, SOM_PROTO, SOM_HIT, SOM_RESIDUALS}
		check_type;			/* Check-image type */
  char		check_name[MAXCHAR];		/* Check-image filename */
/* Multithreading */
  int		nthreads;			/* Number of active threads */
/* Misc */
  enum {QUIET, NORM, LOG, FULL}	verbose_type;	/* How much it displays info */
  char		sdate_start[12];		/* SCAMP start date */
  char		stime_start[12];		/* SCAMP start time */
  char		sdate_end[12];			/* SCAMP end date */
  char		stime_end[12];			/* SCAMP end time */
  int		time_diff;			/* Execution time */
  }	prefstruct;

/*------------------------------- preferences -------------------------------*/

prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		endprefs(void),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);


#endif
