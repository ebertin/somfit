 /*
 				define.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMFit
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global definitions.
*
*	Last modify:	15/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define         VERSION         "1.x"
#define         DATE            "2006-03-14"
#define		THREADS_NMAX	16              /* max. number of threads */
#endif


/*------------------------ what, who, when and where ------------------------*/

#define         BANNER		"SOMFit"
#ifdef USE_THREADS
#define         MYVERSION       VERSION "-MP"
#else
#define         MYVERSION       VERSION
#endif
#define         COPYRIGHT	"Emmanuel BERTIN <bertin@iap.fr>"
#define         INSTITUTE	"TERAPIX team at IAP  http://terapix.iap.fr"

/*----------------------------- Physical constants --------------------------*/

#ifndef PI
#define PI      	3.1415926535898
#endif

/*----------------------------- Internal constants --------------------------*/

#define		BIG		1e+30		/* a huge number */
#define		OUTPUT		stdout		/* where all msgs are sent */
#define		MAXCHAR		512		/* max. number of characters */
#define		MAXFILE		32768		/* max number of input files */

/*------------ Set defines according to machine's specificities -------------*/

#if 0
#define	NO_ENVVAR
#endif

/*--------------------- in case of missing constants ------------------------*/

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef EXIT_SUCCESS
#define 	EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define		EXIT_FAILURE	-1
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*------------------------------- Other Macros ------------------------------*/

#define	DEXP(x)	exp(2.30258509299*(x))	/* 10^x */

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (fseek(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(pos, afile, fname) \
		if ((pos=ftell(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ)))) \
		   error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
                    error(EXIT_FAILURE, "Not enough memory for ", \
                        #ptrout " (" #nel " elements) !"); \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ));};;}

#define QPOPEN(file, cmdline, flag) \
		{if (!(file=popen(cmdline, flag))) \
		  error(EXIT_FAILURE, "*Error*: cannot execute ", cmdline);;}

#define	RINT(x)	(int)(floor(x+0.5))

#define	NPRINTF		if (prefs.verbose_type == NORM) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
			  fprintf(w, "\33[1M> %s\n\33[1A",x);}

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#define QPRINTF		if (prefs.verbose_type != QUIET)	fprintf

#define QIPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[7m%s\33[0m\n", x); \
			else if (prefs.verbose_type == LOG) \
				fprintf(w, "%s\n", x);}

#define QBPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[01;31m%s\33[0m\n", x); \
			else if (prefs.verbose_type == LOG) \
				fprintf(w, "%s\n", x);}

