 /*
 				som.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	A program using Self-Organizing Maps (SOMs)
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for Kohonen's Self Organizing Map (V2.0).
*
*	Last modify:	14/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _SOM_H_
#define _SOM_H_

/*--------------------------------- constants -------------------------------*/

#define	INPUT_MAXDIM		9	/* Maximum dimensionality of input */
#define	SOM_MAXDIM		6	/* Maximum dimensionality of the SOM */

/*---------------------------------- macros ---------------------------------*/

#define	SOM_KERNEL(r, h)	((dtemp=(r)/(h))<70.0? exp(-dtemp):0.0)

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  int		inputdim;		/* Dimensionality of input vector */
  int		*inputsize;		/* Dimensions of the input vector */
  int		nextrainput;		/* number of extra inputs */
  int		ninput;			/* Total number of inputs */
  int		neurdim;		/* Dimensionality of the SOM */
  int		*neursize;		/* Dimensions of the SOM */
  int		nneur;			/* Total number of neurons */
  int		*neurstep;		/* Offset for each dimension */
  float		*weight;		/* Weights */
  float		*residuals;		/* Residuals */
  int		nweight;		/* Total number of weights */
  float		learnrate, clearnrate;	/* Starting and current learn. rates */
  float		learndecay;		/* Learning decay rate */
  float		kernw, ckernw;		/* Starting and current kernel width */
  float		kernwdecay;		/* Kernel width decay rate */
  float		xy_stiff;		/* Stiffness of the X/Y mapping */
  double	err;			/* Global variance */
  double	errpix, errextra;	/* Retina and extra variances */
  int		*freq;			/* Number of winning times per node */
  float		*stderror;		/* chi2 array of the linear fits */
  float		*amp;			/* Array of the fitted amplitudes */
  float		*sigamp;		/* stderr array of the fitted ampl. */
  int		ntrain;			/* # of training examples so far */
  int		nsweep;			/* # of sweeps through the whole set */
  }	somstruct;


/*---------------------------------- protos --------------------------------*/

extern somstruct	*som_init(int *inputsize, int inputdim, int nextra),
			*som_load(char *filename);

extern int		som_win(somstruct *som, float *input, float *inputw,
			float *errp, float *amp);

extern void		som_end(somstruct *som),
			som_learn(somstruct *som, float *lpatt, float *wpatt,
			float globweight, int ampflag),
			som_save(somstruct *som, char *filename),
			som_writecheck(somstruct *som, char *filename);

#endif

