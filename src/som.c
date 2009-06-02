  /*
 				som.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	A program using Self-Organizing Maps (SOMs)
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Implementation of Kohonen's Self Organizing Map (V2.0).
*
*	Last modify:	15/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef USE_THREADS
#include <pthread.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "som.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

double	*sx,*sxx,*sxy;
float	*yweight;

/********************************** som_init *********************************/
/*
Initialization of the SOM.
*/
somstruct	*som_init(int *inputsize, int inputdim, int nextra)
  {
   static char	str[MAXCHAR];
   somstruct	*som;
   float	*w, norm;
   int		i;

/* Allocate memory*/
  QCALLOC(som, somstruct, 1);
  som->clearnrate = som->learnrate = prefs.learnrate[0];
  som->learndecay = prefs.learnrate[1];
  som->ckernw = som->kernw = prefs.somkernw[0];
  som->kernwdecay = prefs.somkernw[1];
  som->xy_stiff = prefs.spatial_stiffness[0];
  som->neurdim = prefs.nsomsize;
  if (som->neurdim>SOM_MAXDIM)
    {
    sprintf(str, "%d", SOM_MAXDIM);
    error(EXIT_FAILURE, "*Error*: This package is presently limited to SOMs"
	"with dimensionality less than or equal to ", str);
    }
  QMALLOC(som->neursize, int, SOM_MAXDIM);
  QMALLOC(som->neurstep, int, SOM_MAXDIM+1);
  for (i=0; i<SOM_MAXDIM; i++)
    som->neursize[i] = 1;
  som->neurstep[0] = som->nneur = 1;
  for (i=0; i<som->neurdim; i++)
    {
    som->nneur *= (som->neursize[i] = prefs.somsize[i]);
    som->neurstep[i+1] = som->nneur;
    }
  som->inputdim = inputdim;
  if (som->inputdim>INPUT_MAXDIM)
    {
    sprintf(str, "%d", INPUT_MAXDIM);
    error(EXIT_FAILURE, "*Error*: This package is presently limited to inputs"
	"with dimensionality less or equal to ", str);
    }
  QMALLOC(som->inputsize, int, INPUT_MAXDIM);
  for (i=0; i<INPUT_MAXDIM; i++)
    som->inputsize[i] = 1;
  som->ninput = 1;
  for (i=0; i<inputdim; i++)
    som->ninput *= (som->inputsize[i] = inputsize[i]);
  som->ninput += (som->nextrainput = nextra);
  som->nweight = som->nneur*som->ninput;
  som->err = som->errpix = som->errextra = 0.0;
  som->ntrain = 0;
  QMALLOC(som->weight, float, som->nweight);
  QCALLOC(som->residuals, float, som->nweight);
  QCALLOC(som->freq, int, som->nneur);
  QMALLOC(som->stderror, float, som->nneur);
  QMALLOC(som->amp, float, som->nneur);
  QMALLOC(som->sigamp, float, som->nneur);
/* Locals */
  QMALLOC(sx, double, som->nneur);
  QMALLOC(sxx, double, som->nneur);
  QMALLOC(sxy, double, som->nneur);
  QMALLOC(yweight, float, som->ninput);

/* Initialize weights to random values*/
  srand((int)(time(NULL))%RAND_MAX);
  w = som->weight;
  norm = som->ninput/2.0;
  for (i=som->nweight; i--;)
    *(w++) = rand()/RAND_MAX/norm;

  return som;
  }


/********************************* som_learn *********************************/
/*
Make the SOM learn a new pattern.
*/
void	som_learn(somstruct *som, float *input, float *inputw, float gweight,
		int ampflag)
  {
   static int	xj[SOM_MAXDIM], xwin[SOM_MAXDIM], dist[SOM_MAXDIM];
   double	amp,dsig2, dtemp;
   float	*w, *inputt, h, norm;
   int		i,j,n,nd,ndt, nwin, ninput, diff,
		*nx, *mul,*xjt,*xwint,*dt;

/* Get the winning neuron */
  nwin = som_win(som, input, inputw, NULL, ampflag?NULL:&norm);

/* Express the winner number in SOM-space coordinates */
  mul = som->neurstep;
  nx = som->neursize;
  nd = som->neurdim;
  xwint = xwin;
  for (i=nd; i--;)
    *(xwint++) = (float)((nwin/(*(mul++)))%(*(nx++)));

/* Update learning parameters */
  som->clearnrate = som->learnrate*exp(-som->ntrain/som->learndecay);
  som->ckernw = som->kernw*exp(-som->ntrain/som->kernwdecay);
  amp = gweight*som->clearnrate;
  dsig2 = 2.0*som->ckernw*som->ckernw;
  if (dsig2<=0.0)
    dsig2 = 0.01;

  ninput = som->ninput;
  w = som->weight;
  xj[nd-1] = -1;
  for (j=som->nneur, n=0; j--; n++)
    {
/*-- Update the current SOM-position vector */
    (*xj)++;
    xjt = xj;
    xwint = xwin;
    dt = dist;
    mul=som->neurstep;
    for (ndt = nd; --ndt && !(n%(*(++mul))); (*xjt)++, xwint++, dt++)
      *(xjt++) = 0;

/*-- Update the current distance from the winning node */
    while(ndt<nd)
      {
      diff = (*(xjt--) - *(xwint--));
      *dt = diff*diff;
      if (ndt++)
        *dt += *(dt+1);
      dt--;
      }

/*-- Update weights according to the kernel amplitude */
    h = amp*SOM_KERNEL((float)*dist, dsig2);
    inputt = input;
    for (i=ninput; i--; w++)
      *w += h*(*(inputt++)-*w);
    }

  som->ntrain++;

  return;
  }


/********************************** som_win **********************************/
/*
Return the winner node number in the SOM.
1 degree of freedom is left: the amplitude of the prototype.
*/
int	som_win(somstruct *som, float *input, float *inputw, float *errp,
		float *amp)
  {
   double	sxx,sw,sxy,syy,erri,erre,err,errmin,errimin,erremin, b,bmin;
   float	*wt,*xt,*yt, xi,yi,wi,wxi, diff;
   int		i,j,n, nmin, nima,nextra;

  xt = som->weight;
  nextra = som->nextrainput;
  nima = som->ninput-nextra;
  errmin = erremin = errimin = BIG;
  b = bmin = 0.0;

  nmin = -1;
  n=0;
  for (j=som->nneur; j--; n++)
    {  
/*-- First, the error from the image-fitting */
    yt = input;
    wt = inputw;
    if (amp)
      {
      sxx = sxy = syy = sw = 0.0;
      for (i=nima; i--;)
        {
        sxy += (wxi = (wi=*(wt++))*(xi=*(xt++)))*(yi=*(yt++));
        sxx += wxi*xi;
        syy += wi*yi*yi;
        sw += wi;
        }
      b = sxy/sxx;
      err = erri = nima*(b*b*sxx + syy - 2.0*b*sxy)/sw;
      }
    else
      {
      erri = 0.0;
      for (i=nima; i--;)
        {
        diff = *(yt++) - *(xt++);
        erri += diff*diff*(double)*(wt++);
        }
      err = (erri /= (double)nima);
      }

/*-- Second, the error of non-pixel parameters */
    erre = 0.0;
    for (i=nextra; i--;)
      {
      diff = *(yt++) - *(xt++);
      erre += (diff*diff*(double)*(wt++))/(double)nextra;
      }

    err += erre;
    if (err<errmin)
      {
      nmin = n;
      errmin = err;
      errimin = erri;
      erremin = erre;
      bmin = b;
      }
    }

  if (errimin<0.0)
    errimin = 0.0;
  som->errpix += errimin;
  som->errextra += erremin;
  som->err += errmin;

  if (errp)
    *errp = (float)sqrt(errimin);

  if (amp)
    *amp = (float)bmin;

  return nmin;
  }


/********************************* som_end ***********************************/
/*
Terminate SOM.
*/
void	som_end(somstruct *som)
  {
/* Free memory*/
  free(som->weight);
  free(som->residuals);
  free(som->freq);
  free(som->stderror);
  free(som->amp);
  free(som->sigamp);
  free(som->inputsize);
  free(som->neursize);
  free(som->neurstep);
  free(som);

/* locals */
  free(sx);
  free(sxx);
  free(sxy);
  free(yweight);
  return;
  }


/********************************* som_load **********************************/
/*
Read the SOM weights in a FITS file.
*/
somstruct	*som_load(char *filename)
  {
   somstruct	*som;
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   char		*head, str[80];
   int		i;

/* Open the cat (well it is not a "cat", but simply a FITS file */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: SOM file not found: ", filename);

  if (!(tab = name_to_tab(cat, "SOM", 0)))
    error(EXIT_FAILURE, "*Error*: SOM table not found in catalog ",
	filename);

/* OK, we now allocate memory for the SOM structure itself */
  QCALLOC(som, somstruct, 1);

/* Load important scalars (which are stored as FITS keywords) */
  head = tab->headbuf;

/* Dimensionality of the input */
  if (fitsread(head, "INPNAXIS", &som->inputdim, H_INT, T_LONG) != RETURN_OK)
    goto headerror;
  if (som->inputdim>INPUT_MAXDIM)
    {
    sprintf(str, "%d", INPUT_MAXDIM);
    error(EXIT_FAILURE, "*Error*: This package is presently limited to inputs"
	"with dimensionality less or equal to ", str);
    }
  QMALLOC(som->inputsize, int, INPUT_MAXDIM);
  for (i=0; i<INPUT_MAXDIM; i++)
    som->inputsize[i] = 1;
  som->ninput = 1;
  for (i=0; i<som->inputdim; i++)
    {
    sprintf(str, "INPAXIS%1d", i+1);
    if (fitsread(head, str, &som->inputsize[i], H_INT,T_LONG) != RETURN_OK)
      goto headerror;
    som->ninput *= som->inputsize[i];
    }
    if (fitsread(head,"INPNEXTR",&som->nextrainput,H_INT,T_LONG) != RETURN_OK)
      som->nextrainput = 0;
    som->ninput += som->nextrainput;

/* Dimensionality of the SOM */
  if (fitsread(head, "SOMNAXIS", &som->neurdim, H_INT, T_LONG) != RETURN_OK)
    goto headerror;
  QMALLOC(som->neursize, int, SOM_MAXDIM);
  QMALLOC(som->neurstep, int, SOM_MAXDIM+1);
  for (i=0; i<SOM_MAXDIM; i++)
    som->neursize[i] = 1;
  som->neurstep[0] = som->nneur = 1;
  for (i=0; i<som->neurdim; i++)
    {
    sprintf(str, "SOMAXIS%1d", i+1);
    if (fitsread(head, str, &som->neursize[i], H_INT,T_LONG) != RETURN_OK)
      goto headerror;
    som->nneur *= som->neursize[i];
    som->neurstep[i+1] = som->nneur;
    }

/* Other scalars */
  if (fitsread(head, "SOMLRATE", &som->learnrate,H_FLOAT,T_FLOAT) != RETURN_OK)
    goto headerror;
  som->clearnrate = som->learnrate;
  if (fitsread(head, "SOMKERNW", &som->kernw, H_FLOAT,T_FLOAT) != RETURN_OK)
    goto headerror;
  som->ckernw = som->kernw;
  if (fitsread(head, "SOMNPASS", &som->ntrain , H_INT, T_LONG) != RETURN_OK)
    goto headerror;
  if (fitsread(head, "SOMNSWEE", &som->nsweep , H_INT, T_LONG) != RETURN_OK)
    goto headerror;

  som->nweight = som->nneur*som->ninput;
  QCALLOC(som->residuals, float, som->nweight);
  QCALLOC(som->freq, int, som->nneur);
  QMALLOC(som->stderror, float, som->nneur);
  QMALLOC(som->amp, float, som->nneur);
  QMALLOC(som->sigamp, float, som->nneur);
/* Locals */
  QMALLOC(sx, double, som->nneur);
  QMALLOC(sxx, double, som->nneur);
  QMALLOC(sxy, double, som->nneur);
  QMALLOC(yweight, float, som->ninput);

/* Load the weight vector */
  key = read_key(tab, "WEIGHTS");
  som->weight = key->ptr;

/* But don't touch my arrays!! */
  blank_keys(tab);
  free_cat(&cat, 1);

  return som;

headerror:
  error(EXIT_FAILURE, "*Error*: Incorrect or obsolete SOM data in ", filename);

  return NULL;
  }


/********************************* som_save **********************************/
/*
Save the SOM weights as a FITS file.
*/
void	som_save(somstruct *som, char *filename)
  {
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   char		*head, str[80];
   int		i;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  cat = new_cat(1);
  init_cat(cat);
  tab = new_tab("SOM");
  add_tab(tab, cat, 0);
/* Let's allocate more than strictly necessary to be sure... */
  QREALLOC(tab->headbuf, char, (3+tab->headnblock)*FBSIZE);
  head = tab->headbuf;
/* ... and blank the extra space */
  memset(head+FBSIZE*tab->headnblock, ' ', 3*FBSIZE);

/* Add and write important scalars as FITS keywords */
  fitsadd(head, "INPNAXIS", "Dimensionality of the input vector");
  fitswrite(head, "INPNAXIS", &som->inputdim, H_INT, T_LONG);
  for (i=0; i<som->inputdim; i++)
    {
    sprintf(str, "INPAXIS%1d", i+1);
    fitsadd(head, str, "Number of element along this axis");
    fitswrite(head, str, &som->inputsize[i], H_INT, T_LONG);
    }
  fitsadd(head, "INPNEXTR", "Number of extra parameters");
  fitswrite(head, "INPNEXTR", &som->nextrainput, H_INT, T_LONG);

  fitsadd(head, "SOMNAXIS", "Dimensionality of the SOM");
  fitswrite(head, "SOMNAXIS", &som->neurdim, H_INT, T_LONG);
  for (i=0; i<som->neurdim; i++)
    {
    sprintf(str, "SOMAXIS%1d", i+1);
    fitsadd(head, str, "Number of element along this axis");
    fitswrite(head, str, &som->neursize[i], H_INT, T_LONG);
    }

/* Other scalars */
  fitsadd(head, "SOMLRATE", "Current learning rate");
  fitswrite(head, "SOMLRATE", &som->clearnrate, H_FLOAT, T_FLOAT);
  fitsadd(head, "SOMKERNW", "Current kernel width");
  fitswrite(head, "SOMKERNW", &som->ckernw, H_FLOAT, T_FLOAT);
  fitsadd(head, "SOMNPASS","Number of training passes so far");
  fitswrite(head, "SOMNPASS", &som->ntrain, H_INT, T_LONG);
  fitsadd(head, "SOMNSWEE","# of sweeps through the data set so far");
  fitswrite(head, "SOMNSWEE", &som->nsweep, H_INT, T_LONG);
  fitsadd(head, "SOMSTIFF", "Stiffness of the X/Y mapping");
  fitswrite(head, "SOMSTIFF", &som->xy_stiff, H_FLOAT, T_FLOAT);

/* Create and fill the arrays */
  key = new_key("WEIGHTS");
  key->naxis = som->neurdim+1;
  QMALLOC(key->naxisn, int, key->naxis);
  key->naxisn[0] = som->ninput;
  for (i=0; i<som->neurdim; i++)
    key->naxisn[i+1] = som->neursize[i];
  strcat(key->comment, "Weight vector, neuron per neuron");
  key->htype = H_FLOAT;
  key->ttype = T_FLOAT;
  key->nbytes = som->nweight*t_size[T_FLOAT];
  key->nobj = 1;
  key->ptr = som->weight;
  add_key(key, tab, 0);

/* Find the useful FITS header size */
  tab->headnblock = fitsfind(head, "END     ")/(FBSIZE/80)+1;

/* Then, just save everything and free memory */
  save_cat(cat, filename);

/* But don't touch my arrays!! */
  blank_keys(tab);
  free_cat(&cat, 1);

  return;
  }


/******************************* som_writecheck ******************************/
/*
Write a FITS image for check.
*/
void	som_writecheck(somstruct *som, char *filename)
  {
   catstruct	*cat;
   tabstruct	*tab;
   char		*head;
   float	*pix, *w;
   int		i,j,x,y, iw, nw, nextra, nextraline;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  cat = new_cat(1);
  init_cat(cat);
  tab = cat->tab;
  tab->naxis = 2;	/* This is an image */
  QMALLOC(tab->naxisn, int, tab->naxis);
  fitsremove(tab->headbuf,"HISTORY ");
  fitsremove(tab->headbuf,"EXTEND  ");
  head = tab->headbuf;

  switch(prefs.check_type)
    {
    case SOM_PROTO:
/*---- Weight vectors for all neurons are arranged as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nextra = som->nextrainput;
      nextraline = (nextra+som->inputsize[0]-1)/som->inputsize[0];
      tab->naxisn[0] = som->neursize[0]*som->inputsize[0];
      tab->naxisn[1] = som->neursize[2]*som->neursize[1]
		*(som->inputsize[1]+nextraline);
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix, float, tab->tabsize);
      iw = som->inputsize[0];
      nw = som->ninput*som->neursize[0];
      tab->bodybuf = (char *)pix;
      for (j=0; j<(som->neursize[1]*som->neursize[2]); j++)
        {
        w = som->weight + j*nw;
        for (y=som->inputsize[1]; y--; w += iw - nw)
          for (i=som->neursize[0]; i--; w += som->ninput - iw)
            for (x=iw; x--;)
              *(pix++) = *(w++);
        if (nextra)
          for (i=som->neursize[0]; i--; w += som->ninput - nextra)
            {
            for (x=nextra; x--;)
              *(pix++) = *(w++);
            pix += iw-nextra;
            }
        }
      break;

    case SOM_HIT:
/*---- One pixel per neuron, with the number of hit */
      tab->bitpix =  BP_LONG;
      tab->bytepix = t_size[T_LONG];
      tab->naxisn[0] = som->neursize[0];
      tab->naxisn[1] = som->neursize[2]*som->neursize[1];
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QMALLOC(tab->bodybuf, char, tab->tabsize*sizeof(int));
      memcpy(tab->bodybuf, som->freq, tab->tabsize);
      break;

    case SOM_RESIDUALS:
/*---- Residual vectors for all neurons are arranged as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nextra = som->nextrainput;
      nextraline = (nextra+som->inputsize[0]-1)/som->inputsize[0];
      tab->naxisn[0] = som->neursize[0]*som->inputsize[0];
      tab->naxisn[1] = som->neursize[2]*som->neursize[1]
		*(som->inputsize[1]+nextraline);
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix, float, tab->tabsize);
      iw = som->inputsize[0];
      nw = som->ninput*som->neursize[0];
      tab->bodybuf = (char *)pix;
      for (j=0; j<(som->neursize[1]*som->neursize[2]); j++)
        {
        w = som->residuals + j*nw;
        for (y=som->inputsize[1]; y--; w += iw - nw)
          for (i=som->neursize[0]; i--; w += som->ninput - iw)
            for (x=iw; x--;)
              *(pix++) = *(w++);
        if (nextra)
          for (i=som->neursize[0]; i--; w += som->ninput - nextra)
            {
            for (x=nextra; x--;)
              *(pix++) = *(w++);
            pix += iw-nextra;
            }
        }
      break;

    default:
      error(EXIT_FAILURE, "*Internal Error*: Yet unavailable CHECKIMAGE type",
	"");
    }

/* Then, just save everything and free memory */
  save_cat(cat, filename);
  free_cat(&cat, 1);

  return;
  }

