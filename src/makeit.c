 /*
 				makesomcat.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMFit
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Learning and testing.
*
*	Last modify:	14/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "som.h"
#include "vignet.h"

#define		PHOTOM_EPS	0.005	/* Typical PSF relative accuracy */
#define		LSTACK_DEFSIZE	200	/* default size for the learn-stack */

/******************************** makesomcat *********************************/
/*
*/
void	makeit(void)

  {
   somstruct		*som, *som2;
   catstruct		*incat;
   tabstruct		*tab;
   keystruct		*key;
   static char		str[MAXCHAR];
   char			*head;
   unsigned short	*flags;
   double		dnorm;
   float		*vigstack,*vig,*vigt,
			*inputstack,*input,*inputt,
			*wstack,*w,*wt, *gwstack, *gw,
			*inputswap, *wswap, gwswap, *r,
			*xm, *ym, *vignet, *peakflux, *flux,
			weight, pix, norm,norm2, noise, stiffx2, stiffy2,
			backnoise,backnoise2,gain, errextra, errpix,
			dx,dy, xstep,ystep;
   int			*f,
			c, i,j, n, nstack,nstackmax, imaw, imah,
			retinaw,retinah,retinasize, inputsize, 
			vigw, vigh, vigsize, nobj, nextra, niter, offset,
			wsampsize,hsampsize, ampflag;


  nextra = 0;
/*-------------------------- Avoid gcc -Wall warnings -----------------------*/
  stiffx2 = stiffy2 = 0.0;	
  vigsize = vigw = vigh = retinasize = retinaw = retinah = inputsize
	= nstackmax = nstack = 0;
  vigstack = inputstack = wstack= gwstack = gw = w = input = vig = NULL;
  if (prefs.spatial_mapflag)
    {
    nextra += 2;
    stiffx2 = prefs.spatial_stiffness[0]*prefs.spatial_stiffness[0];
    stiffy2 = prefs.spatial_stiffness[1]*prefs.spatial_stiffness[1];
    }
   ampflag = prefs.amp_flag;

  for (c=0; c<prefs.nfile; c++)
    {
/*-- Read input catalog */
    if (!(incat = read_cat(prefs.file_name[c])))
      error(EXIT_FAILURE, "*Error*: No such catalog: ", prefs.file_name[c]);

    if ((tab = name_to_tab(incat, "LDAC_IMHEAD", 0))
	&& (key=read_key(tab, "Field Header Card")))
      head = key->ptr;
    else if (!(head=incat->tab->headbuf))
      error(EXIT_FAILURE, "*Error*: I expected a FITS catalog format", "");

    if (fitsread(head, "SEXBKDEV", &backnoise,H_FLOAT,T_FLOAT)== RETURN_ERROR)
      error(EXIT_FAILURE, "*Error*: Keyword not found:", "SEXBKDEV");
    backnoise2 = backnoise*backnoise;

    if (fitsread(head, "SEXGAIN", &gain, H_FLOAT, T_FLOAT) == RETURN_ERROR)
      error(EXIT_FAILURE, "*Error*: Keyword not found:", "SEXGAIN");

    if (prefs.spatial_mapflag)
      {
      if (fitsread(head,  "SEXIMASX", &imaw, H_INT, T_LONG)== RETURN_ERROR
	&& fitsread(head, "NAXIS1  ", &imaw, H_INT, T_LONG)== RETURN_ERROR)
        error(EXIT_FAILURE,"*Error*: Image X-size not found","");

      if (fitsread(head,  "SEXIMASY", &imah, H_INT, T_LONG)== RETURN_ERROR
	&& fitsread(head, "NAXIS2  ", &imah, H_INT, T_LONG)== RETURN_ERROR)
        error(EXIT_FAILURE,"*Error*: Image Y-size not found","");
      }

    if (!(tab = name_to_tab(incat, "LDAC_OBJECTS", 0))
	&& !(tab = name_to_tab(incat, "OBJECTS", 0)))
      error(EXIT_FAILURE, "*Error*: OBJECTS table not found in catalog ",
		prefs.file_name[c]);

    if (!(key = read_key(tab, "X_IMAGE")))
      error(EXIT_FAILURE, "*Error*: X_IMAGE parameter not found in catalog ",
		prefs.file_name[c]);
    xm = (float *)key->ptr;
    if (!(key = read_key(tab, "Y_IMAGE")))
      error(EXIT_FAILURE, "*Error*: Y_IMAGE parameter not found in catalog ",
		prefs.file_name[c]);
    ym = (float *)key->ptr;

    if (!(key = read_key(tab, "FLUX_APER")))
      error(EXIT_FAILURE, "*Error*: FLUX_APER parameter not found in catalog ",
		prefs.file_name[c]);
    flux = (float *)key->ptr;

    if (!(key = read_key(tab, "FLUX_MAX")))
      error(EXIT_FAILURE, "*Error*: FLUX_MAX parameter not found in catalog ",
		prefs.file_name[c]);
    peakflux = (float *)key->ptr;

    if (!(key = read_key(tab, "FLAGS")))
      error(EXIT_FAILURE, "*Error*: FLAGS parameter not found in catalog ",
		prefs.file_name[c]);
    flags = (unsigned short *)key->ptr;
    nobj = key->nobj;

    if (!(key = read_key(tab, "VIGNET")))
      error(EXIT_FAILURE,
	"*Error*: VIGNET parameter not found in catalog ", prefs.file_name[c]);
    vignet = (float *)key->ptr;
    if (key->naxis != 2)
      error(EXIT_FAILURE, "*Error*: VIGNET should be a 2D vector", "");

/* Allocate memory for the first shipment */
    if (!c)
      {
      vigw = *(key->naxisn);
      vigh = *(key->naxisn+1);
      vigsize = vigw*vigh;
/*---- The retina cannot be larger than the input images! */
      if (prefs.retisize[0]>vigw)
        prefs.retisize[0] = vigw;
      if (prefs.retisize[1]>vigh)
        prefs.retisize[1] = vigh;
/*---- Set up the SOM structure */
      retinaw = prefs.retisize[0];
      retinah = prefs.retisize[1];
      retinasize = retinaw*retinah;
      inputsize = retinasize + nextra;
      nstackmax = LSTACK_DEFSIZE;
      QMALLOC(vigstack, float, nstackmax*vigsize);
      QMALLOC(inputstack, float, nstackmax*inputsize);
      QMALLOC(wstack, float, nstackmax*inputsize);
      QMALLOC(gwstack, float, nstackmax);
      nstack = 0;
      vig = vigstack;
      input = inputstack;
      w = wstack;
      gw = gwstack;
      }
    else
      {
      if (vigw != *(key->naxisn) || vigh != *(key->naxisn+1))
        error(EXIT_FAILURE, "*Error*: Incompatible VIGNET sizes in ",
		prefs.file_name[c]);
      }

/* Now examine each vector of the shipment */
    for (n=nobj; n--; vignet += vigsize, xm++,ym++,flags++,peakflux++,flux++)
      {
      if (!(n%100))
        {
        sprintf(str,"Object countdown #%d, %d samples stored",n,nstack+1);
        NFPRINTF(OUTPUT, str);
        }
      if (*flags<2 && *flux>prefs.minflux)
        {
/*------ Increase storage space to receive new candidates if needed */
        if (nstack>=nstackmax)
          {
          nstackmax = (int)(1.62*nstackmax);
          QREALLOC(vigstack, float, nstackmax*vigsize);
          QREALLOC(inputstack, float, nstackmax*inputsize);
          QREALLOC(wstack, float, nstackmax*inputsize);
          QREALLOC(gwstack, float, nstackmax);
          vig = vigstack+nstack*vigsize;
          input = inputstack+nstack*inputsize;
          w = wstack+nstack*inputsize;
          gw = gwstack+nstack;
          }

        memcpy(vig, vignet, vigsize*sizeof(float));
        dx = *xm-(int)(*xm+0.49999);
        dy = *ym-(int)(*ym+0.49999);
        vigt = vig;
        for (i=vigsize; i--;)
          if (*(vigt++) < -1e10)
            *(vigt-1) = 0.0;
        copyvignet_center(vig, vigw, vigh, input, retinaw, retinah,
	  dx, dy, prefs.fwhm, dx,dy);
/*------ Copy and normalisation of the useful part for classification */
        *gw = *flux;
        norm = *flux;
        norm2 = norm*norm;
        inputt = input;
        wt = w;
        for (i=retinasize; i--;)
          {
          pix = ampflag? *(inputt++) : (*(inputt++) /= norm);
          noise = backnoise2;
          if (pix>0.0 && gain>0.0)
            noise += pix/gain;
          *(wt++) = ampflag? 1.0/noise : norm2/noise;
          }

/*------ Constraints from the position of the detection */
        if (prefs.spatial_mapflag)
          {
          *(inputt++) = *xm/imaw;
          *(wt++) = 1/stiffx2;
          *(inputt++) = *ym/imah;
          *(wt++) = 1/stiffy2;
          }

        vig += vigsize;
        input += inputsize;
        w += inputsize;
        gw++;
        nstack++;
        }
      }

    free_cat(&incat, 1); 
    }

  if (!nstack)
    error(EXIT_FAILURE, "*Error*: No source selected!!","");

/* Don't spoil memory! */
  if (nstack<nstackmax)
    {
    QREALLOC(vigstack, float, nstack*vigsize);
    QREALLOC(inputstack, float, nstack*inputsize);
    QREALLOC(wstack, float, nstack*inputsize);
    QREALLOC(gwstack, float, nstack);
    }

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "%d samples loaded\n", nstack);

/* Filter sources belonging to classes with not enough members */
  if (prefs.filter_flag)
    {
     float	thresh;
     int	*win, *wint, freqsum;

/*-------------------------- Avoid gcc -Wall warnings -----------------------*/
  freqsum = 0;	

    NFPRINTF(OUTPUT,"Loading the filtering SOM...");
    som = som_load(prefs.filter_name);

    NFPRINTF(OUTPUT,"Filtering vectors...");
/*-- Make the freq histogram */
    f = som->freq;
    memset(f, 0, som->nneur*sizeof(int));
    QMALLOC(win, int, nstack);
    wint = win;
    input = inputstack;
    w = wstack;
    for (i=nstack; i--; input+=nextra, w+=inputsize)
      if (ampflag)
        f[*(wint++) = som_win(som, input, w, NULL, NULL)]++;
      else
        {
        f[*(wint++) = som_win(som, input, w, NULL, &norm)]++;
        for (j=retinasize; j--;)
          *(input++) /= norm;
        }

    thresh = (float)prefs.filter_thresh;
/*-- If thresh is a percentage, convert it to an absolute threshold */
    if (thresh<1.0)
      {
      for (i=som->nneur; i--;)
        freqsum += *(f++);
      thresh *= (float)freqsum;
      }

    n = 0;
    f = som->freq;
    vig = vigstack;
    input = inputstack;
    w = wstack;
    gw = gwstack;
    wint = win;
    for (i=0; i<nstack; i++)
      if ((float)f[*(wint++)]>= thresh)
        {
        if ((n++) != i)
          {
          memcpy(vig, vigstack+i*vigsize, vigsize*sizeof(float));
          memcpy(input, inputstack+i*inputsize, inputsize*sizeof(float));
          memcpy(w, wstack+i*inputsize, inputsize*sizeof(float));
          }
        *(gw++) = gwstack[i];
        vig += vigsize;
        input += inputsize;
        w +=inputsize;
        }

  if (!n)
    error(EXIT_FAILURE, "*Error*: No pattern left after filtering!", "");

  sprintf(str, "%d patterns removed\n", nstack - n);
  NFPRINTF(OUTPUT, str);

/*-- Don't spoil memory! */
    free(win);
    if (n<nstack)
      {
      nstack = n;
      QREALLOC(vigstack, float, nstack*vigsize);
      QREALLOC(inputstack, float, nstack*inputsize);
      QREALLOC(wstack, float, nstack*inputsize);
      QREALLOC(gwstack, float, nstack);
      }

    som_end(som);
    }

/* Let's mix the data to avoid systematic effects */
/*
  NFPRINTF(OUTPUT,"Mixing vectors...");
  QMALLOC(inputswap, float, inputsize);
  QMALLOC(wswap, float, inputsize);
  input = inputstack;
  w = wstack;
  gw = gwstack;
  for (n=nstack; n--;  input+=inputsize, w+=inputsize, gw++)
    {
    i=(int)(rand()/(RAND_MAX+1.0)*nstack);
    inputt = inputstack+i*inputsize;
    wt = wstack+i*inputsize;
    memcpy(inputswap, inputt, inputsize*sizeof(float));
    memcpy(inputt, input, inputsize*sizeof(float));
    memcpy(input, inputswap, inputsize*sizeof(float));
    memcpy(wswap, wt, inputsize*sizeof(float));
    memcpy(wt, w, inputsize*sizeof(float));
    memcpy(w, wswap, inputsize*sizeof(float));
    gwswap = gwstack[i];
    gwstack[i] = *gw;
    *gw = gwswap;
    }

  free(inputswap);
  free(wswap);
*/


/* Now init the real SOM */
  som = som_init(prefs.retisize, 2, nextra);

/* Initialize SOM weights to a noisy "average" value */
  memset(som->weight, 0, som->nweight*sizeof(float));
  w = som->weight;
  for (j=som->nneur; j--; w+= inputsize)
    {
    input = inputstack;
    for (n=nstack; n--;)
      {
      wt = w;
      for (i=inputsize; i--;)
        *(wt++) += *(input++);
      }
    wt = w;

    for (i=inputsize; i--; wt++)
      *wt = *wt/nstack+(0.5-rand()/(RAND_MAX+1.0))/inputsize;
    }

/* The learning itself */
  niter = prefs.niter/nstack;
  if (prefs.niter && !niter)
    niter = 1;
  errpix = errextra = 0.0;
  for (i=niter; i--;)
    {
    input = inputstack;
    w = wstack;
    gw = gwstack;
    som->err = som->errpix = som->errextra = 0.0;
    NFPRINTF(OUTPUT, str);
    for (n=nstack; n--; input += inputsize, w += inputsize)
      som_learn(som, input, w, 1.0, ampflag);
    errpix = sqrt(som->errpix/(nstack*retinasize));
    errextra = nextra?sqrt(som->errextra/(nstack*nextra)):0.0;
    sprintf(str,"Learning countdown: %5d / Pixel RMS: %7.3f"
	" / Extras RMS: %7.3f",
		i, errpix, errextra);
    }

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT,"Centered map reduced errors"
	" -- Pixel: %7.3f / Extras: %7.3f\n",
	errpix, errextra);

/* Make the freq histogram and the residuals map */
  f = som->freq;
  memset(f, 0, som->nneur*sizeof(int));
  input = inputstack;
  w = wstack;
  for (i=nstack; i--; w+=inputsize)
    {
    f[n=som_win(som, input, w, NULL, ampflag?NULL:&norm)]++;
    wt = som->weight + n*inputsize;
    r = som->residuals+n*inputsize;
    for (j=retinasize; j--;)
      {
      errpix = *(input++) - norm**(wt++);
      *(r++) += errpix*errpix;
      }
    for (j=nextra; j--;)
      {
      errextra = *(input++) - *(wt++);
      *(r++) += errextra*errextra;
      }
    }

/* Normalise residuals */
  r = som->residuals;
  f = som->freq;
  for (j=som->nneur; j--; f++)
    for (i=retinasize; i--; r++)
      *r = sqrt(*r/((float)(*f>1?*f:1)));

/* Prepare the final map with shifted profiles */

/* What size is requested to have a PSF accuracy of PHOTOM_EPS? */
  hsampsize = (((int)(2.355/(prefs.fwhm*sqrt(8*PHOTOM_EPS))+2.0))/2)*2;
  if (hsampsize<1)
    hsampsize = 1;
  wsampsize = hsampsize;
  xstep = 1.0/wsampsize;
  ystep = 1.0/hsampsize;
  prefs.somsize[prefs.nsomsize++] = wsampsize;
  prefs.somsize[prefs.nsomsize++] = hsampsize;
  NPRINTF(OUTPUT, "Subsampling: %d x %d prototypes\n", wsampsize, hsampsize);

/* Init the last final SOM */
/*
  som2 = som_init(prefs.retisize, 2, nextra);
  wt = som2->weight;
  for (j=hsampsize, dy=-0.5; j--; dy+=ystep)
    for (i=wsampsize, dx=-0.5; i--; dx+=xstep)
      {
      w = som->weight;
      for (n=som->nneur; n--;)
        {
        copyvignet_center(w, retinaw, retinah, wt, retinaw, retinah, -dx,-dy,
		prefs.fwhm, 0.0,0.0);
        w += retinasize;
        wt += retinasize;
        for (c=nextra; c--;)
          *(wt++) = *(w++);
        }
      }

  NFPRINTF(OUTPUT, "Testing shifted map...");
*/
/* Make the freq histogram and the residuals map for the shifted SOM */
/*
  f = som2->freq;
  memset(f, 0, som2->nneur*sizeof(int));
  memset(som2->residuals, 0, som2->nneur*sizeof(int));
  som2->err = som2->errpix = som2->errextra = 0.0;
  vig = vigstack;
  input = inputstack;
  w = wstack;
  gw = gwstack;
  for (i=nstack; i--; vig+=vigsize, w += inputsize, input += inputsize)
    {
    copyvignet_center(vig, vigw, vigh, input, retinaw, retinah, 0.0,0.0,
	0.0,0.0,0.0);
*/
/*-- Copy and normalisation of the useful part for classification */
/*
    norm = *(gw++);
    norm2 = norm*norm;
    inputt = input;
    wt = w;
    for (j=retinasize; j--;)
      {
      pix = ampflag? *(inputt++) : (*(inputt++) /= norm);
      noise = backnoise2;
      if (pix>0.0 && gain>0.0)
        noise += pix/gain;
      *(wt++) = ampflag? 1.0/noise : norm2/noise;
      }

    f[n=som_win(som2, input, w, NULL, ampflag?NULL:&norm)]++;
    }

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT,"Shifted map reduced errors -- Pixel: %7.3f / Extras: %7.3f\n",
	sqrt(som2->errpix/(nstack*retinasize)),
	nextra?sqrt(som2->errextra/(nstack*nextra)):0.0);
*/
/* Free memory */
  free(vigstack);
  free(inputstack);
  free(wstack);
  free(gwstack);

  prefs.check_type = SOM_PROTO;
  sprintf(prefs.check_name, "proto.fits");
  som_writecheck(som, prefs.check_name);
  prefs.check_type = SOM_PROTO;
/*
  sprintf(prefs.check_name, "proto2.fits");
  som_writecheck(som2, prefs.check_name);
*/
  prefs.check_type = SOM_HIT;
  sprintf(prefs.check_name, "hit.fits");
  som_writecheck(som, prefs.check_name);

  prefs.check_type = SOM_RESIDUALS;
  sprintf(prefs.check_name, "resi.fits");
  som_writecheck(som, prefs.check_name);

  if (!prefs.filter_flag)
    som_save(som, prefs.filter_name);
/*
  som_save(som2, prefs.nnw_name);
*/
  som_end(som);
/*
  som_end(som2);
*/
  NFPRINTF(OUTPUT,"Done.\n");

  return;
  }

