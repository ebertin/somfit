 /*
 				vignet.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMFit
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Function related to vignet manipulations.
*
*	Last modify:	15/03/2006
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
#include "vignet.h"

#define	INTERPW		6	/* Interpolation function range (x) */
#define	INTERPH		6	/* Interpolation function range (y) */

#define	INTERPF(x)	(x==0.0?1.0:sin(PI*x)*sin(PI*x/3.0)/(PI*PI/3.0*x*x))
				/* Lanczos approximation */

static float	interpm[INTERPW*INTERPH];

/**************************** copyvignet_center ******************************/
/*
Copy a small part of the image and recenter it through sinc interpolation.
Image parts which lie outside boundaries are set to 0.
*/
int	copyvignet_center(float *source, int sw, int sh,
		float *dest, int w, int h, float dx, float dy,
		float fwhm, float dxc, float dyc)
  {
   PIXTYPE	*s,*s0, *dt,*dt0,*dt2;
   double	sxx,sxy;
   float	*m,
		ddx0,ddx,ddy0,ddy,ddy2,sum, f, mval, scale, amp;
   int		i,ix,iy, idmx,idmy, mx,my, xmin,ymin,xmin2,x0,y0,y2, w2,h2,
		idx,idy;

  
/* Initialize destination buffer to zero */
  memset(dest, 0, w*h*sizeof(float));

/* Compute the interpolation mask */
  ddx0 = -(idmx=(INTERPW-(dx>0.0?1:0))/2)-dx;
  ddy = -(idmy=(INTERPH-(dy>0.0?1:0))/2)-dy;
  sum = 0.0;
  m = interpm;
  for (my=INTERPH; my--; ddy+=1.0)
    {
    ddx = ddx0;
    f = INTERPF(ddy);
    for (mx=INTERPW; mx--;ddx+=1.0)
      sum += *(m++) = f*INTERPF(ddx);
    }  

/* Normalize it */

  m = interpm;
  for (i=INTERPW*INTERPH; i--;)
    *(m++) /= sum;

  amp = scale = ddy0 = 0.0;	/* To avoid gcc -Wall warnings */
  if (fwhm>0.01)
    {
/*-- Fit the amplitude of the model ... */
    sxx = sxy = 0.0;
    ddx0 = -dxc - (int)(sw/2);
    ddy0 = ddy = -dyc - (int)(sh/2);
    scale = 2.773/(fwhm*fwhm);
    s = source;
    for (iy=sh; iy--; ddy+=1.0)
      {
      ddx = ddx0;
      ddy2 = ddy*ddy;
      for (ix=sw; ix--; ddx+=1.0)
        {
        f = scale*(ddx*ddx+ddy2);
        f = f<70.0?exp(-f):0.0;
        sxx += f*f;
        sxy += f**(s++);
        }
      } 
    amp = sxy/sxx;
/*-- ... an subtract it from the input vignet */
    ddy = ddy0;
    s = source;
    for (iy=sh; iy--; ddy+=1.0)
      {
      ddx = ddx0;
      ddy2 = ddy*ddy;
      for (ix=sw; ix--; ddx+=1.0)
        {
        f = scale*(ddx*ddx+ddy2);
        f = f<70.0?amp*exp(-f):0.0;
        *(s++) -= f;
        }
      }     
    }

/* Do the interpolation */
  m = interpm;
  xmin = sw/2 - w/2 - idmx;
  ymin = sh/2 - h/2 - idmy;
  for (my=INTERPH; my--; ymin++)
    {
/*-- Set the image boundaries in y */
    if (ymin < 0)
      {
      dt0 = dest - w*ymin;
      y0 = 0;
      if ((h2 = h + ymin) < 0)
        h2 = 0;
      }
    else
      {
      dt0 = dest;
      y0 = ymin;
      h2 = h;
      }
    if ((idy = sh - y0) < h2)
      h2 = idy;
    xmin2 = xmin;
    for (mx=INTERPW; mx--; xmin2++)
      {
      mval = *(m++);
/*---- Set the image boundaries in x */
      if (xmin2 < 0)
        {
        dt = dt0 - xmin2;
        x0 = 0;
        if ((w2 = w + xmin2) < 0)
          w2 = 0;
        }
      else
        {
        dt = dt0;
        x0 = xmin2;
        w2 = w;
        }
      if ((idx = sw - x0) < w2)
        w2 = idx;

      if (h2 >= 0 && w2 >= 0)
        {
        s0 = source + x0;
        y2 = y0;
        for (iy=h2; iy--; dt+=w)
          {
          dt2 = dt;
          s = s0+sw*(y2++);
          for (ix=w2; ix--;)
            *(dt2++) += mval**(s++);
          }
        }
      }
    }

  if (fwhm>0.01)
    {
/*-- Add back the shifted model ... */
    ddy = ddy0;
    s = source;
    for (iy=sh; iy--; ddy+=1.0)
      {
      ddx = ddx0;
      ddy2 = ddy*ddy;
      for (ix=sw; ix--; ddx+=1.0)
        {
        f = scale*(ddx*ddx+ddy2);
        f = f<70.0?amp*exp(-f):0.0;
        *(s++) += f;
        }
      }

/*-- Shift the model center */
    ddx0 += (int)(sw/2 - w/2) + dx;
    ddy0 += (int)(sh/2 - h/2) + dy;

/*-- Add back the shifted model ... */
    ddy = ddy0;
    dt = dest;
    for (iy=h; iy--; ddy+=1.0)
      {
      ddx = ddx0;
      ddy2 = ddy*ddy;
      for (ix=w; ix--; ddx+=1.0)
        {
        f = scale*(ddx*ddx+ddy2);
        f = f<70.0?amp*exp(-f):0.0;
        *(dt++) += f;
        }
      }     
    }


  return RETURN_OK;
  }


