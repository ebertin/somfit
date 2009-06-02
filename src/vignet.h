 /*
 				vignet.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	A program using Self-Organizing Maps (SOMs)
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for vignet.c.
*
*	Last modify:	14/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _VIGNET_H_
#define _VIGNET_H_

/*--------------------------------- constants -------------------------------*/
/*---------------------------------- macros ---------------------------------*/
/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/

extern int		copyvignet_center(float *source, int sw, int sh,
				float *dest, int w, int h, float dx, float dy,
				float fwhm, float dxc, float dyc);

#endif

