 /*
 				main.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMFit
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Parsing of the command line.
*
*	Last modify:	15/03/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_PLPLOT
#include <plplot.h>
#endif

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"

#define		SYNTAX \
"somfit catalog1 [catalog2,...][-c <config_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file: " BANNER " -d \n" \
"> to dump a default extended configuration file: " BANNER " -dd \n"

extern const char	notokstr[];

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   static char	prefsname[MAXCHAR];
   char		**argkey, **argval, *str;
   int		a, narg, nim, opt,opt2;

  if (argc<2)
    {
    fprintf(OUTPUT, "\n                %s  Version %s (%s)\n",
		BANNER, MYVERSION, DATE);
    fprintf(OUTPUT, "\nFor information, please contact: %s\n", COPYRIGHT);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

#ifdef HAVE_PLPLOT
  if (argc>2)
    plParseOpts(&argc, argv, PL_PARSE_SKIP);
#endif

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/* Default parameters */
  prefs.nfile = 1;
  prefs.file_name[0] = "catalog";
  strcpy(prefsname, "somfit.conf");
  narg = nim = 0;

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefsname, argv[++a]);
            break;
          case 'd':
            dumpprefs(opt2=='d' ? 1 : 0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(EXIT_SUCCESS);
            break;
          case 'h':
            fprintf(OUTPUT, "\nSYNTAX: %s", SYNTAX);
#ifdef HAVE_PLPLOT
            fprintf(OUTPUT, "\nPLPLOT-specific options:\n");
            plParseOpts(&argc, argv, PL_PARSE_SKIP);
#endif
            exit(EXIT_SUCCESS);
            break;
          default:
            error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
          }
        }
      else
        {
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      {
/*---- The input image filename(s) */
      for(; (a<argc) && (*argv[a]!='-'); a++)
        for (str=NULL; (str=strtok(str?NULL:argv[a], notokstr)); nim++)
          if (nim<MAXFILE)
            prefs.file_name[nim] = str;
          else
            error(EXIT_FAILURE, "*Error*: Too many input catalogs: ", str);
      prefs.nfile = nim;
      a--;
      }
    }

  readprefs(prefsname, argkey, argval, narg);
  useprefs();

  free(argkey);
  free(argval);

  makeit();

  endprefs();

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "> All done (in %d s)\n", prefs.time_diff);

  exit(EXIT_SUCCESS);
  }

