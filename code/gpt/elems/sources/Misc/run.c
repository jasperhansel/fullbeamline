/* run.c - Run external command */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "elem.h"

#define RUNSTRINGLENGTH 1024

void run_init(gptinit *init)
{
  static char runstring[RUNSTRINGLENGTH] ;
  int rc, i ;
  char *p ;  

  /* Build runstring */
  runstring[0] = '\0' ;
  for(i=1 ; i<=gptgetargnum(init) ; i++)
  {
    p = &runstring[strlen(runstring)] ;
    if(p-runstring > RUNSTRINGLENGTH-256)
      gpterror( "String too long\n" ) ;

    switch( gptgetargtype(init,i) )
    {
      case GPTTYPE_INT:
        sprintf( p, "%d ", gptgetargint(init,i) ) ;
        break ;

      case GPTTYPE_DOUBLE:
        sprintf( p, "%.16g ", gptgetargdouble(init,i) ) ;
        break ;

      case GPTTYPE_STRING:
        sprintf( p, "%s ", gptgetargstring(init,i) ) ;
        break ;
    }
  }

  /* Execute runstring using system */
  rc=system( runstring ) ;

  /* Build application name */
  p = runstring ;
  while( isalnum(*p) ) p++ ;
  *p = '\0' ;

  /* Dislay error codes */
  switch( rc )
  {
    case  0: break ;
    case -1: gpterror( "Run error on \"%s\": %s\n", runstring, strerror(errno) ) ;
    default: gptwarning( "Application \"%s\" returned error %d\n", runstring, rc ) ;
  }
}
