/* pp.c - Print parameters */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void pp_init(gptinit *init)
{
  int i ;
  
  for(i=1 ; i<=gptgetargnum(init) ; i++)
  {
    switch( gptgetargtype(init,i) )
    {
      case GPTTYPE_INT:
        fprintf( stderr, "%d ", gptgetargint(init,i) ) ;
        break ;

      case GPTTYPE_DOUBLE:
        fprintf( stderr, "%.16g ", gptgetargdouble(init,i) ) ;
        break ;

      case GPTTYPE_STRING:
        fprintf( stderr, "%s ", gptgetargstring(init,i) ) ;
        break ;
    }
  }
  fprintf( stderr, "\n" ) ;
}
