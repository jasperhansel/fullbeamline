/* outvalue.c - Wirte a single value in the outputfile */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void outputvalue_init(gptinit *init)
{
  char *name ;
  double value ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(name,value)\n", gptgetname(init) ) ;

  name  = gptgetargstring(init,1) ;
  value = gptgetargdouble(init,2) ;

  gptoutputdouble( name, value ) ;
}
