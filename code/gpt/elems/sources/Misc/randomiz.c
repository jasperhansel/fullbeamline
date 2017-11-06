/* randomiz.c: Initialize random generator */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "elem.h"

void randomize_init(gptinit *init)
{
  int numarg, seed ;

  /* Print usage line when the number of parameters is incorrect */
  numarg = gptgetargnum(init) ;
  if( numarg!=0 && numarg!=1 )
    gpterror( "Syntax: %s([seed])\n", gptgetname(init) ) ;

  /* Read seed or calulate from system time */
  if( numarg==0 )
    seed = (unsigned)time(NULL) ;
  else
    seed = gptgetargint(init,1) ;

  /* Initialise randomizer */
  intpprandsetseed(seed) ;
}
