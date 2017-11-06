/* settcopy.c: Copy a particle set in TIME domain */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Initialization routine */
void settcopy_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar par, *pars ;
  char *name ;
  int N, i, j, len, tmp ;
  double dt ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(set,N,dist)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  N    = gptgetargint(init,2) ;
  dt   = gptgetargdouble(init,3) ;

  if( N<1 ) terminate( "%s: N must be 1 or more\n", gptgetname(init)  ) ;

  /* Get particle set particles */
  if( gpttestparset( name )==NULL )
    gptwarning( "%s: The particle set \"%s\" does not exist\n", gptgetname(init), name ) ;
  set = gptgetparset( name ) ;
  pars = gptgetparsetpars( set, &len ) ;

  /* Multi-copy particle set */ 
  for(j=0;j<len;j++)
  {
    par = pars[j] ;

    for(i=1;i<N;i++)
    {
      par.tstart = pars[j].tstart + i*dt ;
	  gptaddparmqnartid(set,par.Wr,par.GBr,par.m,par.q,par.n,par.paxis,par.r,par.tstart,par.ID) ;
      pars = gptgetparsetpars( set, &tmp ) ;
    }
  }
}
