/* setcopy.c: Copy a particle set */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Initialization routine */
void setcopy_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar par, *pars ;
  char *name ;
  int N, i, j, len, tmp ;
  double dz ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(set,N,dist)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  N    = gptgetargint(init,2) ;
  dz   = gptgetargdouble(init,3) ;

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
      par.Wr[2] = pars[j].Wr[2] + i*dz ;
      gptaddparmqn(set,par.Wr,par.GBr,par.m,par.q,par.n) ;
      pars = gptgetparsetpars( set, &tmp ) ;
    }
  }
}
