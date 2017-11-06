/* setmirror.c: Mirror all particles in a paticle set. */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Initialization routine */
void setmirror_init(gptinit *init)
{
  char *name ;
  double s[3], n[3], nabs ;
  int copy ;

  gptparset *set ;
  gptinitpar par, *pars ;
  int i, j, len, newlen ;

  double Wrmir[3], GBrmir[3], inp ;

  if( gptgetargnum(init)!=8 )
    gpterror( "Syntax: %s(set,sx,sy,sz,nx,ny,nz,copy)\n", gptgetname(init) ) ;

  /* Read all parameters as doubles and store them in info structure */
  name = gptgetargstring(init,1) ;
  s[0] = gptgetargdouble(init,2) ;
  s[1] = gptgetargdouble(init,3) ;
  s[2] = gptgetargdouble(init,4) ;
  n[0] = gptgetargdouble(init,5) ;
  n[1] = gptgetargdouble(init,6) ;
  n[2] = gptgetargdouble(init,7) ;
  copy = gptgetargint   (init,8) ;

  /* Normalize n */
  nabs = gptVECLEN(n) ;
  if( nabs==0.0 )
    gpterror( "%s: Normal (nx,ny,nz) has zero length", gptgetname(init) ) ;
  for(j=0 ; j<3 ; j++) n[j] /= nabs ;  

  /* Get particle set particles */
  if( gpttestparset( name )==NULL )
    gptwarning( "%s: The particle set \"%s\" does not exist\n", gptgetname(init), name ) ;
  set = gptgetparset( name ) ;
  pars = gptgetparsetpars( set, &len ) ;

  /* Mirror */
  for(i=0 ; i<len ; i++)
  {
    pars = gptgetparsetpars( set, &newlen ) ;
    par = pars[i] ;

    for(j=0 ; j<3 ; j++) Wrmir[j]  = par.Wr[j]-s[j] ;
                               inp = gptVECINP(Wrmir,n) ;
    for(j=0 ; j<3 ; j++) Wrmir[j] -= 2*inp*n[j] ; ;
    for(j=0 ; j<3 ; j++) Wrmir[j] += s[j] ;

                               inp = gptVECINP(par.GBr,n) ;
    for(j=0 ; j<3 ; j++) GBrmir[j] = par.GBr[j]-2*inp*n[j] ;

	gptaddparmqnartid(set,Wrmir,GBrmir,par.m,par.q,par.n,par.paxis,par.r,par.tstart,par.ID) ;
  }

  /* Optionally remove */
  if( copy ) return ;

  for(i=0 ; i<len ; i++)
    pars = gptremovepar(set,i,&newlen) ;
}
