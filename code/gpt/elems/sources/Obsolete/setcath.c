/* setcath.c - Modify a particle set according to a spherical cathode */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void setcathode_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len, nerror ;
  double Ra, Rc, kT ;
  double rxy, d ;

  /* Get parameters */
  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(set,Ra,Rc,kT)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  Ra   = gptgetargdouble(init,2) ;
  Rc   = gptgetargdouble(init,3) ;
  kT   = gptgetargdouble(init,4) ;

  /* Check parameters */
  if( Ra<=0.0 )
    gpterror( "%s: Ra must be positive\n", gptgetname(init) ) ;
  if( Ra>=fabs(Rc) && Rc!=0.0 )
    gpterror( "%s: Ra must be smaller than Rc\n", gptgetname(init) ) ;
  if( kT!=0.0 )
    gpterror( "%s: The kT parameter not supported yet\n", gptgetname(init) ) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gpterror( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Make spherical cathode in real space */
  nerror = 0 ;
  for( i=0 ; i<len ; i++ )
  {
    rxy = sqrt( par[i].Wr[0]*par[i].Wr[0]+par[i].Wr[1]*par[i].Wr[1] ) ;

    /* Test radius within Ra */
    if( rxy >= Ra )
    {
      nerror++ ;
      par = gptremovepar(set,i--,&len) ;
      continue ;
    }

    /* Modify coordinates */
    if( Rc!=0.0 )
    {
      d = sqrt(Rc*Rc-rxy*rxy) - sqrt(Rc*Rc-Ra*Ra) ; /* Always positive */
      par[i].Wr[2] += (Rc>0 ? d : -d) ;
    }
  }

  if( nerror != 0 )
    gptwarning( "%s: %d particles outside Ra are removed\n", gptgetname(init), nerror ) ;
}
