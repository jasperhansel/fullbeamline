/* setcurvature.c - Modify a particle set according to a spherical cathode  and/or laser front */

#include <stdio.h>
#include <math.h>
#include "elem.h"

void setcurvature_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len, nerror ;
  double Ra, Rc, Lc, La ;
  double rxy, d, e ;

  /* Get parameters */
  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(set,Ra,Rc,Lc)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  Ra   = gptgetargdouble(init,2) ;
  Rc   = gptgetargdouble(init,3) ;
  Lc   = gptgetargdouble(init,4) ;

  /* Check parameters */
  if( Ra<=0.0 )
    gpterror( "%s: Ra must be positive\n", gptgetname(init) ) ;
  if( Ra>=fabs(Rc) && Rc!=0.0 )
    gpterror( "%s: Ra must be smaller than Rc\n", gptgetname(init) ) ;

  /* Calculate laser aperture */
  La = Ra ;
  if( Lc!=0.0 && La>fabs(Lc) ) La = fabs(Lc) ;

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
    if( rxy >= Ra || rxy>=La )
    {
      nerror++ ;
      par = gptremovepar(set,i--,&len) ;
      continue ;
    }

    /* Cathode curvature */
    d=0.0 ;
    if( Rc!=0.0 )
    {
      d = sqrt(Rc*Rc-rxy*rxy) - sqrt(Rc*Rc-Ra*Ra) ; /* Always positive */
      if( Rc<0 ) d=-d ;
    }

    /* Laser curvature */
    e=0.0 ;
    if( Lc!=0.0 )
    {
      e = sqrt(Lc*Lc-rxy*rxy) - sqrt(Lc*Lc-La*La) ; /* Always positive */
      if( Lc<0 ) e=-e ;
    }

    /* Modify coordinates */
    par[i].Wr[2] += d ;
    par[i].tstart += (e-d)/gpt_c ;
  }

  if( nerror != 0 )
    gptwarning( "%s: %d particles outside Ra or Lc are removed\n", gptgetname(init), nerror ) ;
}
