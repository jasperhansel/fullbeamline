/* setellipse.c - Set homogeneous ellipse */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

extern double dblpulsarrand(void) ;

void setellipse_init(gptinit *init)
{
  double a,b,c ;
  double x,y,z ;
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(set,a,b,c)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  a    = gptgetargdouble(init,2) ;
  b    = gptgetargdouble(init,3) ;
  c    = gptgetargdouble(init,4) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Set ellipse */
  for( i=0 ; i<len ; i++ )
  {
    do
    {
      /* Uniform in box between -a and a, -b and b, -c and c */
      x = a*(2*dblpprand()-1) ;
      y = b*(2*dblpprand()-1) ;
      z = c*(2*dblpprand()-1) ;
    } 
    while( x*x/(a*a)+y*y/(b*b)+z*z/(c*c) >= 1 ) ;

    par[i].Wr[0] = x ;
    par[i].Wr[1] = y ;
    par[i].Wr[2] = z ;
  }
}
