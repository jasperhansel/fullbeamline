/* setscale.c: Scale initial particle coordinates */

#include <stdio.h>
#include <math.h>
#include "elem.h"

void setscale_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int numarg, i, len ;
  double scalex, scaley, scalez ;
  double scaleGBx, scaleGBy, scaleGBz ;
  double scalet ;
  double scalen, scaler ;

  /* Check syntax */
  numarg = gptgetargnum(init) ;
  if( numarg!=7 && numarg!=8 && numarg!=10 )
    gpterror( "Syntax: %s(set,x,y,y,GBx,GBy,GBz,[t,[nmacro,rmacro]])\n", gptgetname(init) ) ;

  /* Get particle set and scala parameters */
  name     = gptgetargstring(init, 1) ;
  scalex   = gptgetargdouble(init, 2) ;
  scaley   = gptgetargdouble(init, 3) ;
  scalez   = gptgetargdouble(init, 4) ;
  scaleGBx = gptgetargdouble(init, 5) ;
  scaleGBy = gptgetargdouble(init, 6) ;
  scaleGBz = gptgetargdouble(init, 7) ;
  scalet   = 1.0 ;
  scalen   = 1.0 ;
  scaler   = 1.0 ;

  if( numarg==8 || numarg==10 )
    scalet = gptgetargdouble(init, 8) ;

  if( numarg==10 )
  {
    scalen = gptgetargdouble(init, 9) ;
    scaler = gptgetargdouble(init,10) ;
  }

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Scale initial distribution */
  for( i=0 ; i<len ; i++ )
  {
    par[i].Wr[0]  *= scalex ;
    par[i].Wr[1]  *= scaley ;
    par[i].Wr[2]  *= scalez ;
    par[i].GBr[0] *= scaleGBx ;
    par[i].GBr[1] *= scaleGBy ;
    par[i].GBr[2] *= scaleGBz ;

    par[i].tstart *= scalet ;

    par[i].n      *= scalen ;
    par[i].r      *= scaler ;
  }
}
