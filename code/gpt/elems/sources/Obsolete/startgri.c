/* startgrid.c - Start with a rectangular grid */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "elem.h"


void startgrid_init(gptinit *init)
{
  double Lx,dx,Ly,dy,Lz,dz ;
  int i,j,k,im,jm,km ;
  double r[3], Br[3] ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=6 && numarg!=9 )
    gpterror( "Syntax: %s(ECS,Lx,dx,Ly,dy,Lz,dz,[Bx,By,Bz])\n", gptgetname(init) ) ;

  /* Read parameters */
  Lx=gptgetargdouble(init,1) ;
  dx=gptgetargdouble(init,2) ;
  Ly=gptgetargdouble(init,3) ;
  dy=gptgetargdouble(init,4) ;
  Lz=gptgetargdouble(init,5) ;
  dz=gptgetargdouble(init,6) ;

  if( numarg==9 )
  {
    Br[0]=gptgetargdouble(init,7) ;
    Br[1]=gptgetargdouble(init,8) ;
    Br[2]=gptgetargdouble(init,9) ;
  } else Br[0]=Br[1]=Br[2]=0.0 ;

  /* Guard the algorith from dividing zero by zero */
  if( Lx==0.0 ) dx=1.0 ;
  if( Ly==0.0 ) dy=1.0 ;
  if( Lz==0.0 ) dz=1.0 ;

  /* Calculate boundaries */
  im=(int)((1+10*DBL_EPSILON)*Lx/dx/2) ;
  jm=(int)((1+10*DBL_EPSILON)*Ly/dy/2) ;
  km=(int)((1+10*DBL_EPSILON)*Lz/dz/2) ;

  /* Check proper use of parameters */
  if( im<0 || jm<0 || km<0 )
    gpterror( "Invalid arguments\n" ) ;

  /* Position particles */
  for(k=-km ; k<=km ; k++)
  {
    r[2]=k*dz ;
    for(j=-jm ; j<=jm ; j++)
    {
      r[1]=j*dy ;
      for(i=-im ; i<=im ; i++)
      {
        r[0]=i*dx ;
        gptaddparticle( init, r ,Br ) ;
      }
    }
  }
}
