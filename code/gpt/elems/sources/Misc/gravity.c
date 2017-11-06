/* gravity.c: Calculate gravitational forces */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct gravity_info
{
  double gx ;
  double gy ;
  double gz ;
} ;

static void gravity_sim( double t, double *x, double *p, void *info ) ;

void gravity_init(gptinit *init)
{
  struct gravity_info *info ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(gx,gy,gz)\n", gptgetname(init) ) ;

  info = (struct gravity_info *)gptmalloc( sizeof(struct gravity_info) ) ;

  info->gx = gptgetargdouble(init,1) ;
  info->gy = gptgetargdouble(init,2) ;
  info->gz = gptgetargdouble(init,3) ;

  /* Register gravity field function */
  odeaddfprfunction( ODEFNC_INT, gravity_sim,info ) ;
}


static void gravity_sim( double t, double *x, double *p, void *vinfo )
{
  struct par *par ;
  double Br[3], beta2, BinpG ;
  double gx, gy, gz ;
  struct gravity_info *info = (struct gravity_info *)vinfo ;

  gx = info->gx ;
  gy = info->gy ;
  gz = info->gz ;

  for(int i=0 ; i<numpar ; i++)
  {
    par = &pars[i] ;

    Br[0] = par->GBr[0]/par->G ;
    Br[1] = par->GBr[1]/par->G ;
    Br[2] = par->GBr[2]/par->G ;

    beta2 = Br[0]*Br[0] + Br[1]*Br[1] + Br[2]*Br[2] ;
    BinpG = Br[0]*gx + Br[1]*gy + Br[2]*gz ;

    par->WF[0] += par->G*par->m * ((1+beta2)*gx - BinpG*Br[0] ) ;
    par->WF[1] += par->G*par->m * ((1+beta2)*gy - BinpG*Br[1] ) ;
    par->WF[2] += par->G*par->m * ((1+beta2)*gz - BinpG*Br[2] ) ;
 
  }
}
