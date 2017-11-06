/* stdxyzmax.c: Remove deviated particles */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct stdxyzmax_info
{
  double Nstdx ;
  double Nstdy ;
  double Nstdz ;
} ;

static void stdxyzmax_end(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo) ;
static double stdsqrt(double x) ;

void stdxyzmax_init(gptinit *init)
{
  struct stdxyzmax_info *info ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(Nstdx,Nstdy,Nstdz)\n", gptgetname(init) ) ;

  info = (struct stdxyzmax_info *)gptmalloc( sizeof(struct stdxyzmax_info) ) ;

  info->Nstdx = gptgetargdouble(init,1) ;
  info->Nstdy = gptgetargdouble(init,2) ;
  info->Nstdz = gptgetargdouble(init,3) ;

  if( info->Nstdx<0 || info->Nstdy<0 || info->Nstdz<0 )
    gpterror( "Nstdx, Nstdy and Nstdz must be equal to or greater than zero.\n") ;

  /* Register the routine removing the appropriate particles after each succesful step */
  odeaddendfunction( ODEFNC_USR-1, (odeendfnc)stdxyzmax_end, info ) ;
}


/* The following routine removes all particles outside x,y,z standard deviations */
static void stdxyzmax_end(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo)
{
  struct stdxyzmax_info *info = (struct stdxyzmax_info *)vinfo ;

  double Nstdx, Nstdy, Nstdz ;
  double sumx, sumx2, avgx, stdx ;
  double sumy, sumy2, avgy, stdy ;
  double sumz, sumz2, avgz, stdz ;

  if(numpar<3) return ;

  Nstdx = info->Nstdx ;
  Nstdy = info->Nstdy ;
  Nstdz = info->Nstdz ;

  sumx = sumx2 = 0 ;
  sumy = sumy2 = 0 ;
  sumz = sumz2 = 0 ;
  unsigned int N = 0 ;

  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    sumx += pars[i].Wr[0] ; sumx2 += pars[i].Wr[0]*pars[i].Wr[0] ;
    sumy += pars[i].Wr[1] ; sumy2 += pars[i].Wr[1]*pars[i].Wr[1] ;
    sumz += pars[i].Wr[2] ; sumz2 += pars[i].Wr[2]*pars[i].Wr[2] ;
    N++ ;
  }
  if( N==0 ) return ;

  avgx = sumx/N ;  
  avgy = sumy/N ;  
  avgz = sumz/N ;  

  stdx = stdsqrt( (sumx2 - sumx*avgx)/N ) ;
  stdy = stdsqrt( (sumy2 - sumy*avgy)/N ) ;
  stdz = stdsqrt( (sumz2 - sumz*avgz)/N ) ;

  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    if(Nstdx>0 && pars[i].Wr[0] < avgx - Nstdx*stdx) gptremoveparticle(&pars[i]) ;
    if(Nstdx>0 && pars[i].Wr[0] > avgx + Nstdx*stdx) gptremoveparticle(&pars[i]) ;
    if(Nstdy>0 && pars[i].Wr[1] < avgy - Nstdy*stdy) gptremoveparticle(&pars[i]) ;
    if(Nstdy>0 && pars[i].Wr[1] > avgy + Nstdy*stdy) gptremoveparticle(&pars[i]) ;
    if(Nstdz>0 && pars[i].Wr[2] < avgz - Nstdz*stdz) gptremoveparticle(&pars[i]) ;
    if(Nstdz>0 && pars[i].Wr[2] > avgz + Nstdz*stdz) gptremoveparticle(&pars[i]) ;
  }
}

static double stdsqrt(double x)
{
  if(x<0) return(0) ;
  return( sqrt(x) ) ;
}
