/* dtmaxpars.c: Start a maximum of N particles in each step */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct dtmaxpars_info
{
  int Nmax ;
} ;

static int dtmaxpars_out(double t, double *dt, double *x, void *vinfo) ;

void dtmaxpars_init(gptinit *init)
{
  struct dtmaxpars_info *info ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(Nmax)\n", gptgetname(init) ) ;

  info = (struct dtmaxpars_info *)gptmalloc( sizeof(struct dtmaxpars_info) ) ;

  info->Nmax = gptgetargint(init,1) ;

  if( info->Nmax!=1 )
    gpterror( "Only Nmax=1 is implemented so far. Sorry." ) ;

  odeaddoutfunction( ODEFNC_USR, dtmaxpars_out, info ) ;
  gptaddmainfunction( GPTMAINFNC_TER, gptfree, info ) ;
}

static int dtmaxpars_out(double t, double *dt, double *x, void *vinfo)
{
//  struct dtmaxpars_info *info = (struct dtmaxpars_info *)vinfo ; // Not used to far

  double dt1 = *dt ;
  double dt2 = *dt ;

  for(int i=0 ; i<numpar ; i++)
  {
    if( pars[i].persistent && !pars[i].alive ) // Particle not added to simulation yet
    {
        if( t + *dt + *dt > pars[i].tstart )
        {
          double newdt = 0.5*(pars[i].tstart-t) ;
          if( dt2>newdt ) dt2=newdt ;
        }

      if( t + *dt > (1.0-DBL_EPSILON)*pars[i].tstart )
      {
        double newdt = (1.0+DBL_EPSILON)*(pars[i].tstart-t) ;
        if( dt1>newdt ) 
        {
//        gptwarning("particle start-time %le", pars[i].tstart) ;
          dt1=newdt ;
        }
      }
    }
  }

  if( *dt>dt1 )
  {
//  gptwarning("Timestep reduced from %le to %le by nearest particle", *dt, dt1 ) ;
    *dt=dt1 ;
  } else if( *dt>dt2 )
  {
//  gptwarning("Timestep optimized from %le to %le by nearest particle", *dt, dt1 ) ;
    *dt=dt2 ;
  }
  
  return( 0 ) ;
}
