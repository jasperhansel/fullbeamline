/* Startcat.c - Start particles emerging from a cathode */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct startcat_info
{
  gpttransform tf ;
  gptparset *ppp ;
  double Ra, Rc ;
  double La, Lc ;
  double dzdt ;
  double dtmax ;
  int len ;
  gptinitpar *pars ;
} ;

static int parssortz(const void *a, const void *b) ;
static void startcat_endusr(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo) ;
static int startcat_outusr(double t, double *dt, double *x, void *vinfo) ;

void startcathode_init(gptinit *init)
{
  struct startcat_info *info ;
  char *name ;
  double Ra, Rc, La, Lc, dzdt, dtmax ;
  double rxy, d, e ;
  gptparset *ppp ;
  int i, len, nerror ;
  gptinitpar *pars ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=6 )
    gpterror( "Syntax: %s(ECS,set,Ra,Rc,Lc,dzdt,dtmax)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  Ra   = gptgetargdouble(init,2) ;
  Rc   = gptgetargdouble(init,3) ;
  Lc   = gptgetargdouble(init,4) ;
  dzdt = gptgetargdouble(init,5) ;
  dtmax= gptgetargdouble(init,6) ;


  /* Check parameters */
  if( gpttestparset( name )==NULL )
    gpterror( "The particle set \"%s\" does not exist\n", name ) ;
  if( Ra<=0.0 )
    gpterror( "%s: Ra must be positive\n", gptgetname(init) ) ;
  if( Ra>=fabs(Rc) && Rc!=0.0 )
    gpterror( "%s: Ra must be smaller than Rc\n", gptgetname(init) ) ;
  if( dzdt==0.0 )
    gpterror( "%s: dzdt must be nonzero\n", gptgetname(init) ) ;
  if( dtmax<0.0 )
    gpterror( "%s: dtmax must be zero or positive\n", gptgetname(init) ) ;

  /* Calculate laser aperture */
  La = Ra ;
  if( Lc!=0.0 && La>fabs(Lc) ) La = fabs(Lc) ;

  /* Allocate and fill info structure */
  info = (struct startcat_info *)gptmalloc( sizeof(struct startcat_info) ) ;
  info->Ra   = Ra ;
  info->Rc   = Rc ;
  info->La   = La ;
  info->Lc   = Lc ;
  info->dzdt = dzdt ;
  info->dtmax= dtmax ;
  gptconcattransform( &info->tf, &init->paxis->a, &init->e ) ;

  /* Undo spherical cathode in real space, compensate for laser arrival time */
  ppp = gptgetparset( name ) ;
  pars = gptgetparsetpars( ppp,&len ) ;
  nerror = 0 ;
  for( i=0 ; i<len ; i++ )
  {
    rxy = sqrt( pars[i].Wr[0]*pars[i].Wr[0]+pars[i].Wr[1]*pars[i].Wr[1] ) ;

    /* Test radius within Ra and La */
    if( rxy>=Ra || rxy>=La )
    {
      nerror++ ;
      pars = gptremovepar(ppp,i--,&len) ;
      continue ;
    }

    /* Modify coordinates */
    d=0.0 ;
    if( Rc!=0.0 )
    {
      d = sqrt(Rc*Rc-rxy*rxy) - sqrt(Rc*Rc-Ra*Ra) ; /* Always positive */
      d *= (Rc>0 ? 1 : -1) ; /* Same sign as "setcathode", but subtracted */
      pars[i].Wr[2] -= d ;
    }

    /* Get laser curvature */
    e=0.0 ;
    if( Lc!=0.0 )
    {
      e = sqrt(Lc*Lc-rxy*rxy) - sqrt(Lc*Lc-La*La) ; /* Always positive */
      e *= (Lc>0 ? 1 : -1) ;
    }
    /* Compesate for laser arrival time */
    pars[i].Wr[2] -= (e-d)*dzdt/gpt_c ;
  }
  if( nerror != 0 )
    gptwarning( "%s: %d particles outside Ra and La are removed\n", gptgetname(init), nerror ) ;

  /* Move all particles in set to buffer */
  pars = gptgetparsetpars(ppp,&len) ;
  info->len = len ;
  info->pars = (struct gptinitpar *)gptmalloc(len*sizeof(*info->pars)) ;
  for( i=0 ; i<len ; i++ ) info->pars[i] = pars[i] ;
  gptremoveparset(ppp) ;
  info->ppp = gptcreateparset( gptgetargstring(init,1) ) ;

  /* Convert z[m]->t[s], sort all particles, shift max to z=0 */
  for(i=0 ; i<len ; i++) info->pars[i].Wr[2] /= dzdt ;
  qsort( info->pars,len,sizeof(*info->pars),parssortz) ;
  for(i=0 ; i<len ; i++) info->pars[i].Wr[2] -= info->pars[len-1].Wr[2] ;

  /* Start all partcles with z=0 and register END function */
  startcat_endusr(0.0,0.0,NULL,NULL,NULL,info,NULL) ;
  odeaddendfunction( ODEFNC_USR, startcat_endusr, info ) ;
  odeaddoutfunction( ODEFNC_USR, startcat_outusr, info ) ;
}

/* Sort particles by z-coordinate in ascending order */
static int parssortz(const void *a, const void *b)
{
  gptinitpar *para = (gptinitpar *)a ;
  gptinitpar *parb = (gptinitpar *)b ;

  if( para->Wr[2] > parb->Wr[2] ) return( 1) ;
  if( para->Wr[2] < parb->Wr[2] ) return(-1) ;
  return( 0 ) ;
}

/* Add due particles to the simulation, Unit of Wr[2] starts as [s]! */
static void startcat_endusr(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo)
{
  struct startcat_info *info = (struct startcat_info *)vinfo ;

  gptinitpar *par ;
  double WCSWr[3], WCSWGBr[3] ;
  double G, deltat ;
  double rxy, d, e ;
  int i ;

  if( info->len==0 ) return ;

  do
  {
    par = &info->pars[info->len-1] ;
    if( par->Wr[2]+tend>=0 )
    {
      /* Extrapolate position to t=t, convert back to [m] */
      G = sqrt(1.0+gptVECSQR(par->GBr)) ; /* Direction is irrelevant for G */
      deltat = tend+par->Wr[2] ;
      par->Wr[2] = 0.0 ;
      for(i=0 ; i<3 ; i++) par->Wr[i] += deltat*gpt_c*par->GBr[i]/G ;

      /* Particle distance from axis */
      rxy = sqrt( par->Wr[0]*par->Wr[0]+par->Wr[1]*par->Wr[1] ) ;

      /* Redo shperical cathode in real space */
      d=0.0 ;
      if( info->Rc!=0.0 )
      {
        d = sqrt(info->Rc*info->Rc-rxy*rxy) - sqrt(info->Rc*info->Rc-info->Ra*info->Ra) ; /* Always positive */
        d *= (info->Rc>0 ? 1 : -1) ; /* Same sign as "setcathode" */
        par->Wr[2] += d ;
      }

      /* Get laser curvature */
      e=0.0 ;
      if( info->Lc!=0.0 )
      {
        e = sqrt(info->Lc*info->Lc-rxy*rxy) - sqrt(info->Lc*info->Lc-info->La*info->La) ; /* Always positive */
        e *= (info->Lc>0 ? 1 : -1) ;
      }
      /* Compesate for laser arrival time */
      par->Wr[2] += (e-d)*info->dzdt/gpt_c ;

      /* Transform particle coordinates to WCS */
      gpttoWCS(&info->tf,par->Wr,WCSWr) ;
      gptdirectiontoWCS(&info->tf,par->GBr,WCSWGBr) ;

      /* Add particle to simulation and update buffer statistics */
      gptaddparmqnar( info->ppp, WCSWr, WCSWGBr, par->m, par->q, par->n, par->paxis, par->r ) ;
      info->len-- ;
    } else break ;
  } while(info->len!=0) ;
}

static int startcat_outusr(double t, double *dt, double *x, void *vinfo)
{
  struct startcat_info *info = (struct startcat_info *)vinfo ;
  gptinitpar *par ;

/* Adapt GPT timestep to next particle release */
  if( info->len==0 ) return( 0 ) ;
  par = &info->pars[info->len-1] ;
  if( *dt > (info->dtmax-par->Wr[2]-t) )
    *dt = (info->dtmax-par->Wr[2]-t)*(1+DBL_EPSILON) ;

  return(0) ;
}
