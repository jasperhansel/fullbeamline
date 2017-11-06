/* spacecharge3Dtree.c - Barnes&Hut relativistic treecode */

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "elem.h"

#include "../../kernel/BHtree.h"

#define TREEMINPAR 64
struct sctree_info
{
  /* BH parameters */
  double theta ;

  /* Zmin, zmax not implemented (yet) */
  double zmin ;
  double zmax ;

  /* Communication with restframeparams() */
  double Gmax ;
  double GBo[3] ;
  double Go ;
} ;

/* Forward declarations */
static int restframeparams(struct sctree_info *info) ;
static int testrestframe_out(double t, double *dt, double *x, struct sctree_info *info) ;
static void sctree_sim( double t, double *x, double *p, void *info ) ;

void spacecharge3Dtree_init(gptinit *init)
{
  struct sctree_info *info ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(theta)\n", gptgetname(init) ) ;

  info = (struct sctree_info *)gptmalloc( sizeof(struct sctree_info) ) ;

  info->theta   = gptgetargdouble(init,1) ;

  info->zmin =-DBL_MAX ;
  info->zmax = DBL_MAX ;
  info->Gmax = 2 ; // For now, we just set Gmax=2

  odeaddfprfunction( ODEFNC_USR, sctree_sim, info ) ;
  odeaddoutfunction( ODEFNC_USR, (odeoutfnc)testrestframe_out, info ) ;
}

/* Main routine to calculate B&H interactions */
static void sctree_sim( double t, double *x, double *p, void *vinfo )
{
  struct sctree_info *info = (struct sctree_info *)vinfo ;

// Get rest-frame parameters
  if( restframeparams( info )) return ;

// Get number of paticles in simulation
  numalivepar = 0 ;
  for(int i=0 ; i<numpar ; i++) if( pars[i].alive ) numalivepar++ ;

// Defer to point-to-point method if N<10
  if( numalivepar<TREEMINPAR )
  {
    sc3Dp2p(t,x,p,info->zmin,info->zmax, 0,1) ;
    return ;
  }

// Get bunch center of particles (not corrected for n), and rmacro
  double avgr[3] = {0,0,0} ;
  double eps2 = pars[0].r2 ; // OK, len!=0 as just checked above

  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    for(int j=0 ; j<3 ; j++) 
	  avgr[j] += pars[i].Wr[j] ;

    if( pars[i].r2 != eps2 )
		gpterror( "Particle radii must be identical. Use 'setrmacrodist(\"%s\",\"u\",rmacro,0)." ) ;
  }

  for(int j=0 ; j<3 ; j++)
    avgr[j] /= numalivepar ;

  if( eps2==0 )
  {
	  gptwarning( "Particles have zero radius resulting in singularities in the field.\n" ) ;
	  gpterror( "Use 'setrmacrodist(\"set\",\"u\",rmacro,0)' with positive rmacro.\n" ) ;
  }

// Create particle list (single core)
  pointchargearray charges(eps2) ;
  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    // Get rest frame position
	double ro[3], rc[3] ;
    for(int j=0 ; j<3 ; j++) rc[j] = pars[i].Wr[j] - avgr[j] ;
    double alpha = gptVECINP(rc,info->GBo)/(info->Go + 1) ;
    for(int j=0 ; j<3 ; j++ ) ro[j] = rc[j] + alpha*info->GBo[j] ;

	// Insert charge
    pointcharge q(ro,pars[i].n*pars[i].q) ;
    charges.insert(q) ;
  }

// Create BHtree
  auto_ptr<electrostaticfield> BHtree = BHnode::createBHelectrostaticfield(charges,info->theta) ;

// Store results and transform back to lab frame
#pragma omp parallel for
  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    // Get rest frame position (again, not so clever)
	double ro[3], rc[3] ;
    for(int j=0 ; j<3 ; j++) rc[j] = pars[i].Wr[j] - avgr[j] ;
    double alpha = gptVECINP(rc,info->GBo)/(info->Go + 1) ;
    for(int j=0 ; j<3 ; j++ ) ro[j] = rc[j] + alpha*info->GBo[j] ;

	// Get field from BH
	vec3 Wr(ro) ;
	vec3 E = BHtree->getEfield(Wr) ;

    double cdot = gptVECINP(E,info->GBo) / (info->Go+1.0) ;
 
    pars[i].WE[0] += info->Go*E[0] - cdot*info->GBo[0] ;
    pars[i].WE[1] += info->Go*E[1] - cdot*info->GBo[1] ;
    pars[i].WE[2] += info->Go*E[2] - cdot*info->GBo[2] ;

    pars[i].WB[0] += (info->GBo[1]*E[2] - info->GBo[2]*E[1]) / gpt_c ;
    pars[i].WB[1] += (info->GBo[2]*E[0] - info->GBo[0]*E[2]) / gpt_c ;
    pars[i].WB[2] += (info->GBo[0]*E[1] - info->GBo[1]*E[0]) / gpt_c ;
  }
}


/* Calculate rest frame parameters */
static int restframeparams(struct sctree_info *info)
{
  double snGB[3], Bo[3] ;
  static int sn0warning=0 ;

  /* Calculate sumG and sumGB, weighted by nmacro */
  double sn  = 0 ;
  double snG = 0 ;
  for(int j=0 ; j<3 ; j++) snGB[j] = 0 ;
  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    sn  += pars[i].n ;
    snG += pars[i].n * pars[i].G ;
    for(int j=0 ; j<3 ; j++) snGB[j] += pars[i].n * pars[i].GBr[j] ;
  }

  // Return when no charge is present
  if(sn==0) 
  {
    if(sn0warning==0)
	{
      gptwarning( "No charge present, space-charge is ignored\n") ;
      sn0warning = 1 ;
	}
    return(1) ;
  }

  // Resume when charge 'appears' and disable future warnings
  if( sn0warning==1 ) 
  {
    gptwarning( "Resuming space-charge calculations\n" ) ;
    sn0warning = 2 ;
  }

  /* Calculate rest frame parameters */
  for(int j=0 ; j<3 ; j++) Bo[j] = snGB[j]/snG ; 
  for(int j=0 ; j<3 ; j++) info->GBo[j] = Bo[j]/sqrt(1-gptVECSQR(Bo)) ;
  info->Go = sqrt(1 + gptVECSQR(info->GBo) ) ;

  return(0) ;
}


/* Test rest frame velocities. Called using OUT interface at beginning of timestep */
static int testrestframe_out(double t, double *dt, double *x, struct sctree_info *info)
{
  int nnewmaxspeed ;
  static int noldmaxspeed=0 ;

  double Gio ;

  /* Calculate number of speed offenders in rest frame */
  nnewmaxspeed=0 ;
  for(int i=0 ; i<numpar; i++ ) if( pars[i].alive )
  {
    Gio = pars[i].G*info->Go-gptVECINP(pars[i].GBr,info->GBo) ;

    if( Gio>info->Gmax ) nnewmaxspeed++ ; 
  }

  /* Print warning in verbose mode or when the number of speed offenders varies by 20% */
  if( (verbose && nnewmaxspeed) || nnewmaxspeed<noldmaxspeed*10/12 || nnewmaxspeed>noldmaxspeed*12/10 )
  {
    gptwarning("Spacecharge3Dtree: %d particles with Lorentz factor in rest frame over %f at t=%e\n", nnewmaxspeed, info->Gmax, t) ;
    noldmaxspeed = nnewmaxspeed ;
  }

  return(0) ;
}
