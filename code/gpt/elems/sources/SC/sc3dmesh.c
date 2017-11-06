/* sc3dmesh.c: Driver routine for the Rostock 3D multigrid space-charge solver
 *
 * Authors: Bas van der Geer, Marieke de Loos ; Pulsar Physics
 * 
 * Changes:
 * 22-Apr-02: SCDEBUG definitions
 * 22-Apr-02: stdx,y,z -> stdmax
 * 21-Jun-02: Equidistant grid when fn=0
 * 21-Jun-02: Return immediately with 1 particle
 * ??-???-??: Fixed interpolation error near boundary
 * 16-Sep-02: Optional cathode fieldline at z=0
 * 17-Sep-02: Boundary conditions to poisson solver
 * ??-???-??: Beam scaling
 * 01-Oct-02: logBox/Cube boundary conditions, Verror parameter
 * 10-Jan-03: Aditional mesh smoothing, removed because it didn't help
 * 01-feb-03: Variable number of meshlines between bunch and bounding box
 * 03-feb-03: All paramters optional
 * 12-feb-03: New interface with Poisson solver
 * 03-Mar-03: Previous solution to Poisson solver
 * 12-Mar-03: Pass solvertype
 * 07-May-03: Maximum aspect ratio
 * 23-May-03: Block memory management for Poisson solver
 * 27-May-03: Getindex buffering
 * 08-Aug-03: Much faster 'testrestframe'
 * 28-Aug-03: DirichletZero
 * 16-Oct-03: Bounding box rewrite
 * 16-Oct-03: Dirichlet boundary conditions for 'brick' and 'sphere' cases
 * 16-Apr-04: TEST: Wavelet filter for charge density
 * 26-May-04: Optionally read meshline positions from file
 * 29-Jul-04: Specify boundary conditions individually
 * 15-Aug-04: Cathode boundary condition correct
 * 24-aug-04: Rest frame error only for succesful timesteps (via OUT interface)
 * 04-Jan-08: Zero charge warning corrected
 * 28-Jan-08: Defer to p2p for few particles
 * 29-jun-09: Rotate mesh into direction of propagation
 * 05-nov-09: Meshing robustified by skipping a fraction 1/2N
 * 05-nov-09: Mixed charge capabilities
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include "elem.h"
// #include "poisson.h"

//#define SCDEBUG_MESHLINES				/* Write mesh iteration information */
//#define SCDEBUG_RESTSTATS				/* Write rest-frame information */
//#define SCDEBUG_MESHSTATS				/* Write histogram of particles/meshnode */
//#define SCDEBUG_MEMORY				/* Debug memory management */

//#define SCDEBUG_MESHFILE "mesh3d.scr"			/* ACAD script to plot mesh and particles */
//#define SCDEBUG_CHARFILE "charge3d.scr"		/* ACAD script to plot charge-density */
//#define SCDEBUG_POTEFILE "potential3d.scr"		/* ACAD script to plot potential */

//#define SCDEBUG_READMESHFROMFILE "grid_new.dat"		/* Read meshline positions from file */

#define SC3DMESH_POTENTIALBRICK				/* 'brick' expression to estimate potential with Dirichlet boundary conditions */
//#define SC3DMESH_POTENTIALSPHERE			/* 'sphere' expression */

//#define WAVELET_INFILE  "charge.fits"
//#define WAVELET_OUTFILE "filter.fits"
//#define WAVELET_COMMAND "\\progra~1\\wavelets\\mr3d_pfilter -m 1 charge.fits filter.fits"

#define MESHPREC 1e-6					/* Precision of meshlines position */
#define MAXMESHIT 1000					/* Maximum number of mesh iterations */


struct parmeshinfo
{
  unsigned int jprev[3] ;
  int alignment ;
} ;

struct sc3dmesh_info
{
  /* Commandline parameters */
  unsigned int Nc[3] ;
  double nstd ;
  double macc ;
  double fn ;
  int bmt ;
  double mfac, mpow ;
  double MeshMaxAspect ;
  double Gmax ;

  int MemOvershoot ;

  struct MGinfo bc ;
  double periodic[3] ;
  double pipe[3] ;
  int cathode ;			/* 'Leftmost' fieldline at z=0 */
  enum MGboundary zloworig ;	/* Original boundary condition */
  int cathodeorig ;		/* Original cathode */

  int Np2p ;			/* Defer to p2p when numaliveparticles < Np2p */
  int isp2p ;			/* Flag is set if p2p */
  double zmin ;			/* Not used */
  double zmax ;			/* Not used */

  /* Rest frame parameters */
  double GBo[3] ;
  double Go ;
  double (*ro)[3] ;
  gptsyncarray *sabuff ;

  /* Rotation in rest frame (before scaling) */
  gpttransform rotation ;

  /* Scaling parameters */
  int bscaling ;
  double scalesmall, scalelarge ;
  double scale[3] ;

  /* Beam statistics */
  unsigned int numparalive ;
  double beamavg[3], beamstd[3], beamstdmax ;
  double beamwavg[3], beamwstd[3], Q ;
  double beammin[3], beammax[3] ;

  /* Final Mesh */
  double mstart[3], mend[3] ;	// Mesh box dimensions
  unsigned int N[3], Nprev[3] ;
  double *meshlines[3] ;

  /* Charge */
  double *rho ;

  /* Potential */
  double Verror ;
  int bSolStartZero ;
  double *V ;
} ;


static void sc3dmesh_sim( double t, double *x, double *p, void *vinfo ) ;
static int restframeparams( struct sc3dmesh_info *info ) ;
static void transformtorestframe(struct sc3dmesh_info *info) ;
static int testrestframe_out(double t, double *dt, double *x, struct sc3dmesh_info *info) ;
static void rotaterestframe(struct sc3dmesh_info *info) ;
static void beamstatistics(struct sc3dmesh_info *info) ;
static void scalebeam(struct sc3dmesh_info *info) ;
static void meshbox(struct sc3dmesh_info *info) ;
static void meshlines(struct sc3dmesh_info *info) ;
static void readmeshfromfile(struct sc3dmesh_info *info, char *fn) ;
static void charge(struct sc3dmesh_info *info) ;
static void filtercharge(struct sc3dmesh_info *info, char *fn) ;
static void readcharge(struct sc3dmesh_info *info, char *fn) ;
static void potentialboundary(struct sc3dmesh_info *info) ;
static void printmeshstats(struct sc3dmesh_info *info) ;
static void interpolatetransform(struct sc3dmesh_info *info) ;

/* Helper routines */
static double Vbrick(double x, double y, double z, double sx, double sy, double sz) ;
static double Vsphere(double x, double y, double z, double sx, double sy, double sz) ;
static double gdfasqrt(double x) ;

/* ACAD output routines */
static void writeacadmesh(char *filename,struct sc3dmesh_info *info) ;
static void writeacadcharge(char *filename,struct sc3dmesh_info *info) ;
static void writeacadpotential(char *filename,struct sc3dmesh_info *info) ;

/* 3D mesh helper functions */
static void parmeshnotify(gptsyncarray *psa, gptsyncstatus stat, int n, struct sc3dmesh_info *info) ;
static int getindex(unsigned int *j, double *ro,double **meshlines,unsigned int *N) ;
static void linearinterpolate3D(double *V,unsigned int *j,double *r,double **meshlines,unsigned int *N, double *E) ;
static void quadraticinterpolate3D(double *V,unsigned int *j,double *r,double **meshlines,unsigned int *N, double *E) ;
static void cubicinterpolate3D(double *V,unsigned int *j,double *r,double **meshlines,unsigned int *N, double *E) ;


/* Initialization routine: Read command line parameters into info structure */
void spacecharge3Dmesh_init(gptinit *init)
{
  struct sc3dmesh_info *info ;
  int i ;
  unsigned int j ;
  char *option, *name ;

  /* Allocate info structure and set defaults */
  info = (struct sc3dmesh_info *)gptmalloc( sizeof(struct sc3dmesh_info) ) ;
  for(i=0 ; i<3 ; i++) info->Nc   [i] = 0 ;
  for(i=0 ; i<3 ; i++) info->Nprev[i] = 0 ;

  info->nstd  = 5 ;
  info->macc  = 0.3 ;
  info->fn    = 1.5 ;
  info->bmt   = 0 ;
  info->Gmax  = 2.0 ;

  info->mfac  = 1 ;
  info->mpow  = 1.0/3.0 ;

  info->MemOvershoot = 10 ;

  info->MeshMaxAspect = DBL_MAX ;

  info->Verror= 0.01 ;
  info->bSolStartZero = 0 ;

  info->scalesmall = 1e-4 ;
  info->scalelarge = 1e+4 ;

  for(j=0 ; j<3 ; j++) info->bc.low[j] = DirichletZero ;
  for(j=0 ; j<3 ; j++) info->bc.upp[j] = DirichletZero ;
  info->cathode = 0 ;

  info->bc.solvertype = SolverMGCG ;

  info->Np2p = 100 ;
  info->zmin  =-DBL_MAX ;
  info->zmax  = DBL_MAX ;

  for(i=1 ; i<=gptgetargnum(init) ; i++)
  {
    option = gptgetargstring(init,i) ;

    if( !strcmp(option,"MeshNtotal") )
    {
      info->bmt   = 1 ;
      if( i==gptgetargnum(init) || gptgetargtype(init,i+1)==GPTTYPE_STRING ) continue ;

      info->Nc[0] = gptgetargint(init,++i) ;
      info->Nc[1] = gptgetargint(init,++i) ;
      info->Nc[2] = gptgetargint(init,++i) ;
    } else

    if( !strcmp(option,"MeshNbunch") )
    {
      info->bmt   = 0 ;
      info->Nc[0] = gptgetargint(init,++i) ;
      info->Nc[1] = gptgetargint(init,++i) ;
      info->Nc[2] = gptgetargint(init,++i) ;
    } else

    if( !strcmp(option,"MeshNfac") )
    {
      info->mfac = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"MeshNpow") )
    {
      info->mpow = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"MemOvershoot") )
    {
      info->MemOvershoot = gptgetargint(init,++i) ;
    } else

    if( !strcmp(option,"MeshMaxAspect") )
    {
      info->MeshMaxAspect = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"MeshAdapt") )
    {
      info->fn    = gptgetargdouble(init,++i)+1 ;
    } else

    if( !strcmp(option,"MeshBoxSize") )
    {
      info->nstd  = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"MeshBoxAccuracy") )
    {
      info->macc  = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"RestMaxGamma") )
    {
      info->Gmax  = gptgetargdouble(init,++i) ;
    } else


    if( !strcmp(option,"BeamScale") )
    {
      info->scalesmall = gptgetargdouble(init,++i) ;
      info->scalelarge = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"Cathode") )
    {
      info->cathode = 1 ;	// Copied into cathodeorig below
    } else

    if( !strcmp(option,"Boundaries") )
    {
      name = gptgetargstring(init,++i) ;

      if( strlen(name)!=6 )
        gpterror( "All 6 boundary conditions need to be specified" ) ;

      info->bc.low[0] = (enum MGboundary)name[0] ;
      info->bc.upp[0] = (enum MGboundary)name[1] ;
      info->bc.low[1] = (enum MGboundary)name[2] ;
      info->bc.upp[1] = (enum MGboundary)name[3] ;
      info->bc.low[2] = (enum MGboundary)name[4] ;
      info->bc.upp[2] = (enum MGboundary)name[5] ;
    } else

    if( !strcmp(option,"BoundaryOpen") )
    {
      for(j=0 ; j<3 ; j++) info->bc.low[j] = Open ;
      for(j=0 ; j<3 ; j++) info->bc.upp[j] = Open ;
    } else

    if( !strcmp(option,"BoundaryXperiodic") )
    {
      info->bc.low[0] = info->bc.upp[0] = Periodic ;
      info->periodic[0] = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"BoundaryYperiodic") )
    {
      info->bc.low[1] = info->bc.upp[1] = Periodic ;
      info->periodic[1] = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"BoundaryZperiodic") )
    {
      info->bc.low[2] = info->bc.upp[2] = Periodic ;
      info->periodic[2] = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"BoundaryXYpipe") )
    {
      info->bc.low[0] = info->bc.upp[0] = Ellipse ;
      info->bc.low[1] = info->bc.upp[1] = Ellipse ;

      info->pipe[0] = gptgetargdouble(init,++i) ;
      info->pipe[1] = gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"SolverAcc") )
    {
      info->Verror= gptgetargdouble(init,++i) ;
    } else

    if( !strcmp(option,"SolverMethod") )
    {
      name = gptgetargstring(init,++i) ;

           if( !strcmp(name,"MG") )
        info->bc.solvertype = SolverMG ;
      else if( !strcmp(name,"MGCG") )
        info->bc.solvertype = SolverMGCG ;
      else if( !strcmp(name,"CG") )
        info->bc.solvertype = SolverCG ;
      else if( !strcmp(name,"SOR") )
        info->bc.solvertype = SolverSOR ;
      else
        gpterror( "Unknown solver method: %s\n", name ) ;
    } else

    if( !strcmp(option,"SolverStartZero") )
    {
      info->bSolStartZero = 1 ;
    } else

    if( !strcmp(option,"Np2p") )
    {
      info->Np2p = gptgetargint(init,++i) ;
    } else

    gpterror( "Unknown option: \"%s\"\n", option) ;
  }
  
  /* Store original cathode, and zlow boundary condition */
  info->cathodeorig = info->cathode ;
  info->zloworig = info->bc.low[2] ;

  /* Error conditions */
  if( info->nstd < 1 )
    gpterror("MeshBoxSize must be at least 1\n") ;
  if( info->fn < 1 )
    gpterror("MeshAdapt must be positive\n") ;

  /* Memory management */
  info->ro   = NULL ;
  for(j=0 ; j<3 ; j++)
    info->meshlines[j] = NULL ;
  info->rho = NULL ;
  info->V   = NULL ;

  gptmemblockinit(&info->bc.mem,65536,info->MemOvershoot) ;
  info->sabuff=NULL ;

   /* Register simulation function */
  odeaddfprfunction( ODEFNC_USR, sc3dmesh_sim,info ) ;
  odeaddoutfunction( ODEFNC_USR, (odeoutfnc)testrestframe_out, info ) ;
}


/* Simulation function: calculate E and B fields due to 3D space-charge with Rostock multigrid poisson solver */
static void sc3dmesh_sim( double t, double *x, double *p, void *vinfo )
{
  struct sc3dmesh_info *info = (struct sc3dmesh_info *)vinfo ;
  int i, numalivepar ;
  double Verror ;
  enum MGreturn rc ;

/* Get number of paticles in simulation */
  numalivepar = 0 ;
//#pragma omp parallel for reduction(+:numalivepar) if(GPTOMPSIMP(numpar))
  for(i=0 ; i<numpar ; i++) if( pars[i].alive ) numalivepar++ ;

/* Special cases: return on 1 particle, Point-to-point if N<Np2p */
  info->isp2p = 1 ;
  if( numalivepar<=1) return ;
  if( numalivepar<info->Np2p )
  {
    sc3Dp2p(t,x,p,info->zmin,info->zmax, 0,1) ;
    return ;
  }
  info->isp2p = 0 ;

/* Calculate rest frame parameters, return immediately when no charge present */
  if( restframeparams( info )) return ;

/* Transform to rest frame */
  transformtorestframe( info ) ;
  rotaterestframe(info) ;

/* Calculate beam statistics in rest frame */
  beamstatistics( info ) ;
  if( info->numparalive<=1 ) return ;
  if( info->beamstdmax==0 ) return ;
  scalebeam( info ) ;

/* Calculate meshlines */
#ifndef SCDEBUG_READMESHFROMFILE
  meshbox( info ) ;
  meshlines( info ) ;
#else
  readmeshfromfile( info, SCDEBUG_READMESHFROMFILE ) ;
#endif

#ifdef SCDEBUG_MESHSTATS
  printmeshstats(info) ;
#endif

/* Calculate density on meshpoints */
  info->rho = (double *)gptrealloc(info->rho,info->N[0]*info->N[1]*info->N[2]*sizeof(double)) ;
  info->V   = (double *)gptrealloc(info->V  ,info->N[0]*info->N[1]*info->N[2]*sizeof(double)) ;

  if( !info->sabuff )
//  info->sabuff = gptcreatesyncarray(gptparsyncset,sizeof(struct parmeshinfo),parmeshnotify,info) ;
    info->sabuff = gptcreatesyncarray(gptparsyncset,sizeof(struct parmeshinfo),(void (*)(gptsyncarray *,gptsyncstatus,int,void *))parmeshnotify,info) ;

  charge( info ) ;

#ifdef WAVELET_COMMAND
  filtercharge( info, WAVELET_INFILE ) ;
  if( system(WAVELET_COMMAND)==-1 )
    gpterror("Error running: " WAVELET_COMMAND "\n" ) ;
  readcharge( info, WAVELET_OUTFILE ) ;
#endif

  potentialboundary( info ) ;

/* Solve Poisson's equaitions with the Rostock multigrid solver */
  Verror = info->Verror ;

  do
  {
    gptmemblockresize(&info->bc.mem) ;

    switch( rc=(enum MGreturn)multigridpoisson283((int *)info->N,info->meshlines,info->rho,info->V,&Verror,&info->bc) )
    {
      case NoMemory:
#ifdef SCDEBUG_MEMORY
        fprintf( stderr, "Resizing block memory: %d bytes out of %d needed.\n",
          info->bc.mem.memnext-info->bc.mem.memstart,
          info->bc.mem.memend -info->bc.mem.memstart ) ;
        fflush(stdout) ;
#endif
      break ;

      case NoConvergence:
        gptwarning( "No convergence in Poisson solver" ) ;
        break ;

      case InvalidBoundaries:
        gpterror( "Invalid boundary conditions" ) ;

      case NotImplementedYet:
        gpterror( "Requested feature is not implemented yet\n" ) ;

      case NoGrids:
        gpterror( "No coarser grids in Multigrid Poisson solver\n" ) ;

	  case NoError:
        break ;
    }
  } while( rc==NoMemory ) ;

/* Calculate E-field from potential and transform to lab frame */
  interpolatetransform(info) ;

#if defined(SCDEBUG_MESHLINES) || defined(SCDEBUG_RESTSTATS) || defined(SCDEBUG_MESHSTATS)
  fflush(stderr) ;
  fflush(stdout) ;
#endif

/* ACAD Debug scripts: Write mesh, charge and potential */
#ifdef SCDEBUG_MESHFILE
  writeacadmesh(SCDEBUG_MESHFILE, info) ;
#endif

#ifdef SCDEBUG_CHARFILE
  writeacadcharge(SCDEBUG_CHARFILE, info) ;
#endif

#ifdef SCDEBUG_POTEFILE
  writeacadpotential(SCDEBUG_POTEFILE, info) ;
#endif
}


/* Calculate rest frame parameters */
static int restframeparams(struct sc3dmesh_info *info)
{
  double sn, snG ;
  double snGB[3], Bo[3] ;
  int i ;
  unsigned int j ;
  static int sn0warning=0 ;

  /* Calculate sumG and sumGB, weighted by nmacro */
  sn  = 0 ;
  snG = 0 ;
  for(j=0 ; j<3 ; j++) snGB[j] = 0 ;
  for(i=0 ; i<numpar ; i++) if( pars[i].alive )
  {
    sn  += pars[i].n ;
    snG += pars[i].n * pars[i].G ;
    for(j=0 ; j<3 ; j++) snGB[j] += pars[i].n * pars[i].GBr[j] ;
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
  for(j=0 ; j<3 ; j++) Bo[j] = snGB[j]/snG ; 
  for(j=0 ; j<3 ; j++) info->GBo[j] = Bo[j]/sqrt(1-gptVECSQR(Bo)) ;
  info->Go = sqrt(1 + gptVECSQR(info->GBo) ) ;

#ifdef SCDEBUG_RESTSTATS
  fprintf( stderr, "Rest frame GBo=(%f,%f,%f), Go=%f\n", info->GBo[0], info->GBo[1], info->GBo[2], info->Go ) ;
#endif

  return(0) ;
}


/* Transform particles to rest frame */
static void transformtorestframe(struct sc3dmesh_info *info)
{
  int i ;
  unsigned int j ;

  info->ro = (double (*)[3])gptrealloc(info->ro,numpar*sizeof(*info->ro) ) ;

  for(i=0 ; i<numpar; i++ ) if( pars[i].alive )
  {
    double alpha = gptVECINP(pars[i].Wr,info->GBo)/(info->Go + 1) ;
    for(j=0 ; j<3 ; j++ ) info->ro[i][j] = pars[i].Wr[j] + alpha*info->GBo[j] ;
  }
}


/* Test rest frame velocities. Called using OUT interface at beginning of timestep */
static int testrestframe_out(double t, double *dt, double *x, struct sc3dmesh_info *info)
{
  int i ;
  int nnewmaxspeed ;
  static int noldmaxspeed=0 ;
  double Gio ;

  /* No rest frame check in p2p mode */
  if( info->isp2p ) return(0) ;

  /* Calculate number of speed offenders in rest frame */
  nnewmaxspeed=0 ;
  for(i=0 ; i<numpar; i++ ) if( pars[i].alive )
  {
    Gio = pars[i].G*info->Go-gptVECINP(pars[i].GBr,info->GBo) ;

    if( Gio>info->Gmax ) nnewmaxspeed++ ; 
  }

  /* Print warning in verbose mode or when the number of speed offenders varies by 20% */
  if( (verbose && nnewmaxspeed) || nnewmaxspeed<noldmaxspeed*10/12 || nnewmaxspeed>noldmaxspeed*12/10 )
  {
    gptwarning("Spacecharge3Dmesh: %d particles with Lorentz factor in rest frame over %f at t=%e\n", nnewmaxspeed, info->Gmax, t) ;
    noldmaxspeed = nnewmaxspeed ;
  }

  return(0) ;
}

/* Calculate and transform beam rotation in rest frame */
static void rotaterestframe(struct sc3dmesh_info *info)
{
  unsigned int i,j ;

  if( info->cathode )
  {
    // No rotation in case of a cathode
    createidentitytransform(&info->rotation) ;
  } else
  {
    createdirectionaltransform(&info->rotation,info->GBo) ;
	for(i=0 ; i<numpar; i++ ) if( pars[i].alive )
	{
	  double r[3] ;
	  gpttoUCS(&info->rotation,info->ro[i],r) ;
	  for(j=0 ; j<3 ; j++) info->ro[i][j] = r[j] ;
	}
  }
}

/* Calculate beam statistics in rest frame */
static void beamstatistics(struct sc3dmesh_info *info)
{
  int i ;
  unsigned int j, nrj ;
  double srj, srj2 ;
  double wrj, wsrj, wsrj2 ;

  info->beamstdmax = 0 ;

  /* Loop over x, y, z using j */
  for(j=0 ; j<3 ; j++)
  {
    /* Determine min and max */
    info->beammin[j] = +DBL_MAX ;
    info->beammax[j] = -DBL_MAX ;
    for(i=0 ; i<numpar ; i++) if( pars[i].alive )
    {
      if( info->ro[i][j] < info->beammin[j] ) info->beammin[j] = info->ro[i][j] ;
      if( info->ro[i][j] > info->beammax[j] ) info->beammax[j] = info->ro[i][j] ;
    }

    /* Determine avg and std, without nmacro */
    nrj = 0 ;
	srj = srj2 = 0 ;
//#pragma omp parallel for reduction(+:nrj,srj,srj2) if(GPTOMPSIMP(numpar))
    for(i=0 ; i<numpar ; i++) if( pars[i].alive )
    {
      nrj++ ;
      srj  += info->ro[i][j] ;
      srj2 += info->ro[i][j]*info->ro[i][j] ;
    }
    info->numparalive = nrj ;
    if( nrj==0 ) return ;

    info->beamavg[j] = srj/nrj ;
    info->beamstd[j] = gdfasqrt((srj2-srj*srj/nrj)/nrj) ;

    /* Determine avg and std, with nmacro */
    wrj = wsrj = wsrj2 = 0 ;
//#pragma omp parallel for reduction(+:wrj,wsrj,wsrj2) if(GPTOMPSIMP(numpar))
    for(i=0 ; i<numpar ; i++) if( pars[i].alive )
    {
	  double w = fabs(pars[i].n*pars[i].q) ;
      wrj   += w ;
      wsrj  += w * info->ro[i][j] ;
      wsrj2 += w * info->ro[i][j]*info->ro[i][j];
    }

    if( wrj!=0 )
    {
      info->beamwavg[j] = wsrj/wrj ;
      info->beamwstd[j] = gdfasqrt((wsrj2-wsrj*wsrj/wrj)/wrj) ;
    } else
    {
      info->beamwavg[j] = srj/nrj ;
      info->beamwstd[j] = gdfasqrt((srj2-srj*srj/nrj)/nrj) ;
    }
    info->Q = wrj ;

    if( info->beamstd[j] > info->beamstdmax )
      info->beamstdmax = info->beamstd[j] ;
  }

  /* Bye bye with 1 particle */
  if( info->numparalive<=1 ) return ;

  /* Fix zero-sized dimension */
  for(j=0 ; j<3 ; j++)
    if( info->beamstd[j]==0 )
    {
      gptwarning("Infinitely small %c-dimension\n", j+'x') ;
      info->beamstd[j] = info->beamstdmax * info->scalesmall ;
    }

#ifdef SCDEBUG_RESTSTATS
  for(j=0 ; j<3 ; j++)
    fprintf( stderr, "pavg%c=%f, pstd%c=%f. ", 'x'+j, info->beamavg[j], 'x'+j,info->beamstd[j] ) ;
  fprintf( stderr, "\n") ;
  for(j=0 ; j<3 ; j++)
    fprintf( stderr, "wavg%c=%f, wstd%c=%f. ", 'x'+j, info->beamwavg[j], 'x'+j,info->beamwstd[j] ) ;
  fprintf( stderr, "\n") ;
  fprintf( stderr, "Q: %e\n", info->Q) ;
#endif
}


/* Scale beam when one dimension is much smaller/larger than the other two */
static void scalebeam(struct sc3dmesh_info *info)
{
  int i ;
  unsigned int j ;
  double scale,tmp ;

  info->bscaling = 0 ;

  for(j=0 ; j<3 ; j++)
  {
    info->scale[j] = 1 ;

    scale = info->beamstd[j]/info->beamstd[(j+1)%3] ;
    if( (tmp=info->beamstd[j]/info->beamstd[(j+2)%3]) > scale ) scale=tmp ;
    if( scale<info->scalesmall )
    {
      scale=info->scalesmall/scale ;

      if( verbose )
        gptwarning( "Scaling %c coordinates by %lf\n", 'x'+j, scale ) ;

      for(i=0 ; i<numpar ; i++) if( pars[i].alive )
      {
         info->ro[i][j] *= scale ;
      }
      info->scale[j] *= scale ;
      info->bscaling = 1 ;
    }

    scale = info->beamstd[j]/info->beamstd[(j+1)%3] ;
    if( (tmp=info->beamstd[j]/info->beamstd[(j+2)%3]) < scale ) scale=tmp ;
    if( scale>info->scalelarge )
    {
      scale = info->scalelarge/scale ;
      gptwarning( "Scaling %c coordinates by %lf\n", 'x'+j, scale ) ;
      for(i=0 ; i<numpar ; i++) if( pars[i].alive )
      {
         info->ro[i][j] *= scale ;
      }
      info->scale[j] *= scale ;
      info->bscaling = 1 ;
    }
  }

  /* Recalculate */
  if( info->bscaling )
    beamstatistics(info) ;
}


/* Bounding box based on mirror field being close enough to zero */
#define BOXSTEPSIZE 1.1
static void meshbox(struct sc3dmesh_info *info)
{
  double Vorigin ;
  unsigned int j ;
  double avg ;
  double newbox[3] ;
  double sx, sy, sz ;

  sx = info->beamwstd[0] ;
  sy = info->beamwstd[1] ;
  sz = info->beamwstd[2] ;

  newbox[0] = info->nstd*sx ;
  newbox[1] = info->nstd*sy ;
  newbox[2] = info->nstd*sz ;

  Vorigin = Vbrick(0,0,0,sx,sy,sz) ;


// Mirror bunch, derivative in cross direction
  while( Vbrick(2*newbox[0],0,0,sx,sy,sz)-Vbrick(2*newbox[0],sy,sz,sx,sy,sz) > info->macc*(Vorigin-Vbrick(0,sy,sz,sx,sy,sz)) )
    newbox[0] *= BOXSTEPSIZE ;
  while( Vbrick(0,2*newbox[1],0,sx,sy,sz)-Vbrick(sx,2*newbox[1],sz,sx,sy,sz) > info->macc*(Vorigin-Vbrick(sx,0,sz,sx,sy,sz)) )
    newbox[1] *= BOXSTEPSIZE ;
  while( Vbrick(0,0,2*newbox[2],sx,sy,sz)-Vbrick(sx,sy,2*newbox[2],sx,sy,sz) > info->macc*(Vorigin-Vbrick(sx,sy,0,sx,sy,sz)) )
    newbox[2] *= BOXSTEPSIZE ;

// Mirror bunch, derivative
//  while( Vbrick(2*newbox[0],0,0,sx,sy,sz)-Vbrick(2*newbox[0]+sx,0,0,sx,sy,sz) > info->macc*(Vorigin-Vbrick(sx,0,0,sx,sy,sz)) )
//   newbox[0] *= BOXSTEPSIZE ;
//  while( Vbrick(0,2*newbox[1],0,sx,sy,sz)-Vbrick(0,2*newbox[1]+sy,0,sx,sy,sz) > info->macc*(Vorigin-Vbrick(0,sy,0,sx,sy,sz)) )
//    newbox[1] *= BOXSTEPSIZE ;
//  while( Vbrick(0,0,2*newbox[2],sx,sy,sz)-Vbrick(0,0,2*newbox[2]+sz,sx,sy,sz) > info->macc*(Vorigin-Vbrick(0,0,sz,sx,sy,sz)) )
//    newbox[2] *= BOXSTEPSIZE ;


#ifdef SCDEBUG_MESHLINES
  fprintf( stderr, "Bounding box: (%.1f,%.1f,%.1f) standard deviations\n", newbox[0]/info->beamstd[0], newbox[1]/info->beamstd[1], newbox[2]/info->beamstd[2] ) ;
#endif

  for(j=0 ; j<3 ; j++)
  {
    /* Determine meshline positions */
    avg = info->beamavg[j] ;

    info->mstart[j] = avg-newbox[j] ;
    info->mend[j]   = avg+newbox[j] ;

    /* Overwrite for periodic and ellipse boundary conditions */
    if( info->bc.low[j]==Periodic ) info->mstart[j] = avg-info->periodic[j]/2 ;
    if( info->bc.upp[j]==Periodic ) info->mend[j]   = avg+info->periodic[j]/2 ;
    if( info->bc.low[j]==Ellipse  ) info->mstart[j] = avg-info->pipe[j]/2 ;
    if( info->bc.upp[j]==Ellipse  ) info->mend[j]   = avg+info->pipe[j]/2 ;
  }

  /* Optional cathode at z=0 */
  if( info->cathodeorig )
  {
    if( info->mstart[2] < 0 )
    {
      info->mstart[2] = 0 ; // -DBL_EPSILON ;

      if( verbose && info->bc.low[2]!=DirichletZero )
        gptwarning( "Image charges at cathode enabled" ) ;

      info->bc.low[2] = DirichletZero ;
	  info->cathode = 1 ;
    } else
    {
      if( verbose && info->bc.low[2]!=info->zloworig )
        gptwarning( "Image charges at cathode disabled, continuing with %c boundary condition", info->zloworig ) ;

      info->bc.low[2] = info->zloworig ;
	  info->cathode = 0 ;
    }
  }
}

/* Calculate mesh lines positions */
#define FNREDUCE 0.8
static void meshlines(struct sc3dmesh_info *info)
{
  int i,ll,lr,Nj,Njold ; //TODO: ll,lr etc. should be unsigned
  unsigned int j,k,bCont ;
  double *meshlinesj, *meshtemp, mstart, mend ;
  double fnll, fnlr, fni ;
  double mmin, mmax, mnew, mmove ;
  double fns ;
  double hmin, hmax, dh, aspect ;

  unsigned int Nhis, n ;
  double *his ;
  double hismin[3], hismax[3] ;

  /* Initialization: his */
#define HISREDUCE 0.25
  Nhis = 1+(int)(sqrt((double)info->numparalive)*HISREDUCE) ;
  his = (double *)gptmalloc((Nhis+1)*sizeof(*his)) ;

  for(j=0 ; j<3 ; j++) hismin[j] = (info->beammin[j] >= info->mstart[j] ? info->beammin[j] : info->mstart[j] ) ;
  for(j=0 ; j<3 ; j++) hismax[j] = (info->beammax[j] <= info->mend  [j] ? info->beammax[j] : info->mend  [j] ) ;

  /* Main loop over x, y, z using j */
  for( j=0 ; j<3  ; j++ )
  {
    mstart = info->mstart[j] ;
    mend   = info->mend[j] ;

    /* Normalized charge density histogram */
    for(k=0 ; k<=Nhis ; k++) his[k] = 0 ;
    for(i=0 ; i<numpar ; i++) if( pars[i].alive )
    {
      k = (int)(Nhis*(info->ro[i][j]-hismin[j])/(hismax[j]-hismin[j])) ;
//    if( k<0     ) k=0 ;
//	  if( k>=Nhis ) k=Nhis-1 ;
	  if( k<0     ) continue ;
	  if( k>=Nhis ) continue ;
      his[1+k]++ ;
    }

    for(k=0 ; k<=Nhis ; k++) his[k] = pow((double)his[k],1.0/3.0) ;
	for(k=1 ; k<=Nhis ; k++) his[k] += 1.0/(double)Nhis ; // Enforce increasing integral by adding a small constant value
    for(k=1 ; k<=Nhis ; k++) his[k] += his[k-1] ;
    for(k=0 ; k<=Nhis ; k++) his[k] /= his[Nhis] ;

    Nj         = info->Nc[j] ;
    if( Nj==0 ) Nj = (int)(0.5+info->mfac*pow(info->numparalive,info->mpow)) ;
    if( Nj <4 ) Nj = 4 ;

    meshtemp   = NULL ;

    /* Reduced fn with FNREDUCE safety margin */
    fns = 1+FNREDUCE*(info->fn-1) ;

    /* Double !!SLOW!! iterative procedure */
    ll=lr=1 ; /* Start with 1 additional meshline on both sides */
    do
    {
      /* Memory management */
      meshlinesj = info->meshlines[j] = (double *)gptrealloc(info->meshlines[j],(Nj+0)*sizeof(*meshlinesj)) ;
      meshtemp   =                      (double *)gptrealloc(meshtemp,          (Nj+1)*sizeof(*meshtemp  )) ;

      /* his: Meshlines  following beam density distribution */
      n=1 ;
      for(i=ll ; i<Nj-lr ; i++)
      {
//      double f = (double)(i-ll)/(Nj-ll-lr-1) ;	// Equidistant, endpoints INCLUDED
        double f = (double)(i+0.5-ll)/(Nj-ll-lr) ;	// Equidistant, endpoints NOT included

        while( his[n]<f && n!=Nhis ) n++ ;
        meshtemp[i] = hismin[j]+((f-his[n-1])/(his[n]-his[n-1])+n-1)*(hismax[j]-hismin[j])/Nhis ;
      }

      /* Additional left meshlines exponentially decreasing in size */
      fnll = pow(fns,(double)ll) ;
      fni=1 ;
      for(i=0 ; i<ll ; i++) 
      {
        meshtemp[     i] = mstart + (meshtemp[     ll]-mstart) * fnll*(1-fni)/(fnll-1) ;
        fni /= fns ;
      }

      /* Additional right meshlines exponentially increasing in size */
      fnlr = pow(fns,(double)lr) ;
      fni=1 ;
      for(i=0 ; i<lr ; i++) 
      {
        meshtemp[Nj-1-i] = mend   + (meshtemp[Nj-1-lr]-mend  ) * fnlr*(1-fni)/(fnlr-1) ;
        fni /= fns ;
      }

      /* First guess for meshlines */
      if(info->fn==1 ) for(i=0 ; i<Nj ; i++) meshtemp[i] = mstart+i*(mend-mstart)/(Nj-1) ;
      for(i=0 ; i<Nj ; i++) meshlinesj[i] = meshtemp[i] ; 

//      fprintf( stderr, "%c in (%lf,%lf)\n", 'x'+j, mstart, mend ) ;
//      for(i=0 ; i<Nj ; i++) fprintf( stderr, "% .3lf ", 1000*meshtemp[i] ) ; fprintf( stderr, "\n" ) ;

      /* Adapt meshlines positions to enfore maximum size differenence between neighbouring lines */
      do
      {
        k=0 ;
        mmove=0 ;

        for(i=1 ; i<Nj-1 ; i++)
        {
          mmin = (meshlinesj[i+1]+meshlinesj[i-1]*info->fn)/(1+info->fn) ;
          mmax = (meshlinesj[i-1]+meshlinesj[i+1]*info->fn)/(1+info->fn) ;

          if( meshtemp[i]<mmin ) mnew = mmin ; else
          if( meshtemp[i]>mmax ) mnew = mmax ; else
          mnew = meshtemp[i] ;

          mmove += fabs(mnew-meshlinesj[i]) ;
          meshlinesj[i] = mnew ;
        }

//        for(i=0 ; i<Nj ; i++) fprintf( stderr, "% .3lf ", 1000*meshlinesj[i] ) ; fprintf( stderr, "\n" ) ;

        k++ ;
      } while( mmove > MESHPREC*Nj*(mend-mstart) && k<MAXMESHIT) ;
      if( k==MAXMESHIT )
        gptwarning( "Non-converging mesh optimization. avg%c=%f, std%c=%f. %d lines between (%f,%f).", 'x'+j,info->beamavg[j], 'x'+j,info->beamstd[j], Nj, mstart, mend ) ;

      /* Calculate aspect ratio */
      hmin = DBL_MAX ;
      hmax = 0 ;
      for(i=0 ; i<Nj-1 ; i++)
      {
        dh = meshlinesj[i+1]-meshlinesj[i] ;
        if( dh<hmin ) hmin=dh ;
        if( dh>hmax ) hmax=dh ;
      }
      aspect = hmax/hmin ;

      bCont=0 ;
      Njold = Nj ;
      if( (meshlinesj[ll+1]-meshlinesj[ll]) < (meshlinesj[ll]-meshlinesj[ll-1]) )
      {
        ll++ ;
        if( !info->bmt ) Nj++ ;

        bCont=1 ;
      }

      if( (meshlinesj[Njold-1-lr]-meshlinesj[Njold-1-lr-1]) < (meshlinesj[Njold-lr]-meshlinesj[Njold-1-lr]) )
      {
        lr++ ;
        if( !info->bmt ) Nj++ ;

        bCont=1 ;
      }

      if( !bCont && aspect>info->MeshMaxAspect )
      {
        fns = 1+FNREDUCE*(fns-1) ;
        bCont=1 ;
#ifdef SCDEBUG_MESHLINES
        fprintf( stderr, "A%c=%lf. Reducing fn to %lf\n", 'x'+j, aspect, fns) ;
#endif
      }

    } while( bCont && Nj-ll-lr>1 ) ;

   /* Store final number of meshlines */
   info->N[j] = Nj ;

    /* Debug: Print statistics */
#ifdef SCDEBUG_MESHLINES
    fprintf( stderr, "A%c=%lf. %d lines between (%f,%f). ", 'x'+j, aspect, Nj, mstart, mend ) ;
    fprintf( stderr, "Created %d+%d additional meshlines with fn=%lf. \n", ll, lr, fns-1 ) ;
#endif

    /* Consistency check */
    for(i=0 ; i<info->N[j]-1 ; i++)
      if( ! (info->meshlines[j][i+1] > info->meshlines[j][i]) )
      {
        gptwarning("Beammin[%d]=%le, mstart[%d]=%le, hismin[%d]=%le\n", j,info->beammin[j], j,info->mstart[j], j,hismin[j] ) ;
        gptwarning("Beammax[%d]=%le, mend  [%d]=%le, hismax[%d]=%le\n", j,info->beammax[j], j,info->mend  [j], j,hismax[j] ) ;

        for(i=0 ; i<info->N[j] ; i++)
          gptwarning("Mesh[%d]=%le\n", i, info->meshlines[j][i] ) ;

        gpterror("Internal mesh inconsistency: Meshline positions (%c) are not monotonically increasing", j+'x') ;
      }
  }

  /* Memory management */
  gptfree(meshtemp) ;
  gptfree(his) ;
}


/* Read meshlines from file */
static void readmeshfromfile(struct sc3dmesh_info *info, char *fn)
{
  FILE *fp ;
  unsigned int i, j ;

  gptwarning( "Reading meshlines from file \"%s\"...", fn ) ;

  if( (fp=fopen(fn,"r"))==NULL )
    gpterror("%s: %s", fn, strerror(errno) ) ;

  fscanf(fp,"%d%d%d", &info->N[0], &info->N[1], &info->N[2]) ;

  for(j=0 ; j<3 ; j++)
  {
    info->meshlines[j] = (double *)gptrealloc(info->meshlines[j],info->N[j]*sizeof(double)) ;

    for(i=0 ; i<info->N[j] ; i++)
      if( fscanf(fp,"%lf",&info->meshlines[j][i])!=1 )
        gpterror("%s: %s\n", fn, strerror(errno)) ;

    info->mstart[j] = info->meshlines[j][0] ;
    info->mend  [j] = info->meshlines[j][info->N[j]-1] ;
  }

  fclose(fp) ;

#ifdef SCDEBUG_MESHLINES
    for(j=0 ; j<3 ; j++)
      printf( "%c: %d lines between (%f,%f)\n", 'x'+j, info->N[j], info->mstart[j], info->mend[j] ) ;
#endif
}


/* Calculate charge */
static void charge(struct sc3dmesh_info *info)
{
  int i ;
  unsigned int *j,k, N1,N2,N3,N4,N5,N6,N7,N8 ;
  double ft,fu,fv, gt,gu,gv, f1,f2,f3,f4,f5,f6,f7,f8, rhoi ;
  struct parmeshinfo *pm ;

/* Set density to zero */
  for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++) info->rho[i] = 0 ;

/* Fix periodic boundaries */
  for(k=0 ; k<3 ; k++)
  {
    if( info->bc.low[k]==Periodic ) for(i=0 ; i<numpar; i++ ) if( info->ro[i][k]<info->meshlines[k][           0] ) info->ro[i][k] += (info->meshlines[k][info->N[k]-1]-info->meshlines[k][0]) ;
    if( info->bc.upp[k]==Periodic ) for(i=0 ; i<numpar; i++ ) if( info->ro[i][k]<info->meshlines[k][info->N[k]-1] ) info->ro[i][k] -= (info->meshlines[k][info->N[k]-1]-info->meshlines[k][0]) ;
  }

/* Calculate the total charge at the corners of every meshbox */
  pm = (struct parmeshinfo *)info->sabuff->memory ;
  for(i=0 ; i<numpar; i++ ) if( pars[i].alive )
  {
    j = pm[i].jprev ;
    if( getindex(j,info->ro[i],info->meshlines, info->N)!=0 )
    {
      N1 = (j[2]*info->N[1]+j[1])*info->N[0]+j[0] ;
      N2 = N1+1 ;
      N3 = N1+info->N[0]+1 ;
      N4 = N1+info->N[0] ;
      N5 = N1+info->N[0]*info->N[1] ;
      N6 = N1+info->N[0]*info->N[1]+1 ;
      N7 = N1+info->N[0]*(info->N[1]+1)+1 ;
      N8 = N1+info->N[0]*(info->N[1]+1) ;

      ft= (info->ro[i][0]-info->meshlines[0][j[0]])/(info->meshlines[0][j[0]+1] - info->meshlines[0][j[0]]) ;
      fu= (info->ro[i][1]-info->meshlines[1][j[1]])/(info->meshlines[1][j[1]+1] - info->meshlines[1][j[1]]) ;
      fv= (info->ro[i][2]-info->meshlines[2][j[2]])/(info->meshlines[2][j[2]+1] - info->meshlines[2][j[2]]) ;
      gt = 1-ft ;
      gu = 1-fu ;
      gv = 1-fv ;

      f1=gt*gu*gv ; f2=ft*gu*gv ; f3=ft*fu*gv ; f4=gt*fu*gv ; 
      f5=gt*gu*fv ; f6=ft*gu*fv ; f7=ft*fu*fv ; f8=gt*fu*fv ; 

      rhoi = pars[i].n*pars[i].q / gpt_eps0 ;

// Max only
//      {
//        double f0 = 0 ;
//        if( f1>f0 ) f0=f1 ;
//        if( f2>f0 ) f0=f2 ;
//        if( f3>f0 ) f0=f3 ;
//        if( f4>f0 ) f0=f4 ;
//        if( f5>f0 ) f0=f5 ;
//        if( f6>f0 ) f0=f6 ;
//        if( f7>f0 ) f0=f7 ;
//        if( f8>f0 ) f0=f8 ;

//        if( f1==f0 ) f1=1 ; else f1=0 ;
//        if( f2==f0 ) f2=1 ; else f2=0 ;
//        if( f3==f0 ) f3=1 ; else f3=0 ;
//        if( f4==f0 ) f4=1 ; else f4=0 ;
//        if( f5==f0 ) f5=1 ; else f5=0 ;
//        if( f6==f0 ) f6=1 ; else f6=0 ;
//        if( f7==f0 ) f7=1 ; else f7=0 ;
//        if( f8==f0 ) f8=1 ; else f8=0 ;
//      }

      info->rho[N1] += f1*rhoi ;
      info->rho[N2] += f2*rhoi ;
      info->rho[N3] += f3*rhoi ;
      info->rho[N4] += f4*rhoi ;
      info->rho[N5] += f5*rhoi ;
      info->rho[N6] += f6*rhoi ;
      info->rho[N7] += f7*rhoi ;
      info->rho[N8] += f8*rhoi ;
    } else
      j[0]=-1 ;
  }

  /* Don't start from previous solution when mesh-size differes */
  if( info->bSolStartZero || info->Nprev[0]!=info->N[0] || info->Nprev[1]!=info->N[1] || info->Nprev[2]!=info->N[2] )
  {
#ifdef SCDEBUG_MESHLINES
    printf( "Number of meshlines changed from (%d,%d,%d) to (%d,%d,%d). Setting potiential to zero\n", info->Nprev[0], info->Nprev[1], info->Nprev[2], info->N[0], info->N[1], info->N[2] ) ;
#endif
    for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++) info->V[i] = 0 ;
  }
  for(i=0 ; i<3 ; i++) info->Nprev[i] = info->N[i] ;
}


static void swap(void *x,int n)
{
  char tmp, *arr = (char *)x ;
  int i ;

  for(i=0 ; i<n/2 ; i++)
  {
    tmp=arr[i] ; arr[i]=arr[n-i-1] ; arr[n-i-1]=tmp ;
  }
}

#define FITSBLOCK 2880
/* Write charge density to FITS file */
static void filtercharge(struct sc3dmesh_info *info, char *fn)
{
  FILE *fp ;
  int i, Ntotal ;
  float *rhofloat ;
  off64_t pos ;

  if( info->N[0]!=30 || info->N[1]!=30 || info->N[2]!=30 )
    gpterror("Only 30 meshlines supported by demo version\n" ) ;

  Ntotal = info->N[0]*info->N[1]*info->N[2] ;

  if( (fp=fopen(fn,"wb"))==NULL )
    gpterror("%s: %s",fn,strerror(errno)) ;
  fprintf(fp, "SIMPLE  =                    T / THIS IS REGULAR FITS                           " ) ;
  fprintf(fp, "BITPIX  =                  -32 /                                                " ) ;
  fprintf(fp, "NAXIS   =                    3 /                                                " ) ;
  fprintf(fp, "NAXIS1  =                   30 /                                                " ) ;
  fprintf(fp, "NAXIS2  =                   30 /                                                " ) ;
  fprintf(fp, "NAXIS3  =                   30 /                                                " ) ;
  fprintf(fp, "HISTORY GPT charge density filter test                                          " ) ;
  fprintf(fp, "END                                                                             " ) ;

  pos = ftello(fp) ;
  for(i=0 ; i<(FITSBLOCK-(pos%FITSBLOCK))%FITSBLOCK ; i++) fprintf(fp," " ) ;

  rhofloat = (float *)gptmalloc(Ntotal*sizeof(*rhofloat)) ;
  for(i=0 ; i<Ntotal ; i++) rhofloat[i] = (float)(info->rho[i]) ;

  for(i=0 ; i<Ntotal ; i++) rhofloat[i] /= (float)(pars[0].q*pars[0].n / gpt_eps0) ;

  for(i=0 ; i<Ntotal ; i++) swap(&(((char *)rhofloat)[i*sizeof(*rhofloat)]),sizeof(*rhofloat)) ;
  fwrite(rhofloat,Ntotal,sizeof(*rhofloat),fp) ;

  pos = ftello(fp) ;
  for(i=0 ; i<(FITSBLOCK-(pos%FITSBLOCK))%FITSBLOCK ; i++) fputc('\0',fp) ;

  fclose(fp) ;

  gptfree(rhofloat) ;
}


static void readcharge(struct sc3dmesh_info *info, char *fn)
{
  FILE *fp ;
  float *rhofloat ;
  int i, Ntotal ;

  gptwarning( "Continue from stored charge...\n" ) ;

  Ntotal = info->N[0]*info->N[1]*info->N[2] ;

  if( (fp=fopen(fn,"rb"))==NULL )
    gpterror("%s: %s",fn,strerror(errno)) ;

  rhofloat = (float *)gptmalloc(Ntotal*sizeof(*rhofloat)) ;

  fseek(fp,FITSBLOCK,SEEK_SET) ;
  fread(rhofloat,Ntotal,sizeof(*rhofloat),fp) ;
  for(i=0 ; i<Ntotal ; i++) swap(&(((char *)rhofloat)[i*sizeof(*rhofloat)]),sizeof(*rhofloat)) ;

  for(i=0 ; i<Ntotal ; i++) info->rho[i] = rhofloat[i]*pars[0].q*pars[0].n / gpt_eps0 ;

  gptfree(rhofloat) ;
  fclose(fp) ;
}


#ifdef SC3DMESH_POTENTIALSPHERE
#define VFUNCTION Vsphere
#endif
#ifdef SC3DMESH_POTENTIALBRICK
#define VFUNCTION Vbrick
#endif

/* Estimate the potential at the boundary for Dirichlet boundary conditions */
static void potentialboundary(struct sc3dmesh_info *info)
{
  unsigned int i, j, k, ip, Nx, Ny, Nz ;
  double xx, yy, zz ;

  Nx=info->N[0];
  Ny=info->N[1];
  Nz=info->N[2];

  /* X-dimension */
  if( info->bc.low[0]==Dirichlet || info->bc.upp[0]==Dirichlet )
  {
    for(k=0; k<info->N[2]; k++)
      for(j=0; j<info->N[1]; j++)
    {
      yy=info->meshlines[1][j]-info->beamwavg[1] ;
      zz=info->meshlines[2][k]-info->beamwavg[2] ;

      if( info->bc.low[0]==Dirichlet )
      {
        i=0 ;
        xx=info->meshlines[0][i]-info->beamwavg[0] ;
        ip=i+Nx*j+Nx*Ny*k;
        info->V[ip]=info->rho[ip] = info->Q * VFUNCTION(xx,yy,zz,info->beamwstd[0],info->beamwstd[1],info->beamwstd[2]) ;
      }

      if( info->bc.upp[0]==Dirichlet )
      {
        i=Nx-1 ;
        xx=info->meshlines[0][i]-info->beamwavg[0] ;
        ip=i+Nx*j+Nx*Ny*k;
        info->V[ip]=info->rho[ip] = info->Q * VFUNCTION(xx,yy,zz,info->beamwstd[0],info->beamwstd[1],info->beamwstd[2]) ;
      }
    }
  }

  /* Y-dimension */
  if( info->bc.low[1]==Dirichlet || info->bc.upp[1]==Dirichlet )
  {
    for(k=0; k<info->N[2]; k++)
      for(i=0; i<info->N[0]; i++)
    {
      xx=info->meshlines[0][i]-info->beamwavg[0] ;
      zz=info->meshlines[2][k]-info->beamwavg[2] ;

      if( info->bc.low[1]==Dirichlet )
      {
        j=0 ;
        yy=info->meshlines[1][j]-info->beamwavg[1] ;
        ip=i+Nx*j+Nx*Ny*k;
        info->V[ip]=info->rho[ip] = info->Q * VFUNCTION(xx,yy,zz,info->beamwstd[0],info->beamwstd[1],info->beamwstd[2]) ;
      }

      if( info->bc.upp[1]==Dirichlet )
      {
        j=Ny-1 ;
        yy=info->meshlines[1][j]-info->beamwavg[1] ;
        ip=i+Nx*j+Nx*Ny*k;
        info->V[ip]=info->rho[ip] = info->Q * VFUNCTION(xx,yy,zz,info->beamwstd[0],info->beamwstd[1],info->beamwstd[2]) ;
      }
    }
  }

  /* Z-dimension */
  if( info->bc.low[2]==Dirichlet || info->bc.upp[2]==Dirichlet )
  {
    for(j=0; j<info->N[1]; j++)
      for(i=0; i<info->N[0]; i++)
    {
      xx=info->meshlines[0][i]-info->beamwavg[0] ;
      yy=info->meshlines[1][j]-info->beamwavg[1] ;

      if( info->bc.low[2]==Dirichlet )
      {
        k=0 ;
        zz=info->meshlines[2][k]-info->beamwavg[2] ;
        ip=i+Nx*j+Nx*Ny*k;
        info->V[ip]=info->rho[ip] = info->Q * VFUNCTION(xx,yy,zz,info->beamwstd[0],info->beamwstd[1],info->beamwstd[2]) ;
      }

      if( info->bc.upp[2]==Dirichlet )
      {
        k=Nz-1 ;
        zz=info->meshlines[2][k]-info->beamwavg[2] ;
        ip=i+Nx*j+Nx*Ny*k;
        info->V[ip]=info->rho[ip] = info->Q * VFUNCTION(xx,yy,zz,info->beamwstd[0],info->beamwstd[1],info->beamwstd[2]) ;
      }
    }
  }
}


/* Print histogram of number of particles per mesh */
static void printmeshstats(struct sc3dmesh_info *info)
{
  int i ;
  unsigned int *n, nmax, ntotal ;
  unsigned int *nhis ;
  unsigned int j[3] ;

  /* Allocate memory for n array */
  n=(unsigned int *)gptmalloc(info->N[0]*info->N[1]*info->N[2]*sizeof(*n)) ;
  for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++) n[i]=0 ;

  /* Count number of particles in mesh */
  for(i=0 ; i<numpar; i++ )
    if( pars[i].alive )
      if( getindex(j,info->ro[i],info->meshlines, info->N)!=0 )
        n[(j[2]*info->N[1]+j[1])*info->N[0]+j[0]]++ ;

  /* Determine maximum and allocate nhis array */
  nmax=0 ;
  for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++)
    if( n[i]>nmax ) nmax=n[i] ;

  /* Allocate memory and calculate histogram */
  nhis=(unsigned int *)gptmalloc((nmax+1)*sizeof(*nhis)) ;
  for(i=0 ; i<=nmax ; i++) nhis[i]=0 ;
  for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++) nhis[n[i]]++ ;

  /* Print histogram */
  printf( "Meshnode usage\n" ) ;
  for(i=0 ; i<=nmax ; i++)
    printf( "%d: %d\n", i, nhis[i] ) ;

  /* Calculate total */
  ntotal=0 ;
  for(i=0 ; i<=nmax ; i++) ntotal += i*nhis[i] ;
  printf( "Total: %d\n", ntotal ) ;

  /* Free memory */
  gptfree(nhis) ;
  gptfree(n) ;
}

static double Vsphere(double x, double y, double z, double sx, double sy, double sz )
{
  double r=sqrt(x*x+y*y+z*z);

  return( 1/(4*gpt_pi*gpt_eps0*r) ) ;
}

static double Vbricktmp(double x, double y, double z)
{
  double r=sqrt(x*x+y*y+z*z) ;

  return(
    -y*z*log(x+r)+0.5*x*x*atan((y*z)/(x*r))
    -z*x*log(y+r)+0.5*y*y*atan((z*x)/(y*r))
    -x*y*log(z+r)+0.5*z*z*atan((x*y)/(z*r)) ) ;
}

static double Vbrick( double x, double y, double z, double sx, double sy, double sz )
{
  double tmp ;

  sx *= sqrt(3.0) ;
  sy *= sqrt(3.0) ;
  sz *= sqrt(3.0) ;

  tmp =
    +Vbricktmp(x-sx,y-sy,z-sz)
    -Vbricktmp(x-sx,y-sy,z+sz)
    -Vbricktmp(x-sx,y+sy,z-sz)
    +Vbricktmp(x-sx,y+sy,z+sz)
    -Vbricktmp(x+sx,y-sy,z-sz)
    +Vbricktmp(x+sx,y-sy,z+sz)
    +Vbricktmp(x+sx,y+sy,z-sz)
    -Vbricktmp(x+sx,y+sy,z+sz) ;
    
  return( tmp/(8*sx*sy*sz*4*gpt_pi*gpt_eps0) ) ;
}


/* Interpolate electric field in rest frame and transform back to lab frame */
static void interpolatetransform(struct sc3dmesh_info *info)
{
  int i ;
  struct parmeshinfo *pm = (struct parmeshinfo *)info->sabuff->memory ;

//#pragma omp parallel for if(GPTOMPSIMP(numpar))
  for(i=0 ; i<numpar; i++ ) if( pars[i].alive )
  {
    unsigned int *j, k ;
    double E[3], Erotated[3], cdot ;

    j = pm[i].jprev ;    
    if( j[0]!=-1 )
    {
      /* Interpolate V to E field in rest frame */
//      linearinterpolate3D(info->V,j,info->ro[i],info->meshlines,info->N, Erotated) ;
      quadraticinterpolate3D(info->V,j,info->ro[i],info->meshlines,info->N, Erotated) ;

      /* Unscale */
      if( info->bscaling )
      {
        for(k=0 ; k<3 ; k++) if( info->scale[ k     ]<1 ) Erotated[k]*=info->scale[ k     ] * info->scale[k] ;
        for(k=0 ; k<3 ; k++) if( info->scale[(k+1)%3]<1 ) Erotated[k]*=info->scale[(k+1)%3] ;
        for(k=0 ; k<3 ; k++) if( info->scale[(k+2)%3]<1 ) Erotated[k]*=info->scale[(k+2)%3] ;
      }

      /* Unrotate: Erotated->E */
      gptdirectiontoWCS(&info->rotation,Erotated,E) ;

      /* Transform to rest frame */
      cdot = gptVECINP(E,info->GBo) / (info->Go+1.0) ;
 
      pars[i].WE[0] += info->Go*E[0] - cdot*info->GBo[0] ;
      pars[i].WE[1] += info->Go*E[1] - cdot*info->GBo[1] ;
      pars[i].WE[2] += info->Go*E[2] - cdot*info->GBo[2] ;

      pars[i].WB[0] += (info->GBo[1]*E[2] - info->GBo[2]*E[1]) / gpt_c ;
      pars[i].WB[1] += (info->GBo[2]*E[0] - info->GBo[0]*E[2]) / gpt_c ;
      pars[i].WB[2] += (info->GBo[0]*E[1] - info->GBo[1]*E[0]) / gpt_c ;
    }
  }
}


/* Helper routines                                                                      */
/* ========================================================================================= */
static double gdfasqrt(double x)
{
  if( x<0 ) return(0) ;
  return( sqrt(x) ) ;
}


/* ACAD output routines                                                                      */
/* ========================================================================================= */

/* Write the mesh and all particles in rest frame to AutoCAD script */
static void writeacadmesh(char *filename,struct sc3dmesh_info *info)
{
  FILE *fp ;
  int i ;
  double **ml = info->meshlines ;

  if( (fp=fopen(filename,"w"))==NULL )
    gpterror( "%s: %s\n", filename, strerror(errno)) ;

  fprintf(fp,"osnap off\n" ) ;

  /* Particles */
  for(i=0 ; i<numpar ; i++) if(pars[i].alive )
    fprintf(fp,"point %f,%f,%f\n", info->ro[i][0], info->ro[i][1], info->ro[i][2] ) ;

  /* Mesh faces in x-direction */
  for(i=0 ; i<info->N[0] ; i++)
    fprintf( fp, "3dface %f,%f,%f %f,%f,%f %f,%f,%f %f,%f,%f \n",
      ml[0][i],ml[1][           0],ml[2][           0],
      ml[0][i],ml[1][info->N[1]-1],ml[2][           0],
      ml[0][i],ml[1][info->N[1]-1],ml[2][info->N[2]-1],
      ml[0][i],ml[1][           0],ml[2][info->N[2]-1] ) ;

  /* Mesh faces in y-direction */
  for(i=0 ; i<info->N[1] ; i++)
    fprintf( fp, "3dface %f,%f,%f %f,%f,%f %f,%f,%f %f,%f,%f \n",
      ml[0][           0],ml[1][i],ml[2][           0],
      ml[0][info->N[0]-1],ml[1][i],ml[2][           0],
      ml[0][info->N[0]-1],ml[1][i],ml[2][info->N[2]-1],
      ml[0][           0],ml[1][i],ml[2][info->N[2]-1] ) ;

  /* Mesh faces in z-direction */
  for(i=0 ; i<info->N[2] ; i++)
    fprintf( fp, "3dface %f,%f,%f %f,%f,%f %f,%f,%f %f,%f,%f \n",
      ml[0][           0],ml[1][           0],ml[2][i],
      ml[0][info->N[0]-1],ml[1][           0],ml[2][i],
      ml[0][info->N[0]-1],ml[1][info->N[1]-1],ml[2][i],
      ml[0][           0],ml[1][info->N[1]-1],ml[2][i] ) ;

  fprintf(fp,"zoom e\n") ;

  fclose(fp) ;
}





/* Write charge to a stack of surface plots as AutoCAD script */
static void writeacadcharge(char *filename,struct sc3dmesh_info *info)
{
  FILE *fp ;
  unsigned int i,j,k ;
  double rhomax, meshmin ;
  double **ml = info->meshlines ;
  double *rho = info->rho ;

  /* Calculate maximum charge density */
  rhomax=0 ;
  for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++)
    if( fabs(rho[i])>rhomax ) rhomax=fabs(rho[i]) ;

  /* Calculate minimum meshline distance in z-direction */
  meshmin=DBL_MAX ;
  for(i=0 ; i<info->N[2]-1 ; i++)
    if( ml[2][i+1]-ml[2][i]<meshmin ) meshmin=ml[2][i+1]-ml[2][i] ;


  /* Write file */
  if( (fp=fopen(filename,"w"))==NULL )
    gpterror( "%s: %s\n", filename, strerror(errno)) ;

  fprintf(fp,"osnap off\n" ) ;

  for(k=0 ; k<info->N[2] ; k++)
  {
    fprintf(fp, "3dmesh %d %d ", info->N[0], info->N[1] ) ;

    for(i=0 ; i<info->N[0] ; i++)
      for(j=0 ; j<info->N[1] ; j++)
        fprintf(fp, "%f,%f,%f ",
          ml[0][i],
          ml[1][j],
          ml[2][k]+fabs(info->rho[(k*info->N[1]+j)*info->N[0]+i])*2*meshmin/rhomax
        ) ;
  }

  fprintf(fp,"zoom e\n") ;

  fclose(fp) ;

  /* Debug: maximum rho and minimum meshsize */
  printf( "rhomax=%f, meshmin=%f\n", rhomax, meshmin ) ;
}


/* Write potential to a stack of surface plots as AutoCAD script */
static void writeacadpotential(char *filename,struct sc3dmesh_info *info)
{
  FILE *fp ;
  unsigned int i,j,k ;
  double Vmax, meshmin ;
  double **ml = info->meshlines ;
  double *V = info->V ;

  /* Calculate maximum charge density */
  Vmax=0 ;
  for(i=0 ; i<info->N[0]*info->N[1]*info->N[2] ; i++)
    if( fabs(V[i])>Vmax ) Vmax=fabs(V[i]) ;

  /* Calculate minimum meshline distance in z-direction */
  meshmin=DBL_MAX ;
  for(i=0 ; i<info->N[2]-1 ; i++)
    if( ml[2][i+1]-ml[2][i]<meshmin ) meshmin=ml[2][i+1]-ml[2][i] ;


  /* Write file */
  if( (fp=fopen(filename,"w"))==NULL )
    gpterror( "%s: %s\n", filename, strerror(errno)) ;

  fprintf(fp,"osnap off\n" ) ;

  for(k=0 ; k<info->N[2] ; k++)
  {
    fprintf(fp, "3dmesh %d %d ", info->N[0], info->N[1] ) ;

    for(i=0 ; i<info->N[0] ; i++)
      for(j=0 ; j<info->N[1] ; j++)
        fprintf(fp, "%f,%f,%f ",
          ml[0][i],
          ml[1][j],
          ml[2][k]+fabs(V[(k*info->N[1]+j)*info->N[0]+i])*2*meshmin/Vmax
        ) ;
  }

  fprintf(fp,"zoom e\n") ;

  fclose(fp) ;

  /* Debug: maximum rho and minimum meshsize */
  printf( "Vmax=%f, meshmin=%f\n", Vmax, meshmin ) ;
}


/* 3D mesh helper functions                                                                        */
/* =============================================================================================== */

/* Helper routine, initialize new parmesh buffer items to zero */
static void parmeshnotify(gptsyncarray *psa, gptsyncstatus stat, int n, struct sc3dmesh_info *info)
{
  struct parmeshinfo *pm ;
  int i ;

  switch(stat)
  {
    case gptcreatearray:
      pm = (struct parmeshinfo *)psa->memory ;
      for(i=0 ; i<psa->itemcount ; i++)
        pm[i].jprev[0] = pm[i].jprev[1] = pm[i].jprev[2] = -1 ;
      break ;

    case gptadditem:
      pm = (struct parmeshinfo *)psa->memory ;
      pm[n].jprev[0] = pm[n].jprev[1] = pm[n].jprev[2] = -1 ;
      break ;

    default: ;
  }
}

/* Helper routine: Get meshbox index from particle coordinate */
static int getindex(unsigned int *j, double *ro,double **meshlines,unsigned int *N)
{
  int i, jl, jh, jm ;
  int inc ;

  for(i=0 ; i<3 ; i++)
  {
    jl = j[i] ;

    if(jl<0 || jl>=(int)N[i])
    {
      jl=-1 ;
      jh=N[i] ;
    } else
    {
      inc=1 ;
      if( ro[i]>=meshlines[i][jl] )
      {
        if(jl==(int)N[i]-1) return(0) ;
        jh=jl+1 ;
        while(ro[i]>=meshlines[i][jh])
        {
          jl=jh ;
          inc *= 2 ;
          jh=jl+inc ;
          if(jh>(int)N[i]-1)
          {
            jh=N[i] ;
            break ;
          }
        }
      } else
      {
        if(jl==0) return(0) ;
        jh=jl-- ;
        while(ro[i]<meshlines[i][jl])
        {
          jh=jl ;
          inc *= 2 ;
          jl=jh-inc ;
          if(jl<0)
          {
            jl=0 ;
            break ;
          }
        }
      }
    }

    while( jh-jl>1 )
    {
      jm = (jh+jl)/2 ;
      if( ro[i] > meshlines[i][jm] )
        jl = jm ;
      else
        jh = jm ;
    }

    if( jl==((int)N[i]-1) || jl==-1 ) return(0) ;

    j[i] = jl ;
  }

  return(1) ;
}


/* Trilinear interpolation, derivative only */
static void linearinterpolate3D(double *V,unsigned int *j,double *r,double **meshlines,unsigned int *N, double *E)
{
  unsigned int N1,N2,N3,N4,N5,N6,N7,N8 ; 
  double dx,dy,dz ;
  double V1,V2,V3,V4,V5,V6,V7,V8 ;
  double ft,fu,fv, gt,gu,gv ;

  N1 = (j[2]*N[1]+j[1])*N[0]+j[0] ;
  N2 = N1+1 ;
  N3 = N1+N[0]+1 ;
  N4 = N1+N[0] ;
  N5 = N1+N[0]*N[1] ;
  N6 = N1+N[0]*N[1]+1 ;
  N7 = N1+N[0]*(N[1]+1)+1 ;
  N8 = N1+N[0]*(N[1]+1) ;

  dx = (meshlines[0][j[0]+1] - meshlines[0][j[0]]) ;
  dy = (meshlines[1][j[1]+1] - meshlines[1][j[1]]) ;
  dz = (meshlines[2][j[2]+1] - meshlines[2][j[2]]) ;

  ft= (r[0]-meshlines[0][j[0]]) / dx ;
  fu= (r[1]-meshlines[1][j[1]]) / dy ;
  fv= (r[2]-meshlines[2][j[2]]) / dz ;
  gt = 1-ft ;
  gu = 1-fu ;
  gv = 1-fv ;

  V1 = V[N1] ; V2 = V[N2] ; V3 = V[N3] ; V4 = V[N4] ;
  V5 = V[N5] ; V6 = V[N6] ; V7 = V[N7] ; V8 = V[N8] ;

  E[0] = (gv*(gu*(V1-V2) - fu*(V3-V4)) + fv*(gu*(V5-V6) - fu*(V7-V8))) / dx ;
  E[1] = (gv*(gt*(V1-V4) + ft*(V2-V3)) + fv*(gt*(V5-V8) + ft*(V6-V7))) / dy ;
  E[2] = (gt*(gu*(V1-V5) + fu*(V4-V8)) + ft*(gu*(V2-V6) + fu*(V3-V7))) / dz ;
}

static void quadraticinterpolate3D(double *V,unsigned int *j,double *r,double **meshlines,unsigned int *N, double *E)
{
  double ft,fu,fv, gt,gu,gv ;
  int N000 ;
  double Ex0,Ex1, Ey0, Ey1, Ez0, Ez1 ;

  /* Test if particle is within a mesh-box near the boundary */
  if( j[0]==0 || j[1]==0 || j[2]==0 || j[0]==N[0]-2 || j[1]==N[1]-2 || j[2]==N[2]-2 )
  {
    linearinterpolate3D(V,j,r,meshlines,N,E) ;
    return ;
  }

  ft= (r[0]-meshlines[0][j[0]]) / (meshlines[0][j[0]+1] - meshlines[0][j[0]]) ;
  fu= (r[1]-meshlines[1][j[1]]) / (meshlines[1][j[1]+1] - meshlines[1][j[1]]) ;
  fv= (r[2]-meshlines[2][j[2]]) / (meshlines[2][j[2]+1] - meshlines[2][j[2]]) ;
  gt = 1-ft ;
  gu = 1-fu ;
  gv = 1-fv ;

  N000 = (j[2]*N[1]+j[1])*N[0]+j[0] ;

#define VV(j0,j1,j2) V[N000+((j2)*N[1]+(j1))*N[0]+(j0)]
  Ex0 =-((gu*gv*VV(+1,0,0)+gu*fv*VV(+1,0,1)+fu*gv*VV(+1,1,0)+fu*fv*VV(+1,1,1)) -
         (gu*gv*VV(-1,0,0)+gu*fv*VV(-1,0,1)+fu*gv*VV(-1,1,0)+fu*fv*VV(-1,1,1)) ) / (meshlines[0][j[0]+1]-meshlines[0][j[0]-1]) ;
  Ex1 =-((gu*gv*VV(+2,0,0)+gu*fv*VV(+2,0,1)+fu*gv*VV(+2,1,0)+fu*fv*VV(+2,1,1)) -
         (gu*gv*VV( 0,0,0)+gu*fv*VV( 0,0,1)+fu*gv*VV( 0,1,0)+fu*fv*VV( 0,1,1)) ) / (meshlines[0][j[0]+2]-meshlines[0][j[0]  ]) ;

  Ey0 =-((gt*gv*VV(0,+1,0)+gt*fv*VV(0,+1,1)+ft*gv*VV(1,+1,0)+ft*fv*VV(1,+1,1)) -
         (gt*gv*VV(0,-1,0)+gt*fv*VV(0,-1,1)+ft*gv*VV(1,-1,0)+ft*fv*VV(1,-1,1)) ) / (meshlines[1][j[1]+1]-meshlines[1][j[1]-1]) ;
  Ey1 =-((gt*gv*VV(0,+2,0)+gt*fv*VV(0,+2,1)+ft*gv*VV(1,+2,0)+ft*fv*VV(1,+2,1)) -
         (gt*gv*VV(0, 0,0)+gt*fv*VV(0, 0,1)+ft*gv*VV(1, 0,0)+ft*fv*VV(1, 0,1)) ) / (meshlines[1][j[1]+2]-meshlines[1][j[1]  ]) ;

  Ez0 =-((gt*gu*VV(0,0,+1)+gt*fu*VV(0,1,+1)+ft*gu*VV(1,0,+1)+ft*fu*VV(1,1,+1)) -
         (gt*gu*VV(0,0,-1)+gt*fu*VV(0,1,-1)+ft*gu*VV(1,0,-1)+ft*fu*VV(1,1,-1)) ) / (meshlines[2][j[2]+1]-meshlines[2][j[2]-1]) ;
  Ez1 =-((gt*gu*VV(0,0,+2)+gt*fu*VV(0,1,+2)+ft*gu*VV(1,0,+2)+ft*fu*VV(1,1,+2)) -
         (gt*gu*VV(0,0, 0)+gt*fu*VV(0,1, 0)+ft*gu*VV(1,0, 0)+ft*fu*VV(1,1, 0)) ) / (meshlines[2][j[2]+2]-meshlines[2][j[2]  ]) ;
#undef VV
 
  E[0] = gt*Ex0 + ft*Ex1 ;
  E[1] = gu*Ey0 + fu*Ey1 ;
  E[2] = gv*Ez0 + fv*Ez1 ;
}

/* Tricubic inetpolations, derivative only */
static void cubicinterpolate3D(double *V,unsigned int *j,double *r,double **meshlines,unsigned int *N, double *E)
{
/*
  ??
*/
}
