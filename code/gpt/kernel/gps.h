/* gps.h - Main header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

/* Program version */
#define PROGNAME "GPT"
#define PROGMAJ 3
#define PROGMIN 00

/* Useful macro's */
#define gptSQR(x) ((x)*(x))
#define gptVECSQR(x) (gptSQR((x)[0])+gptSQR((x)[1])+gptSQR((x)[2]))
#define gptVECLEN(x) sqrt(gptVECSQR(x))
#define gptVECINP(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

/* Constants */
#define gpt_c  299792458.0
#define gpt_e  2.71828182845904523536
#define gpt_pi 3.14159265358979323846
#define gpt_deg (180.0/gpt_pi)

#define gpt_ma  1.66057e-27
#define gpt_me  9.10953e-31
#define gpt_mp  1.67265e-27
#define gpt_qe -1.6021892e-19

#define gpt_mu0  (1e-7*4.0*gpt_pi)
#define gpt_eps0 (1/(gpt_c*gpt_c*gpt_mu0))


typedef struct gpttransform
{
  double o[3] ;
  double m[3][3] ;
  int mid ;            /* m is Identity matrix */
} gpttransform ;

struct EBF
{
  double E[3] ;
  double B[3] ;
  double F[3] ;
} ;


typedef struct par
{
/* Main particle coordinates */
  double *Wr  ;
  int offWr ;  /* Position returned by ODE driver */
  double *GBr ;
  int offGBr ; /* Position returned by ODE driver */

  double m ;
  double q ;
  double n ;
  double r2 ;

  double tstart ; /* Particle start time */

/* Used by elems */
  double r[3] ;
  double E[3]  ;
  double B[3]  ;

/* Derived variables */
  double G ;

/* Final fields in "wcs" */
  double WE[3] ;
  double WB[3] ;
  double WF[3];

/* Coordinate transform support */
  struct axis *paxis ;
  struct axis *newaxis ;

/* Local element support */
  int kl, kc, kr ;	/* Local elements are between kl<=k<kc en tussen kc<=k<kr */
  int el, er ;		/* 'hit' in ALL elements to the left and to the rigth respectively */

/* Miscellaneous */
  int persistent ; /* Particle can not be removed from system, previous alive functionality */
  int alive ;      /* Particle is present in this (sub)step, switched on/off as function of t>tstart */
  int tokill ;
  int ID ;
  struct gptparset *set ;

} gptpar ;


typedef struct gpttrajectory
{
  /* Start and endpoints in ECS */
  double rstart[3], rend[3] ;

  /* Particle direction, normalized direction and lambda: P=lambda*dr */
  double dr[3], ndr[3] ;
  double lambda ;

  /* Intersions point, surface normal and interpolated gamma */
  double P[3], n[3] ;
  double inp, GBint[3], Gint ;

  /* Boundary element */
  struct gptboundaryelem *pelem ;
} gpttrajectory ;


typedef int (*simfn)(gptpar *par,double t,void *info) ;
typedef void (*exitfn)(void *info) ;
typedef void (*gptscatterfnc)(gptpar *par,double t,double dt,struct gpttrajectory *trajectory, void *info) ;
typedef int  (*gptboundaryfnc)(gptpar *par, double t, double dt, gpttrajectory *trajectory, void *info ) ;

struct gptEBelem
{
  gpttransform e ;
  void *info ;
  simfn sim ;
  exitfn exit ;
} ;

struct gptboundaryelem
{
  /* Location */
  gpttransform e ;

  /* Boundary */
  void *info ;
  gptboundaryfnc bound ;

  /* Material */
  gptscatterfnc scatter ;
  void *scatterinfo ;
} ;

struct axis
{
  gpttransform a ;
  char *name ;
  int numlocal ;
  struct gptEBelem *local ;
  int numglobal ;
  struct gptEBelem *global ;
  int numboundary ;
  struct gptboundaryelem *boundary ;
  struct axis *next ;
} ;

struct arg
{
  double d ;
  char *s ;
} ;

typedef struct inputlist
{
  /* Arguments */
  int argc ;
  struct arg *args ;

  /* Coordinate transforms */
  struct axis *paxis ;
  gpttransform e ;

  /* Scattering */
  char *scattername ;
} gptinit ;

struct elementnames
{
  char *name ;
  void (*init)(struct inputlist *p);
} ;

/* Grayscale bitmap support (elem.c) */
struct grayscalebmp
{
  int xpels ;	// Number of horizontal pixels
  int ypels ;	// Number of vertical pixels
  double xres ;	// Horizontal pixels per meter
  double yres ;	// Vertical pixels per meter
  
  unsigned char *pixels ;
} ;
void readgrayscalebmp(struct grayscalebmp *pbmp,char *filename) ;

/* 3D reltivistic point-to-point space-charge (sc.c in gptelems) */
void sc3Dp2p( double t, double *x, double *p, double zmin, double zmax, int nOff, int nCPU ) ;

/* Main */
extern struct elementnames elems[] ;

extern int numpar ;
extern int numalivepar ;
extern int panicing ;
extern struct par *pars ;
extern int verbose ;
extern char *progname ;
extern char *gdfoutfo, *gdfoutfa ;

void terminate(char *fmt,...);

/* Error handling */
#ifndef _ERRMSG_DEFINED
#include "../../utils/utils.h"
#endif
extern struct errmsg gpsem[] ;
extern void (*agpserr)( char *fmt, ... ) ;
extern void (*agpswar)( char *fmt, ... ) ;
extern void (*aclnerr)( char *fmt, ... ) ;
extern void (*averb)( char *fmt, ... ) ;
#define gpterror (**agpserr)
#define gptwarning (**agpswar)
#define gptverbose (**averb)

#include "gpserr.h"
