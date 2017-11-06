/* map2D_E.c - 2D electrostatic field map in the zx-plane */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <algorithm>
using namespace std ;

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "elem.h"


struct map2d_info
{
  /* GDF information */
  struct gdfmem gm ;

  /* Thickness/2 */
  double L ;

  /* Coordinate information, only used for sorting */
  double *x, *z ;

  /* Array size information */
  double zmin, zmax, zdelta ;
  double xmin, xmax, xdelta ;
  int numx, numz ;

  /* Field information */
  double *fx, *fz;
} ;


static int map2d_sim(gptpar *par, double t, struct map2d_info *info) ;
static void map2d_exit(struct map2d_info *info) ;


/* Helper functions for sorting */
static void swapmap2D(struct map2d_info *info, int i, int j) ;
static void qsortmap2D(struct map2d_info *info, int left, int right) ;


void map2Dr_E_init(gptinit *init)
{
  struct map2d_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  int numc ;
  char *fxname, *fzname ;
  char *xname,*zname ;
  double ffac ;

  /* Grid characteristics */
  int i,j,points,N ;

  /* GDF information */
  struct gdff ingdff ;
  char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=7 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,x,z,Ex,Ez,L,Efac)\n", gptgetname(init) ) ;

  info = (struct map2d_info *)gptmalloc( sizeof(struct map2d_info) ) ;

  /* Initialize GDF */
  filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  /* Retrieve x and z arrays */
  xname = gptgetargstring(init,2) ;
  zname = gptgetargstring(init,3) ;

  ds=gptinputgetgroup(&info->gm,filename,xname) ;
  info->x=gptinputdoublearray(ds,xname,&points) ;
  info->z=gptinputdoublearraypoints(ds,zname,points) ;

  /* Retrieve field from datafile */
  fxname   = gptgetargstring(init,4) ;
  fzname   = gptgetargstring(init,5) ;

  info->fx = gptinputdoublearraypoints(ds,fxname  ,points) ;
  info->fz = gptinputdoublearraypoints(ds,fzname  ,points) ;

  /* Get thickness and Efac */
  info->L = gptgetargdouble(init,6)/2 ;
  ffac    = gptgetargdouble(init,7) ;
  for(i=0 ; i<points ; i++) info->fx[i]   *= ffac ;
  for(i=0 ; i<points ; i++) info->fz[i]   *= ffac ;

  /* Sort all datapoints in correct order: X runs fastest */
  for(i=0 ; i<points ; i++) swapmap2D(info,rand()%points,rand()%points) ;
  qsortmap2D(info,0,points-1) ;

  /* Obtain file characteristics */
  info->xmin = info->x[0] ; 
  info->zmin = info->z[0] ; 
  info->xmax = info->x[points-1] ; 
  info->zmax = info->z[points-1] ;
  i=1 ;
  while(i<points && info->x[i-1]<info->x[i]) i++ ;
  info->numx = i ; 
  info->numz = points/info->numx ; 
  if( info->numx<2 ) gpterror( "%s: All X-coordinates are equal to %g\n", filename, info->xmin ) ;
  if( info->numz<2 ) gpterror( "%s: All Z-coordinates are equal to %g\n", filename, info->zmin ) ;
  info->xdelta = (info->xmax-info->xmin)/(info->numx-1) ; 
  info->zdelta = (info->zmax-info->zmin)/(info->numz-1) ;

  /* Test file for correctness */
  for(j=0 ; j<info->numz ; j++) for(i=0 ; i<info->numx ; i++)
  {
    N=j*info->numx+i ;
    if( info->x[N] < info->xmin+i*info->xdelta - info->xdelta/8 ||
        info->x[N] > info->xmin+i*info->xdelta + info->xdelta/8 ||
        info->z[N] < info->zmin+j*info->zdelta - info->zdelta/8 ||
        info->z[N] > info->zmin+j*info->zdelta + info->zdelta/8 )
    {
      gptwarning( "x-direction: %d samples in range [%g,%g] with delta %g\n", info->numx, info->xmin, info->xmax, info->xdelta ) ;
      gptwarning( "z-direction: %d samples in range [%g,%g] with delta %g\n", info->numz, info->zmin, info->zmax, info->zdelta ) ;
      gpterror( "%s: No rectangular grid with the above parameters at (%g,%g)\n", filename, info->x[N], info->z[N] ) ;
    }  
  }

  gptaddEBelement( init, map2d_sim, map2d_exit, GPTELEM_GLOBAL, info ) ;
}


static int map2d_sim(gptpar *par, double t, struct map2d_info *info)
{
  int i,j,N ;
  double ft,fu, gt, gu ;
  int N1,N2,N3,N4 ;
  double f1,f2,f3,f4 ;

  if( Z<info->zmin || Z>info->zmax ) return( 0 ) ;
  if( X<info->xmin || X>info->xmax ) return( 0 ) ;
  if( fabs(Y)>info->L ) return( 0 ) ;

  /* Calculate master index and offset fractions */
  i = (int)((X-info->xmin)/info->xdelta) ;
  j = (int)((Z-info->zmin)/info->zdelta) ;
  N = j*info->numx + i ;
  ft= (X-info->xmin-i*info->xdelta)/info->xdelta ;
  fu= (Z-info->zmin-j*info->zdelta)/info->zdelta ;
  gt = 1-ft ;
  gu = 1-fu ;

  /* Calculate boundary offsets */
  N1 = N ;
  N2 = N+1 ;
  N3 = N+info->numx+1 ;
  N4 = N+info->numx ;

/* Consistency check in case of problems with this element
 * if(i<0 || i>=info->numx-1) gpterror( "Error in x range: %lf maps to %d\n", X, i ) ;
 * if(j<0 || j>=info->numz-1) gpterror( "Error in z range: %lf maps to %d\n", Z, j ) ;
 * if(ft<0 || ft>1+16*DBL_EPSILON ) gpterror( "x fraction is %.16f at x=%f\n", ft, X ) ;
 * if(fu<0 || fu>1+16*DBL_EPSILON ) gpterror( "z fraction is %.16f at z=%f\n", fu, Z ) ;
 */

  f1=gt*gu ; f2=ft*gu ; f3=ft*fu ; f4=gt*fu ;
  EX   = f1*info->fx[N1]   + f2*info->fx[N2]   + f3*info->fx[N3]   + f4*info->fx[N4] ;
  EZ   = f1*info->fz[N1]   + f2*info->fz[N2]   + f3*info->fz[N3]   + f4*info->fz[N4] ;

  return( 1 ) ;
}

static void map2d_exit(struct map2d_info *info)
{
  gdfmfreechilds( &info->gm.ds ) ;
  free(info) ;
}


/* swap two points, keeping all arrays in sync */
static void swapmap2D(struct map2d_info *info, int i, int j)
{
  if( i!=j &&
      info->x[i]==info->x[j] &&
      info->z[i]==info->z[j] )
    gpterror( "Duplicate points at (%g,%g)\n", info->x[i], info->z[i] ) ;

  swap(info->x[i],info->x[j]) ;
  swap(info->z[i],info->z[j]) ;

  swap(info->fx[i],info->fx[j]  ) ;
  swap(info->fz[i],info->fz[j]  ) ;
}


/* Sort all datapoints in correct order: R runs fastest */
static void qsortmap2D(struct map2d_info *info, int left, int right)
{
  int i, last ;

  if( left>=right ) return ;

  swapmap2D(info,left,(left+right)/2) ;
  last = left ;
  for(i=left+1 ; i<=right ; i++)
    if(info->z[i]<info->z[left] ||
      (info->x[i]<info->x[left] && info->z[i]==info->z[left]))
        swapmap2D(info,++last,i) ;
  swapmap2D(info,left,last) ;
  qsortmap2D(info,left,last-1) ;
  qsortmap2D(info,last+1,right) ;
}
