/* map2D_E.c - 2D rectangular, cylindrical symmetric field map 
 *             for the electric field
 */

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

  /* Coordinate information, only used for sorting */
  double *r, *z ;

  /* Array size information */
  double zmin, zmax, zdelta ;
  double rmin, rmax, rdelta ;
  int numr, numz ;

  /* Field information */
  double *fr, *fphi, *fz;
} ;


static int map2d_sim(gptpar *par, double t, struct map2d_info *info) ;
static void map2d_exit(struct map2d_info *info) ;


/* Helper functions for sorting */
static void swapmap2D(struct map2d_info *info, int i, int j) ;
static void qsortmap2D(struct map2d_info *info, int left, int right) ;


void map2D_E_init(gptinit *init)
{
  struct map2d_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  int numc ;
  char *frname, *fzname ;
  char *rname,*zname ;
  double ffac ;

  /* Grid characteristics */
  int i,j,points,N ;

  /* GDF information */
  struct gdff ingdff ;
  char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=6 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,r,z,Er,Ez,Efac)\n", gptgetname(init) ) ;

  info = (struct map2d_info *)gptmalloc( sizeof(struct map2d_info) ) ;

  /* Initialize GDF */
  filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  /* Retrieve r and z arrays */
  rname = gptgetargstring(init,2) ;
  zname = gptgetargstring(init,3) ;

  ds=gptinputgetgroup(&info->gm,filename,rname) ;
  info->r=gptinputdoublearray(ds,rname,&points) ;
  info->z=gptinputdoublearraypoints(ds,zname,points) ;

  /* Retrieve field from datafile */
  frname   = gptgetargstring(init,4) ;
  fzname   = gptgetargstring(init,5) ;

  info->fr = gptinputdoublearraypoints(ds,frname  ,points) ;
  info->fz = gptinputdoublearraypoints(ds,fzname  ,points) ;

  ffac = gptgetargdouble(init,6) ;
  for(i=0 ; i<points ; i++) info->fr[i]   *= ffac ;
  for(i=0 ; i<points ; i++) info->fz[i]   *= ffac ;


  /* Sort all datapoints in correct order: R runs fastest */
  for(i=0 ; i<points ; i++) swapmap2D(info,rand()%points,rand()%points) ;
  qsortmap2D(info,0,points-1) ;

  /* Obtain file characteristics */
  info->rmin = info->r[0] ; 
  info->zmin = info->z[0] ; 
  info->rmax = info->r[points-1] ; 
  info->zmax = info->z[points-1] ;
  i=1 ;
  while(i<points && info->r[i-1]<info->r[i]) i++ ;
  info->numr = i ; 
  info->numz = points/info->numr ; 
  if( info->numr<2 ) gpterror( "%s: All R-coordinates are equal to %g\n", filename, info->rmin ) ;
  if( info->numz<2 ) gpterror( "%s: All Z-coordinates are equal to %g\n", filename, info->zmin ) ;
  info->rdelta = (info->rmax-info->rmin)/(info->numr-1) ; 
  info->zdelta = (info->zmax-info->zmin)/(info->numz-1) ;

  /* Test file for correctness */
  for(j=0 ; j<info->numz ; j++) for(i=0 ; i<info->numr ; i++)
  {
    N=j*info->numr+i ;
    if( info->r[N] < info->rmin+i*info->rdelta - info->rdelta/8 ||
        info->r[N] > info->rmin+i*info->rdelta + info->rdelta/8 ||
        info->z[N] < info->zmin+j*info->zdelta - info->zdelta/8 ||
        info->z[N] > info->zmin+j*info->zdelta + info->zdelta/8 )
    {
      gptwarning( "r-direction: %d samples in range [%g,%g] with delta %g\n", info->numr, info->rmin, info->rmax, info->rdelta ) ;
      gptwarning( "z-direction: %d samples in range [%g,%g] with delta %g\n", info->numz, info->zmin, info->zmax, info->zdelta ) ;
      gpterror( "%s: No rectangular grid with the above parameters at (%g,%g)\n", filename, info->r[N], info->z[N] ) ;
    }  
  }
  for(i=0 ; i<info->numr ; i++ ) if( info->r[i]<0.0 )
    gpterror( "Negative r-coordinates present in \"%s\"\n", filename ) ;

  gptaddEBelement( init, map2d_sim, map2d_exit, GPTELEM_GLOBAL, info ) ;
}


static int map2d_sim(gptpar *par, double t, struct map2d_info *info)
{
  double R ;
  int i,j,N ;
  double ft,fu, gt, gu ;
  int N1,N2,N3,N4 ;
  double f1,f2,f3,f4 ;
  double Er ;

  if( Z<info->zmin || Z>info->zmax ) return( 0 ) ;
  R = sqrt(X*X + Y*Y) ;
  if( R<info->rmin || R>info->rmax ) return( 0 ) ;

  /* Calculate master index and offset fractions */
  i = (int)((R-info->rmin)/info->rdelta) ;
  j = (int)((Z-info->zmin)/info->zdelta) ;
  N = j*info->numr + i ;
  ft= (R-info->rmin-i*info->rdelta)/info->rdelta ;
  fu= (Z-info->zmin-j*info->zdelta)/info->zdelta ;
  gt = 1-ft ;
  gu = 1-fu ;

  /* Calculate boundary offsets */
  N1 = N ;
  N2 = N+1 ;
  N3 = N+info->numr+1 ;
  N4 = N+info->numr ;

/* Consistency check in case of problems with this element
 * if(i<0 || i>=info->numr-1) gpterror( "Error in r range: %lf maps to %d\n", R, i ) ;
 * if(j<0 || j>=info->numz-1) gpterror( "Error in z range: %lf maps to %d\n", Z, j ) ;
 * if(ft<0 || ft>1+16*DBL_EPSILON ) gpterror( "r fraction is %.16f at r=%f\n", ft, R ) ;
 * if(fu<0 || fu>1+16*DBL_EPSILON ) gpterror( "z fraction is %.16f at z=%f\n", fu, Z ) ;
 */

  f1=gt*gu ; f2=ft*gu ; f3=ft*fu ; f4=gt*fu ;
  Er   = f1*info->fr[N1]   + f2*info->fr[N2]   + f3*info->fr[N3]   + f4*info->fr[N4] ;
  EZ   = f1*info->fz[N1]   + f2*info->fz[N2]   + f3*info->fz[N3]   + f4*info->fz[N4] ;

  gptr2carth(Er,X,Y,&EX,&EY) ;

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
      info->r[i]==info->r[j] &&
      info->z[i]==info->z[j] )
    gpterror( "Duplicate points at (%g,%g)\n", info->r[i], info->z[i] ) ;

  swap(info->r[i],info->r[j]) ;
  swap(info->z[i],info->z[j]) ;

  swap(info->fr[i]  ,info->fr[j]  ) ;
  swap(info->fz[i]  ,info->fz[j]  ) ;
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
      (info->r[i]<info->r[left] && info->z[i]==info->z[left]))
        swapmap2D(info,++last,i) ;
  swapmap2D(info,left,last) ;
  qsortmap2D(info,left,last-1) ;
  qsortmap2D(info,last+1,right) ;
}
