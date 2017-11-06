/* map2DTM.c - rectangular, cylindrical symmetric field map 
 *              for standing and traveling wave cavity
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

struct mapcav_info
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

  /* Cavity information */
  double k, w, phi;
} ;


static int mapcav_sim(gptpar *par, double t, struct mapcav_info *info) ;
static void mapcav_exit(struct mapcav_info *info) ;


/* Helper functions for sorting */
static void swapmapcav(struct mapcav_info *info, int i, int j) ;
static void qsortmapcav(struct mapcav_info *info, int left, int right) ;


void map25D_TM_init(gptinit *init)
{
  struct mapcav_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  int numc ;
  char *frname, *fzname, *fphiname ;
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
  if( numc!=10 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,r,z,Er,Ez,Bphi,ffac,k,phi,w)\n", gptgetname(init) ) ;

  info = (struct mapcav_info *)gptmalloc( sizeof(struct mapcav_info) ) ;

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
  fphiname = gptgetargstring(init,6) ;

  info->fr   = gptinputdoublearraypoints(ds,frname  ,points) ;
  info->fz   = gptinputdoublearraypoints(ds,fzname  ,points) ;
  info->fphi = gptinputdoublearraypoints(ds,fphiname,points) ;

  ffac = gptgetargdouble(init,7) ;
  for(i=0 ; i<points ; i++) info->fr[i]   *= ffac ;
  for(i=0 ; i<points ; i++) info->fz[i]   *= ffac ;
  for(i=0 ; i<points ; i++) info->fphi[i] *= ffac ;

  info->k   = gptgetargdouble(init,8) ;
  info->phi = gptgetargdouble(init,9) ;
  info->w   = gptgetargdouble(init,10) ;

  /* Sort all datapoints in correct order: R runs fastest */
  for(i=0 ; i<points ; i++) swapmapcav(info,rand()%points,rand()%points) ;
  qsortmapcav(info,0,points-1) ;

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

  gptaddEBelement( init, mapcav_sim, mapcav_exit, GPTELEM_GLOBAL, info ) ;
}


static int mapcav_sim(gptpar *par, double t, struct mapcav_info *info)
{
  double R ;
  int i,j,N ;
  double ft,fu, gt, gu ;
  int N1,N2,N3,N4 ;
  double f1,f2,f3,f4 ;
  double Er, Bphi ;

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
  Bphi = f1*info->fphi[N1] + f2*info->fphi[N2] + f3*info->fphi[N3] + f4*info->fphi[N4] ;

  Er   *=  cos(info->w*t - info->k*Z + info->phi) ;
  EZ   *=  cos(info->w*t - info->k*Z + info->phi) ;
  Bphi *= -sin(info->w*t - info->k*Z + info->phi) ;

  gptr2carth(Er,X,Y,&EX,&EY) ;
  gptrphi2carth(0.0,Bphi,X,Y,&BX,&BY) ;

  return( 1 ) ;
}


static void mapcav_exit(struct mapcav_info *info)
{
  gdfmfreechilds( &info->gm.ds ) ;
  free(info) ;
}


/* swap two points, keeping all arrays in sync */
static void swapmapcav(struct mapcav_info *info, int i, int j)
{
  if( i!=j &&
      info->r[i]==info->r[j] &&
      info->z[i]==info->z[j] )
    gpterror( "Duplicate points at (%g,%g)\n", info->r[i], info->z[i] ) ;

  swap(info->r[i],info->r[j]) ;
  swap(info->z[i],info->z[j]) ;

  swap(info->fr[i]  ,info->fr[j]  ) ;
  swap(info->fz[i]  ,info->fz[j]  ) ;
  swap(info->fphi[i],info->fphi[j]) ;
}


/* Sort all datapoints in correct order: R runs fastest */
static void qsortmapcav(struct mapcav_info *info, int left, int right)
{
  int i, last ;

  if( left>=right ) return ;

  swapmapcav(info,left,(left+right)/2) ;
  last = left ;
  for(i=left+1 ; i<=right ; i++)
    if(info->z[i]<info->z[left] ||
      (info->r[i]<info->r[left] && info->z[i]==info->z[left]))
        swapmapcav(info,++last,i) ;
  swapmapcav(info,left,last) ;
  qsortmapcav(info,left,last-1) ;
  qsortmapcav(info,last+1,right) ;
}
