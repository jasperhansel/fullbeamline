/* map3D_V.c - 3D rectangular field map for the electric potential */

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


struct map3d_info
{
  /* GDF information */
  struct gdfmem gm ;

  /* Coordinate information, only used for sorting */
  double *x, *y, *z ;

  /* Array size information */
  double xmin, xmax, xdelta ;
  double ymin, ymax, ydelta ;
  double zmin, zmax, zdelta ;
  int numx, numy, numz ;

  /* Field information */
  double *V ;
  char Vtype ;
} ;


static int map3d_sim(gptpar *par, double t, struct map3d_info *info) ;
static void map3d_exit(struct map3d_info *info) ;


/* Helper functions for sorting */
static void swapmap3D(struct map3d_info *info, int i, int j) ;
static void qsortmap3D(struct map3d_info *info, int left, int right) ;


void map3D_V_init(gptinit *init)
{
  struct map3d_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  int numc ;
  char *Vname ;
  char *xname,*yname,*zname ;
  double Vfac ;

  /* Grid characteristics */
  int i,j,k,points,N ;

  /* GDF information */
  struct gdff ingdff ;
  char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  numc=gptgetargnum(init) ;
  if( numc!=6 && numc!=7 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,x,y,z,V,Vfac,[Vtype])\n", gptgetname(init) ) ;

  info = (struct map3d_info *)gptmalloc( sizeof(struct map3d_info) ) ;

  /* Initialize GDF */
  filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  /* Retrieve x, y and z arrays */
  xname = gptgetargstring(init,2) ;
  yname = gptgetargstring(init,3) ;
  zname = gptgetargstring(init,4) ;

  ds=gptinputgetgroup(&info->gm,filename,xname) ;
  info->x=gptinputdoublearray(ds,xname,&points) ;
  info->y=gptinputdoublearraypoints(ds,yname,points) ;
  info->z=gptinputdoublearraypoints(ds,zname,points) ;

  /* Retrieve potential from datafile */
  Vname = gptgetargstring(init,5) ;
  info->V = gptinputdoublearraypoints(ds,Vname,points) ;

  Vfac = gptgetargdouble(init,6) ;
  for(i=0 ; i<points ; i++) info->V[i] *= Vfac ;

  info->Vtype = 'l' ;
  if(numc==7)
  {
    info->Vtype = tolower(gptgetargstring(init,7)[0]) ;
    if(info->Vtype!='l') gpterror("Unknown type for V interpolation: %c\n", info->Vtype ) ;
  }


  /* Sort all datapoints in correct order: x runs fastest */
  for(i=0 ; i<points ; i++) swapmap3D(info,rand()%points,rand()%points) ;
  qsortmap3D(info,0,points-1) ;

  /* Obtain file characteristics */
  info->xmin = info->x[0] ; 
  info->ymin = info->y[0] ; 
  info->zmin = info->z[0] ; 
  info->xmax = info->x[points-1] ; 
  info->ymax = info->y[points-1] ; 
  info->zmax = info->z[points-1] ;
  i=1 ;
  while(i<points && info->x[i-1]<info->x[i]) i++ ;
  info->numx = i ;
  while(i<points && info->y[i-info->numx]<info->y[i]) i+=info->numx ;
  info->numy = i/info->numx ;
  info->numz = points/info->numx/info->numy ;
  if( info->numx<2 ) gpterror( "%s: All X-coordinates are equal to %g\n", filename, info->xmin ) ;
  if( info->numy<2 ) gpterror( "%s: All Y-coordinates are equal to %g\n", filename, info->ymin ) ;
  if( info->numz<2 ) gpterror( "%s: All Z-coordinates are equal to %g\n", filename, info->zmin ) ;
  info->xdelta = (info->xmax-info->xmin)/(info->numx-1) ; 
  info->ydelta = (info->ymax-info->ymin)/(info->numy-1) ; 
  info->zdelta = (info->zmax-info->zmin)/(info->numz-1) ;

  /* Test file for correctness */
  for(k=0 ; k<info->numz ; k++) for(j=0 ; j<info->numy ; j++) for(i=0 ; i<info->numx ; i++)
  { 
    N=(k*info->numy+j)*info->numx+i ;
    if( info->x[N] < info->xmin+i*info->xdelta - info->xdelta/8 ||
        info->x[N] > info->xmin+i*info->xdelta + info->xdelta/8 ||
        info->y[N] < info->ymin+j*info->ydelta - info->ydelta/8 ||
        info->y[N] > info->ymin+j*info->ydelta + info->ydelta/8 ||
        info->z[N] < info->zmin+k*info->zdelta - info->zdelta/8 ||
        info->z[N] > info->zmin+k*info->zdelta + info->zdelta/8 )
    {
      gptwarning( "x-direction: %d samples in range [%g,%g] with delta %g\n", info->numx, info->xmin, info->xmax, info->xdelta ) ;
      gptwarning( "y-direction: %d samples in range [%g,%g] with delta %g\n", info->numy, info->ymin, info->ymax, info->ydelta ) ;
      gptwarning( "z-direction: %d samples in range [%g,%g] with delta %g\n", info->numz, info->zmin, info->zmax, info->zdelta ) ;
      gpterror( "%s: No rectangular grid with the above parameters at (%g,%g,%g)\n", filename, info->x[N], info->y[N], info->z[N] ) ;
    }
  }

  gptaddEBelement( init, map3d_sim, map3d_exit, GPTELEM_GLOBAL, info ) ;
}


static int map3d_sim(gptpar *par, double t, struct map3d_info *info)
{
  int i,j,k,N ;
  double ft,fu,fv, gt,gu,gv ;
  int N1,N2,N3,N4,N5,N6,N7,N8 ;
  double V1,V2,V3,V4,V5,V6,V7,V8 ;

  if( Z<info->zmin || Z>info->zmax ) return( 0 ) ;
  if( Y<info->ymin || Y>info->ymax ) return( 0 ) ;
  if( X<info->xmin || X>info->xmax ) return( 0 ) ;

  /* Calculate master index and offset fractions */
  i = (int)((X-info->xmin)/info->xdelta) ;
  j = (int)((Y-info->ymin)/info->ydelta) ;
  k = (int)((Z-info->zmin)/info->zdelta) ;
  N = (k*info->numy+j)*info->numx+i ;
  ft= (X-info->xmin-i*info->xdelta)/info->xdelta ;
  fu= (Y-info->ymin-j*info->ydelta)/info->ydelta ;
  fv= (Z-info->zmin-k*info->zdelta)/info->zdelta ;
  gt = 1-ft ;
  gu = 1-fu ;
  gv = 1-fv ;

  /* Calculate boundary offsets */
  N1 = N ;
  N2 = N+1 ;
  N3 = N+info->numx+1 ;
  N4 = N+info->numx ;
  N5 = N+info->numx*info->numy ;
  N6 = N+info->numx*info->numy+1 ;
  N7 = N+info->numx*(info->numy+1)+1 ;
  N8 = N+info->numx*(info->numy+1) ;

/* Consistency check in case of problems with this element
 * if(i<0 || i>=info->numx-1) gpterror( "Error in x range: %lf maps to %d\n", X, i ) ;
 * if(j<0 || j>=info->numy-1) gpterror( "Error in y range: %lf maps to %d\n", Y, j ) ;
 * if(k<0 || k>=info->numz-1) gpterror( "Error in z range: %lf maps to %d\n", Z, k ) ;
 * if(ft<0 || ft>1+16*DBL_EPSILON ) gpterror( "x fraction is %.16f at x=%f\n", ft, X ) ;
 * if(fu<0 || fu>1+16*DBL_EPSILON ) gpterror( "y fraction is %.16f at y=%f\n", fu, Y ) ;
 * if(fv<0 || fv>1+16*DBL_EPSILON ) gpterror( "z fraction is %.16f at z=%f\n", fv, Z ) ;
 */

  V1 = info->V[N1] ; V2 = info->V[N2] ; V3 = info->V[N3] ; V4 = info->V[N4] ;
  V5 = info->V[N5] ; V6 = info->V[N6] ; V7 = info->V[N7] ; V8 = info->V[N8] ;
  EX = (gv*(gu*(V1-V2) - fu*(V3-V4)) + fv*(gu*(V5-V6) - fu*(V7-V8))) / info->xdelta ;
  EY = (gv*(gt*(V1-V4) + ft*(V2-V3)) + fv*(gt*(V5-V8) + ft*(V6-V7))) / info->ydelta ;
  EZ = (gt*(gu*(V1-V5) + fu*(V4-V8)) + ft*(gu*(V2-V6) + fu*(V3-V7))) / info->zdelta ;

  return( 1 ) ;
}

static void map3d_exit(struct map3d_info *info)
{
  gdfmfreechilds( &info->gm.ds ) ;
  free(info) ;
}


/* swap two points, keeping all arrays in sync */
static void swapmap3D(struct map3d_info *info, int i, int j)
{
  if( i!=j &&
      info->x[i]==info->x[j] &&
      info->y[i]==info->y[j] &&
      info->z[i]==info->z[j] )
    gpterror( "Duplicate points at (%g,%g,%g)\n", info->x[i], info->y[i], info->z[i] ) ;

  swap(info->x[i],info->x[j]) ;
  swap(info->y[i],info->y[j]) ;
  swap(info->z[i],info->z[j]) ;

  swap(info->V[i],info->V[j]) ;
}


/* Sort all datapoints in correct order: x runs fastest */
static void qsortmap3D(struct map3d_info *info, int left, int right)
{
  int i, last ;

  if( left>=right ) return ;

  swapmap3D(info,left,(left+right)/2) ;
  last = left ;
  for(i=left+1 ; i<=right ; i++)
    if(info->z[i]<info->z[left] ||
      (info->y[i]<info->y[left] && info->z[i]==info->z[left]) ||
      (info->x[i]<info->x[left] && info->z[i]==info->z[left] && info->y[i]==info->y[left]))
        swapmap3D(info,++last,i) ;
  swapmap3D(info,left,last) ;
  qsortmap3D(info,left,last-1) ;
  qsortmap3D(info,last+1,right) ;
}
