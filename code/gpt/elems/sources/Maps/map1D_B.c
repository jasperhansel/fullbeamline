/* map1D_B.c - Reads a 1D table of axial Bz samples and extrapolates these
 *             to a cylindrical symmetric field map for the magnetic field
 */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "elem.h"


struct map1d_info
{
  /* GDF information */
  struct gdfmem gm ;

  /* Coordinate and field information */
  double *z, *fz, *fz2 ;
  int points ;
} ;


static int map1d_sim(gptpar *par, double t, struct map1d_info *info) ;
static void map1d_exit(struct map1d_info *info) ;


void map1D_B_init(gptinit *init)
{
  struct map1d_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  char *zname ;
  char *fzname ;
  double ffac ;

  /* Grid characteristics */
  int i,points ;

  /* GDF information */
  struct gdff ingdff ;
  char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,z,Bz,Bfac)\n", gptgetname(init) ) ;

  info = (struct map1d_info *)gptmalloc( sizeof(struct map1d_info) ) ;

  /* Initialize GDF */
  filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  /* Retrieve z array */
  zname = gptgetargstring(init,2) ;
  ds=gptinputgetgroup(&info->gm,filename,zname) ;
  info->z=gptinputdoublearray(ds,zname,&points) ;
  info->points = points ;

  /* Retrieve field from datafile */
  fzname   = gptgetargstring(init,3) ;
  info->fz = gptinputdoublearraypoints(ds,fzname,points) ;

  /* Multiply field with extra factor */
  ffac = gptgetargdouble(init,4) ;
  for(i=0 ; i<points ; i++) info->fz[i] *= ffac ;

  /* Check ascending z-order */
  for( i=0 ; i<(points-1) ; i++ )
    if( info->z[i]>=info->z[i+1] )
      gpterror( "%s: z array must be in ascending order\n", filename ) ;

  /* Calculate spline information */
  info->fz2 = (double *)gptmalloc(points*sizeof(double)) ;
  if( spline( info->z, info->fz, info->fz2, points, NATSPLINE, NATSPLINE ) )
    gpterror( "%s: Error calculating spline\n", filename ) ;

  gptaddEBelement( init, map1d_sim, map1d_exit, GPTELEM_GLOBAL, info ) ;
}


static int map1d_sim(gptpar *par, double t, struct map1d_info *info)
{
  double fz, dfz, ddfz ;
  double Bror ;

  if( splintdd( info->z, info->fz, info->fz2, info->points, Z,
    &fz, &dfz, &ddfz ) ) return( 0 ) ;

  Bror =-dfz/2 ;
  BX   = X*Bror ;
  BY   = Y*Bror ;
  BZ   = fz - (X*X+Y*Y)*ddfz/4 ;

  return( 1 ) ;
}

static void map1d_exit(struct map1d_info *info)
{
  free(info->fz2) ;
  gdfmfreechilds( &info->gm.ds ) ;
  free(info) ;
}

