/* setfile.c - Add a group in a GDF file to a particle set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "elem.h"

#define BUFSIZE 65536

static char axisname[] = "wcs" ;


void setfile_init(gptinit *init)
{
  struct gdff gdff ;
  struct gdfmem gm ;
  struct gdfdata *ds ;
  char *name, *filename ;
  char groupname[256] ;
  gptparset *set ;
  double *x, *y, *z, *Bx, *By, *Bz, *G, *GBx, *GBy, *GBz ;
  double *m, *q, *nmacro, *rmacro, *tstart, *IDs ;
  double B2, defm, defq, defn, defr, *pdeft, deft, gamma ;
  int allB, allGB ;
  int i, numpar ;
  struct axis *paxis ;
  double r[3], GBr[3] ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,filename)\n", gptgetname(init) ) ;

  name     = gptgetargstring(init,1) ;
  filename = gptgetargstring(init,2) ;

  set = gptgetparset( name ) ;

  /* Open and read GDF file */
  if( (gdff.fp=fopen( filename, "rb" ))==NULL )
    errdisp( gpsem, GPSINPERR_FE, filename, strerror( errno ) ) ;
  setvbuf(gdff.fp,NULL,_IOFBF,BUFSIZE);
  gdfrmem( &gdff, &gm, GDFR_ABORTONERROR ) ;

  /* Select correct group */
  ds=gptinputgetgroup(&gm,filename,"x") ;

  /* Get number of particles */
  x  = gptinputdoublearray( ds, "x", &numpar ) ;

  /* Obtain mandatory coordinates */
  x  = gptinputdoublearraypoints( ds, "x", numpar ) ;
  y  = gptinputdoublearraypoints( ds, "y", numpar ) ;
  z  = gptinputdoublearraypoints( ds, "z", numpar ) ;

  /* Obtain optional coordinates */
  G      = gptinputdoublearraypointsnull( ds, "G", numpar ) ;
  Bx     = gptinputdoublearraypointsnull( ds, "Bx", numpar ) ;
  By     = gptinputdoublearraypointsnull( ds, "By", numpar ) ;
  Bz     = gptinputdoublearraypointsnull( ds, "Bz", numpar ) ;

  GBx    = gptinputdoublearraypointsnull( ds, "GBx", numpar ) ;
  GBy    = gptinputdoublearraypointsnull( ds, "GBy", numpar ) ;
  GBz    = gptinputdoublearraypointsnull( ds, "GBz", numpar ) ;

  m      = gptinputdoublearraypointsnull( ds, "m", numpar ) ;
  q      = gptinputdoublearraypointsnull( ds, "q", numpar ) ;
  nmacro = gptinputdoublearraypointsnull( ds, "nmacro", numpar ) ;
  rmacro = gptinputdoublearraypointsnull( ds, "rmacro", numpar ) ;
  tstart = gptinputdoublearraypointsnull( ds, "t", numpar ) ;

  IDs    = gptinputdoublearraypointsnull( ds, "ID", numpar ) ;

  gdfmgettitle(groupname,ds) ;
  errdisp( gpsem, GPSVERB_START_FILE, numpar, filename, groupname ) ;

  /* Test if all required coordinates are present */
  allB = allGB = 0 ;
  if(  Bx!=NULL &&  By!=NULL &&  Bz!=NULL ) allB  = 1 ;
  if( GBx!=NULL && GBy!=NULL && GBz!=NULL ) allGB = 1 ;

  if( allGB && G!=NULL )
    gptwarning( "%s: G array is ignored\n", filename ) ;

//if( (allB==0 && allGB==0) )
//  gpterror( "%s: Cannot load Bx,y,z or GBx,y,z array\n", filename ) ;

  /* Store defaults */
  gptgetvardouble("m", &defm ) ;
  gptgetvardouble("q", &defq ) ;

  defn = defr = 0.0 ;
  gptgetvardouble("nmacro", &defn ) ;
  gptgetvardouble("rmacro", &defr ) ;

  deft = 0.0 ;
  if( (pdeft=gdfmgetdbl(ds,"time"))!=NULL )
  {
	deft=*pdeft ;
	gptinstalldouble("time",deft) ;
  }

  paxis = getaxis( axisname ) ;

  /* Add particles */
  for(i=0 ; i<numpar ; i++)
  {
    r[0] = x[i] ;
    r[1] = y[i] ;
    r[2] = z[i] ;

    if( allGB )
    {
      GBr[0] = GBx[i] ;
      GBr[1] = GBy[i] ;
      GBr[2] = GBz[i] ;
    } else if( allB )
    {
      B2 = Bx[i]*Bx[i]+By[i]*By[i]+Bz[i]*Bz[i] ;
      if( B2 > 1 )
        gpterror( "%s: Velocity of particle %d (%g,%g,%g)*c is larger than the speed of light\n", filename, i, Bx[i], By[i], Bz[i] ) ;

      gamma = (G!=NULL ? G[i] : 1.0/sqrt(1-B2) ) ;
      GBr[0] = gamma*Bx[i] ;
      GBr[1] = gamma*By[i] ;
      GBr[2] = gamma*Bz[i] ;
    } else
    {
      GBr[0] = 0 ;
      GBr[1] = 0 ;
      GBr[2] = 0 ;
    }

    gptaddparmqnartid( set, r, GBr,
		    m     !=NULL ? m[i]      : defm,
		    q     !=NULL ? q[i]      : defq,
		    nmacro!=NULL ? nmacro[i] : defn,
            paxis,
            rmacro!=NULL ? rmacro[i] : defr,
			tstart!=NULL ? tstart[i] : deft,
			IDs   !=NULL ? (int)(IDs[i]+0.5) : 0 ) ;
  }

  /* Free memory */
  gdfrclose(&gdff) ;
  gdfmfreechilds( &gm.ds ) ;

  return ;
}
