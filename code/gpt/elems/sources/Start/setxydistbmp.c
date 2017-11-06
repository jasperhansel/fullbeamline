/* setxydistbmp.c - Set xy-distribution based on bitmap file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "elem.h"

void setxydistbmp_init(gptinit *init)
{
  char *setname ;
  gptparset *set ;
  gptinitpar *par ;
  char *filename ;

  int argnum, i, len ;
  struct grayscalebmp bmp ;
  unsigned char maxpixel ;

  double xrand, yrand, zrand ;

/* Command-line and read bitmap file */
  argnum = gptgetargnum(init) ;
  if( argnum!=2 && argnum!=4 )
    gpterror( "Syntax: %s(set,filename,[xres,yres])\n", gptgetname(init) ) ;

  setname  = gptgetargstring(init,1) ;
  filename = gptgetargstring(init,2) ;

  readgrayscalebmp(&bmp,filename) ;

/* No black bitmaps allowed */
  maxpixel = 0 ;
  for(i=0 ; i<bmp.xpels*bmp.ypels ; i++) if( bmp.pixels[i]>maxpixel ) maxpixel=bmp.pixels[i] ;
  if(maxpixel==0) gpterror( "%s: Bitmap is entirely black. No particles can be started.\n", filename ) ;

/* Optionally overwerite pixel size */
  if( argnum==4 )
  {
    bmp.xres = gptgetargdouble(init,3) ;
    bmp.yres = gptgetargdouble(init,4) ;
  }
  
/* Get particle set */
  if( gpttestparset( setname )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", setname ) ;
  set = gptgetparset( setname ) ;
  par = gptgetparsetpars( set,&len ) ;

/* Verbose output */
  gptverbose( "Setting xy-distribution based on %dx%d bitmap \"%s\" with a maximum pixel value of %d\n", bmp.xpels, bmp.ypels, filename, maxpixel ) ;

/* Set distribution, rejection method */
  for( i=0 ; i<len ; i++ )
  {
    xrand = bmp.xpels * dblpprand() ;
    yrand = bmp.ypels * dblpprand() ;
    zrand = ((int)maxpixel+1) * dblpprand() ;

    if( bmp.pixels[(((int)yrand)*bmp.xpels + (int)xrand)] > zrand )
    {
      par[i].Wr[0] = (xrand - 0.5*bmp.xpels)/bmp.xres ;
      par[i].Wr[1] = (yrand - 0.5*bmp.ypels)/bmp.yres ;
    } else i-- ;
  }

/* Cleanup */
  gptfree(bmp.pixels) ;
}