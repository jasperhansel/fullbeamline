/* scatterbitmap.c - Scatter function based on a grayscale bitmap */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "elem.h"

#define EXTRAPOLATE 1e-6


struct bmpscat_info
{
  struct grayscalebmp bmp ;
  unsigned char maxpixel ;

  struct scatter_info scatinfo ;
} ;

void scatterbitmap_scat(gptpar *par,double t,double dt, gpttrajectory *trajectory,void *info) ;


void scatterbitmap_init(gptinit *init)
{
  struct bmpscat_info *info ;
  char *name, *filename ;
  int numarg ;
  int i ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=4 )
    gpterror( "Syntax: %s(ECS,name,filename,xres,yres)\n", gptgetname(init) ) ;

  info = (struct bmpscat_info *)gptmalloc( sizeof(struct bmpscat_info) ) ;

  name           = gptgetargstring(init,1) ;
  filename       = gptgetargstring(init,2) ;

  readgrayscalebmp(&info->bmp,filename) ;

/* Overwrite resolution */
  info->bmp.xres = gptgetargdouble(init,3) ;
  info->bmp.yres = gptgetargdouble(init,4) ;

/* No black bitmaps allowed */
  info->maxpixel = 0 ;
  for(i=0 ; i<info->bmp.xpels*info->bmp.ypels ; i++) if( info->bmp.pixels[i]>info->maxpixel ) info->maxpixel=info->bmp.pixels[i] ;
  if(info->maxpixel==0) gpterror( "%s: Bitmap is entirely black. No particles can be started.\n", filename ) ;

  /* Install all functions */
  gptscatterinit(init,&info->scatinfo,name) ;
  gptinstallscatterfnc(name,scatterbitmap_scat,info) ;
}


void scatterbitmap_scat(gptpar *par,double t,double dt, gpttrajectory *ptraj,void *vinfo)
{
  struct bmpscat_info *info = (struct bmpscat_info *)vinfo ;

  double Gnew, extratime ;
  double GBnew2, GBnew[3], vnew[3], rnew[3] ;
  int i ;

  /* Increase output buffer, set Qout and Eout to 0 */
  gptscatterinitpar(&info->scatinfo) ;

  /* Set distribution, rejection method */
  Gnew = ptraj->Gint ;
  GBnew[0] = ptraj->GBint[0] ;
  GBnew[1] = ptraj->GBint[1] ;
  while(1)
  {
    double xrand = info->bmp.xpels * dblpprand() ;
    double yrand = info->bmp.ypels * dblpprand() ;
    double zrand = ((int)info->maxpixel+1) * dblpprand() ;

    if( info->bmp.pixels[(((int)yrand)*info->bmp.xpels + (int)xrand)] > zrand )
    {
      GBnew[0] += (xrand - 0.5*info->bmp.xpels)*info->bmp.xres ;
      GBnew[1] += (yrand - 0.5*info->bmp.ypels)*info->bmp.yres ;
      break ;
    }
  }
  GBnew2 = Gnew*Gnew - GBnew[0]*GBnew[0] - GBnew[1]*GBnew[1] - 1 ;
  if( GBnew2 < 0 )
    gpterror("Too large angle in bitmapscatter element: %e %e sqrt(%e)", GBnew[0], GBnew[1], GBnew2) ;
  GBnew[2] = sqrt(GBnew2) ;

  extratime = (1-ptraj->lambda+EXTRAPOLATE)*dt ;
//  extratime = EXTRAPOLATE*dt ;
//  if( extratime<0 ) extratime=0 ;

  /* Forward scatter */
  for(i=0 ; i<3 ; i++) vnew[i] = gpt_c*GBnew[i]/Gnew ;
  for(i=0 ; i<3 ; i++) rnew[i] = ptraj->P[i]+vnew[i]*extratime ;

  gptscatteraddparmqn(&info->scatinfo, gptgetparset("forwardscatter"), rnew, GBnew, par->m, par->q, par->n ) ;

  /* Remove original particle and exit */
  gptscatterremoveparticle(&info->scatinfo,ptraj,par) ;

  return ;
}
