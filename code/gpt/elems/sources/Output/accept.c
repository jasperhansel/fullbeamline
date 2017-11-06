/* accept.c: Determine beamline acceptance */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


#define CHUNKCOUNT 65536
enum state { PARNEW, PARKIL } ;
struct buffparinfo
{
  enum state state ;
  double  Wr[3],  GBr[3],  G,  time ; /* Start position */
  double kWr[3], kGBr[3], kG, ktime ; /* Kill  position */
} ;

struct acceptance_info
{
  struct buffparinfo *buff ;
  unsigned int lastpar ;
  gpttransform outtf ; /* Transformation to output particles */
} ;

static int  acceptance_out_usr( double t, double *dt, double *x, void *vinfo ) ;
static void acceptance_end_usr( double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo ) ;
static void acceptance_ter( void *vinfo ) ;


void acceptance_init(gptinit *init)
{
  struct acceptance_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=0 )
    gpterror( "Syntax: %s(ECS)\n", gptgetname(init) ) ;

  /* Allocate and initialize info structure */
  info = (struct acceptance_info *)gptmalloc( sizeof(struct acceptance_info) ) ;
  info->buff    = (struct buffparinfo *)gptmallocitemarray(sizeof(struct buffparinfo),CHUNKCOUNT,1) ;
  info->lastpar = 0 ;
  gptconcattransform( &info->outtf, &init->paxis->a, &init->e ) ;

/* Register ODE functions */
  odeaddoutfunction( ODEFNC_USR, acceptance_out_usr, info ) ; /* ODEFNC_USR correct ?? */
  odeaddendfunction( ODEFNC_USR, acceptance_end_usr, info ) ; /* ODEFNC_USR correct ?? */
  gptaddmainfunction( GPTMAINFNC_TER, acceptance_ter, info ) ;
}


static int  acceptance_out_usr( double t, double *dt, double *x, void *vinfo)
{
  struct acceptance_info *info = (struct acceptance_info *)vinfo ;

  for(int i=0 ; i<numpar ; i++)
  {
    /* Store new particle */
    if( (unsigned int)pars[i].ID > info->lastpar )
    {
      info->buff = (struct buffparinfo *)gptmallocitem(info->buff,sizeof(struct buffparinfo),CHUNKCOUNT,info->lastpar) ;

      gpttoUCS         (&info->outtf,pars[i].Wr ,info->buff[pars[i].ID-1].Wr ) ;
      gptdirectiontoUCS(&info->outtf,pars[i].GBr,info->buff[pars[i].ID-1].GBr) ;
      info->buff[pars[i].ID-1].G     = pars[i].G ;
      info->buff[pars[i].ID-1].time  = t ;
	  info->buff[pars[i].ID-1].state = PARNEW ;

      info->lastpar = pars[i].ID ;
    }    
  }

  return(0) ;
}

static void acceptance_end_usr( double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo )
{
  struct acceptance_info *info = (struct acceptance_info *)vinfo ;

  for(int i=0 ; i<numpar ; i++)
  {
    /* Mark deleted particles, and store position */
    if( pars[i].tokill )
    {
      gpttoUCS         (&info->outtf,pars[i].Wr ,info->buff[pars[i].ID-1].kWr ) ;
      gptdirectiontoUCS(&info->outtf,pars[i].GBr,info->buff[pars[i].ID-1].kGBr) ;
      info->buff[pars[i].ID-1].kG    = pars[i].G ;
      info->buff[pars[i].ID-1].ktime = tend ;
	  info->buff[pars[i].ID-1].state = PARKIL ;
    }
  }
}

static void acceptance_ter( void *vinfo )
{
  struct acceptance_info *info = (struct acceptance_info *)vinfo ;
  unsigned int i, j ;
  double *arr ;

  arr = (double *)gptmalloc(info->lastpar*sizeof(double)) ;

  gptoutputdoublegroup( "acceptance", 0 ) ;

  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].Wr[0] ;
  gptoutputdoublearray("OK_x_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].Wr[1] ;
  gptoutputdoublearray("OK_y_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].Wr[2] ;
  gptoutputdoublearray("OK_z_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].GBr[0]/info->buff[i].G ;
  gptoutputdoublearray("OK_Bx_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].GBr[1]/info->buff[i].G ;
  gptoutputdoublearray("OK_By_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].GBr[2]/info->buff[i].G ;
  gptoutputdoublearray("OK_Bz_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].G ;
  gptoutputdoublearray("OK_G_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARNEW) arr[j++] = info->buff[i].time ;
  gptoutputdoublearray("OK_t_start",arr,j) ;

  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].Wr[0] ;
  gptoutputdoublearray("NOK_x_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].Wr[1] ;
  gptoutputdoublearray("NOK_y_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].Wr[2] ;
  gptoutputdoublearray("NOK_z_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].GBr[0]/info->buff[i].G ;
  gptoutputdoublearray("NOK_Bx_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].GBr[1]/info->buff[i].G ;
  gptoutputdoublearray("NOK_By_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].GBr[2]/info->buff[i].G ;
  gptoutputdoublearray("NOK_Bz_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].G ;
  gptoutputdoublearray("NOK_G_start",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].time ;
  gptoutputdoublearray("NOK_t_start",arr,j) ;

  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kWr[0] ;
  gptoutputdoublearray("NOK_x_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kWr[1] ;
  gptoutputdoublearray("NOK_y_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kWr[2] ;
  gptoutputdoublearray("NOK_z_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kGBr[0]/info->buff[i].kG ;
  gptoutputdoublearray("NOK_Bx_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kGBr[1]/info->buff[i].kG ;
  gptoutputdoublearray("NOK_By_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kGBr[2]/info->buff[i].kG ;
  gptoutputdoublearray("NOK_Bz_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].kG ;
  gptoutputdoublearray("NOK_G_end",arr,j) ;
  for(i=j=0 ; i<info->lastpar ; i++) if(info->buff[i].state==PARKIL) arr[j++] = info->buff[i].ktime ;
  gptoutputdoublearray("NOK_t_end",arr,j) ;

  gptoutputendgroup() ;

  gptfree(arr) ;
  gptfreeitemarray(info->buff,sizeof(struct buffparinfo),CHUNKCOUNT,info->lastpar) ;
  gptfree(info) ;
}
