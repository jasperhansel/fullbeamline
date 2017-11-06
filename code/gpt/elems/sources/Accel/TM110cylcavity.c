/* TM110cylcavity.c - Cylindrical TM110 mode resonant cavity */

/* Jacco Nohlmans, Bas van der Geer
 * 29-sep-06
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define x11 3.83170597020751231561

struct tmcav_info
{
  double d, radius ;
  double fac, phi ;
  double w110 ;
} ;

static int tmcav_sim(gptpar *par, double t, struct tmcav_info *info) ;


void TM110cylcavity_init(gptinit *init)
{
  struct tmcav_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,d,R,const,phi)\n", gptgetname(init) ) ;

  info = (struct tmcav_info *)gptmalloc( sizeof(struct tmcav_info) ) ;

  info->d      = gptgetargdouble(init,1)/2 ;
  info->radius = gptgetargdouble(init,2) ;
  info->fac    = gptgetargdouble(init,3) ;
  info->phi    = gptgetargdouble(init,4) ;

  if( info->d<=0 || info->radius<=0 )
    gpterror( "Cavity dimensions must be larger than zero\n" ) ;

  if( info->fac==0 )
    gptwarning( "Const is equal to zero, no power in cavity\n" ) ;

  info->w110 = x11*gpt_c/info->radius ;

  gptaddEBelement( init, tmcav_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int tmcav_sim(gptpar *par, double t, struct tmcav_info *info)
{
    double r ;
    double fac, radius, w110, phi ;
    double Br,Bphi ;  
    
    if( fabs(Z) > info->d ) return( 0 ) ;
    r = sqrt(X*X + Y*Y) ;
    radius = info->radius ;
    if( r > radius ) return( 0 ) ;
    
    fac    = info->fac ;
    w110   = info->w110 ;
    phi    = info->phi ;
    
    if( x11*r/radius < 1e-6 )
    {
        double xovr = X*x11/radius ;
        double yovr = Y*x11/radius ;
        
        EZ = fac*xovr*sin(w110*t + phi)/2 ;
        BX = fac*(-xovr*yovr*cos(w110*t + phi))/(8*gpt_c) ;
        BY = fac*(4*(yovr*yovr-8) - xovr*xovr*(yovr*yovr-12))*cos(w110*t+phi)/(64*gpt_c) ;
    } else
    {
        double sinxy = Y/r ; /* sin(angle in xy-plane) */
        double cosxy = X/r ;
        
        EZ   =  fac * cosxy * sin(w110*t + phi) * bessj1(x11*r/radius) ;
        Br   = -fac * sinxy * cos(w110*t + phi) * bessj1(x11*r/radius) / (w110*r) ;
        Bphi = -fac * cosxy * cos(w110*t + phi) *(bessj0(x11*r/radius)*x11/radius - bessj1(x11*r/radius)/r)/w110 ;
        
        gptrphi2carth(Br,Bphi,X,Y,&BX,&BY) ;
    }
    
    return( 1 ) ;
}
