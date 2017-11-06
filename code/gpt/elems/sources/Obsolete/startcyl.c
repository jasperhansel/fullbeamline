/* startcyl.c - Start with cilindrical beam 
 * =========
 * Derived from genpar7, Class C modulation included
 *
 * This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"


struct Par
{
  double x, y, z, bx, by, bz, E;
} ;

void startcyl_init(gptinit *init)
{
  struct Par *Pars ;
  int    i, nps, mid;
  int numarg ;
  double E0, sigE, R, Divspr, Div, epsn, gam, betaz, sigz, psi, omega=1.0 ;
  double Erest ;

  /* modulate */    double dz, z, ddE, dE, kbeam=1. ;
                    char   Emodtype, Dmodtype ;

  /* fill rings */  int	    nr, nrs, npr, np ; 
                    double  phi	  ;

  /* shuffle */     int j,k ;
                    double t ;
  /* correct bz */  double bet2, betxy2, betz=0. ;
  /* addparr */     double r[3], Br[3] ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=9 && numarg!=10 )
    gpterror( "Syntax: %s(ECS,nps,E0,Emod,sigE,R,div,divspread,Densmod,sigz/psi,[omega]);\n",
               gptgetname(init) ) ;

  /* Read parameters */
  nps      = gptgetargint(init,1) ;
  E0       = gptgetargdouble(init,2);
  Emodtype = toupper(gptgetargstring(init,3)[0]) ;
  sigE     = gptgetargdouble(init,4) ;
  R        = gptgetargdouble(init,5) ;
  Div      = gptgetargdouble(init,6) ;
  Divspr   = gptgetargdouble(init,7) ;
  Dmodtype = toupper(gptgetargstring(init,8)[0]) ;
  sigz=psi = gptgetargdouble(init,9) ; 

/* Calculate Erest=|m c^2/q| */
  {
    double m, q ;
    gptgetvardouble("m",&m) ;
    gptgetvardouble("q",&q) ;
    Erest = fabs(m*gpt_c*gpt_c/q) ;
  }

  gam	   = 1.0 + E0/Erest ;
  betaz    = sqrt(1.0 - 1.0/(gam*gam) ) ;
  epsn     = betaz*gam*gpt_pi*R*Divspr;   

  if(numarg==10)
  {
    omega    = gptgetargdouble(init,10) ;
    kbeam    = omega/(betaz*gpt_c) ; /* radians gun per meter beam for Class C */
  }

/* Guard the algorithm */
  if( nps <= 0 || E0 <=0.0 || sigE<0.0 || R<0.0  || Divspr<0.0
               || sigz<0.0 || omega<0. || psi<0. || psi>gpt_pi )
    gpterror( "Impossible input parameters in %s \n",  gptgetname(init) ) ;

 
/* reduce nps to fit whole rings */
   nrs = 1 ;
   i = 1 ;
   while ( (i + 6*nrs) <= nps )
   {
      i += 6*nrs ;
      nrs++ ;
   }
   nps = i ;


/* allocate memory */
   Pars=(struct Par *)gptmalloc( nps * sizeof(struct Par) ) ;

/*  Density and Energy modulation
 *  Sigma is interpreted as total bunchlength reps total spread in case "Lin"
 *  and as half conduction angle, psi, for "Class_C"
 */
 
/* Density modulation */
   z = dz = 0.0 ;
   mid = nps/2 ;
   for ( i = 0 ;  i <= mid ;  i++ )
   {
      Pars[mid-i].z = -z ;
      Pars[mid+i].z =  z ;

      if (nps>1) switch( Dmodtype )
      {
        case 'U' : dz = sigz/nps ; break ;
        case 'G' : dz = ( sigz*sqrt(2.0*gpt_pi)*exp( z*z/(2.0*sigz*sigz)) )/nps ;  break ;
        case 'C' : dz = (1./kbeam)*2.0*( sin(psi) - psi*cos(psi) ) / (nps*( cos(kbeam*z)-cos(psi))); break ;
        default  : gpterror( "Dmod in %s should be \"Uniform\", \"Gaussian\" or \"Class C\"\n", gptgetname(init) ) ;
      }
      z += dz ;
   }
   
    
/* Energy modulation */
   sigE *= E0 ;
   dE = ddE  = 0.0 ;
   for( i = 0 ; i <= mid ; i++)
   {
      Pars[mid-i].E = E0 - dE ;
      Pars[mid+i].E = E0 + dE ;

      if (nps>1) switch( Emodtype )
      {
        case 'U' : ddE = sigE/nps ; break ;
        case 'G' : ddE = sigE*sqrt(2.0*gpt_pi)*exp(dE*dE/(2.0*sigE*sigE)) /nps ; break ;
        default: gpterror( "Emod in %s should be \"Uniform\" or \"Gaussian\"\n", gptgetname(init) ) ;
      }
      dE += ddE ;
   }

   for(i=0 ; i<nps ; i++)
   {
     gam = 1.0 + Pars[i].E/Erest ;
     Pars[i].bz = sqrt(1.0 - 1.0/(gam*gam));
   }


/* fill rings */
   np = 0 ;
   for (nr = 0 ;  nr < nrs ; nr++ )
   {
      double fac = nrs/(nrs-0.7) ; /* Empirical correction function, Pulsar */
      npr=0 ;
      do
      {  
         phi	     = (nr==0) ? 0.0 : npr*2*gpt_pi/(6*nr) ;
	 Pars[np].x  = ((nrs==1) ? 0.0 : fac*nr*R*cos(phi)/nrs) ;
	 Pars[np].y  = ((nrs==1) ? 0.0 : fac*nr*R*sin(phi)/nrs) ;
	 Pars[np].bx = ((nrs==1) ? 0.0 : fac*nr*Divspr*betaz*cos(phi)/nrs) ;
	 Pars[np].by = ((nrs==1) ? 0.0 : fac*nr*Divspr*betaz*sin(phi)/nrs) ;
	 np++ ;
	 npr++ ;
      } while( npr<6*nr ) ;
    }


/* shuffle */
  for(i=0 ; i<10*nps ; i++)
  {
    j = rand()%nps ; k = rand()%nps ;
    t = Pars[j].x  ; Pars[j].x  = Pars[k].x  ; Pars[k].x  = t ;
    t = Pars[j].y  ; Pars[j].y  = Pars[k].y  ; Pars[k].y  = t ;
  }

  for(i=0 ; i<10*nps ; i++)
  { 
    j = rand()%nps ; k = rand()%nps ;
    t = Pars[j].bx ; Pars[j].bx = Pars[k].bx ; Pars[k].bx = t ;
    t = Pars[j].by ; Pars[j].by = Pars[k].by ; Pars[k].by = t ;
  }

  for(i=0 ; i<10*nps ; i++)
  {
    j = rand()%nps ; k = rand()%nps ;
    t = Pars[j].bz ; Pars[j].bz = Pars[k].bz ; Pars[k].bz = t ;
    t = Pars[j].E  ; Pars[j].E  = Pars[k].E  ; Pars[k].E  = t ;
  }


/* add divergence */
  for (np = 0 ;  np < nps ; np++ )
  {
     Pars[np].bx += (R>0? Div*betaz*Pars[np].x/R : Div*betaz) ;
     Pars[np].by += (R>0? Div*betaz*Pars[np].y/R : Div*betaz) ;
  }


/* correct bz such that E is total energy */
  for(i=0 ; i<nps ; i++)
  {
       gam = 1.0 + fabs( Pars[i].E) / Erest ;
       bet2 = 1.0 - 1.0/(gam*gam);
       betxy2 = Pars[i].bx*Pars[i].bx + Pars[i].by*Pars[i].by ;
       if ( bet2 <= betxy2 )
       {
         Pars[i].bx *= sqrt(bet2/betxy2) ;
         Pars[i].by *= sqrt(bet2/betxy2) ;
         Pars[i].bz  = 0.0 ;
       }
       else
       {
         betz = sqrt( bet2-betxy2 );
         Pars[i].bz = betz * (Pars[i].E >= 0.0 ? 1.0 : -1.0 ) ;
       }
  }


/* Add particles to kernel */
  for(i=0; i<nps; i++)
  {
    r[0]  = Pars[i].x  ;
    r[1]  = Pars[i].y  ;
    r[2]  = Pars[i].z  ;
    Br[0] = Pars[i].bx ;
    Br[1] = Pars[i].by ;
    Br[2] = Pars[i].bz ;
   
    gptaddparticle( init, r ,Br ) ;
  }

/* Free memory */
  gptfree(Pars) ;
}
