/* spacecharge3Dtree.c - Barnes&Hut relativistic treecode */

#define _CRT_SECURE_NO_DEPRECATE

#include <vector>

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#include "../../kernel/BHtree.h"

/* Forward declarations */
static void sc_sim( double t, double *x, double *p, void *vinfo, int nOff, int nCPU ) ;
static void sc_tout( double t, double *dt, double *x, void *info ) ;
static void sc_screen(struct screenparinfo *pars, int len, void *info) ;

struct pointchargeset_info
{
  auto_ptr<electrostaticfield> BHtree ;
} ;

void pointchargeset_init(gptinit *init)
{
// Command-line parameters
  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,theta)\n", gptgetname(init) ) ;

  char *name   = gptgetargstring(init,1) ;
  double theta = gptgetargdouble(init,2) ;

// Get particle set
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set %s does not exist\n", name ) ;
  gptparset *set = gptgetparset( name ) ;
  int len ;
  gptinitpar *par = gptgetparsetpars( set,&len ) ;

// No particles, nothing to do
  if( len==0 )
  {
	  gptwarning( "%s: No particles in set \"%s\". Element ignored.\n", gptgetname(init), name ) ;
	  return ;
  }

// Get rmacro (must be identical for all particles for B&H efficiency reasons)
  double eps = par[0].r ; // OK, len!=0 as just checked above
  for(int i=1 ; i<len ; i++)
  {
    if( par[i].r != eps )
		gpterror( "%s: Particle radii must be identical: Use 'setrmacrodist(\"%s\",\"u\",rmacro,0)'\n", gptgetname(init), name ) ;
  }
  if( eps==0 )
  {
	  gptwarning( "%s: Particles have zero radius resulting in singularities in the field.\n", gptgetname(init) ) ;
	  gptwarning( "%s: Use 'setrmacrodist(\"%s\",\"u\",rmacro,0)' with positive rmacro.\n",  gptgetname(init), name ) ;
  }

// Create particle array
  pointchargearray charges(eps*eps) ;

  for(int i=0 ; i<len ; i++)
  {
    vec3 Wr(par[i].Wr) ;
    pointcharge q(Wr,par[i].n*par[i].q) ;
    charges.insert(q) ;
  }

// Dump particles
  gptremoveparset(set) ;

// Create BHtree
  struct pointchargeset_info *info = new(pointchargeset_info) ;
  info->BHtree = BHnode::createBHelectrostaticfield(charges,theta) ;

// Register callback functions
  odemtaddfprfunction( ODEFNC_USR, sc_sim,info ) ;
  toutaddfnc( sc_tout, info ) ;
  screenaddfnc( sc_screen, info ) ;
}

/* Main routine to calculate B&H interactions */
static void sc_sim( double t, double *x, double *p, void *vinfo, int nOff, int nCPU )
{
  struct pointchargeset_info *info = (struct pointchargeset_info *)vinfo ;

  for(int i=nOff ; i<numpar ; i+=nCPU) if( pars[i].alive )
  {
    vec3 Wr(pars[i].Wr[0],pars[i].Wr[1],pars[i].Wr[2]) ;
    vec3 Ebem = info->BHtree->getEfield(Wr) ;
    pars[i].WE[0] += Ebem[0] ;
    pars[i].WE[1] += Ebem[1] ;
    pars[i].WE[2] += Ebem[2] ;
  }
}

// Add kinetic energy
template<typename _pars>
static void Ekin_min_V(_pars particles, int len, std::vector<double> &V)
{
  int i, j ;

  for(i=j=0 ; i<len ; i++) if( particles[i].alive )
  {
    double G    = sqrt(1.0+gptVECSQR(particles[i].GBr)) ;
	double Ekin = (G-1.0)*particles[i].m*gpt_c*gpt_c/(-gpt_qe) ;
	double Vpar = V[j]*particles[i].q/(-gpt_qe) ;
	V[j++] = Ekin + Vpar ;
  }
}

/* Tout output function, outputs electrostatic potential due to the BEM particles */
static void sc_tout( double t, double *dt, double *x, void *vinfo )
{
// No particles, no output
  if( numalivepar==0 ) return ;

// Init
  struct pointchargeset_info *info = (struct pointchargeset_info *)vinfo ;
  std::vector<double> V ;
  V.reserve(numalivepar) ;

// Get potentials
  for(int i=0 ; i<numpar ; i++) if( pars[i].alive )
	V.push_back(info->BHtree->getV(vec3(pars[i].Wr))) ;

// Output
  gptoutputdoublearray("Vbem",&V[0],numalivepar) ;

// Kinetic energy minus BEM potential
  Ekin_min_V(pars,numpar,V) ;
  gptoutputdoublearray("EkinVbem",&V[0],numalivepar) ;
}

/* Screen output function, outputs electrostatic potential due to the BEM particles */
void sc_screen(struct screenparinfo *pars, int len, void *vinfo)
{
 // No particles, no output
  if( len==0 ) return ;

// Init
  struct pointchargeset_info *info = (struct pointchargeset_info *)vinfo ;
  std::vector<double> V ;
  V.reserve(len) ;

// Get potentials
  for(int i=0 ; i<len ; i++) if( pars[i].alive )
	V.push_back(info->BHtree->getV(vec3(pars[i].Wr))) ;

// Combined into V1
  gptoutputdoublearray("Vbem",&V[0],len) ;

// Kinetic energy minus BEM potential
  Ekin_min_V(pars,len,V) ;
  gptoutputdoublearray("EkinVbem",&V[0],len) ;
}
