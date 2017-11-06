/* rectcoil.c - Model for thick current sheet */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define SMALLNZ 1e-4

struct rectcoil_info
{
  double rstart, rend ;
  double a ;
  double I ;
} ;

static int rectcoil_sim(gptpar *par, double t, struct rectcoil_info *info) ;
static double integ0(double zz, double rx) ;
static double integ1(double zz, double rx) ;
static double integ2(double zz, double rx) ;
static double integ3(double zz, double rx) ;
static double integ4(double zz, double rx) ;
static double B0(double z, double a, double rstart, double rend ) ;
static double B1(double z, double a, double rstart, double rend ) ;
static double B2(double z, double a, double rstart, double rend ) ;
static double B3(double z, double a, double rstart, double rend ) ;
static double B4(double z, double a, double rstart, double rend ) ;


void rectcoil_init(gptinit *init)
{
  struct rectcoil_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,r1,r2,L,I)\n", gptgetname(init) ) ;

  info = (struct rectcoil_info *)gptmalloc( sizeof(struct rectcoil_info) ) ;

  info->rstart = gptgetargdouble(init,1) ;
  info->rend   = gptgetargdouble(init,2) ;
  info->a      = gptgetargdouble(init,3)/2 ;
  info->I      = gptgetargdouble(init,4) ;

  if( info->rstart<=0 || info->rend<=0 )
    gpterror( "Coil dimensions must be larger than zero\n" ) ;
  if( info->rend<info->rstart )
    gpterror( "Outer radius of coil must be larger than inner radius\n" ) ;
  if( info->rstart==info->rend )
    gpterror( "Outer radius of coil cannot be equal to inner radius\n" ) ;

  gptaddEBelement( init, rectcoil_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int rectcoil_sim(gptpar *par, double t, struct rectcoil_info *info)
{
  double rstart, rend, a, I ;
  double r2 ;
  double Bror ;    /* Br/r */

  rstart = info->rstart ;
  rend   = info->rend ;
  a      = info->a ;
  I      = info->I ;

/*  if( fabs(Z) > ?????? ) return( 0 ) ; */

  r2 = X*X+Y*Y ;

  BZ   = gpt_mu0*I*(B0(Z,a,rstart,rend) - r2*B2(Z,a,rstart,rend)/4 + 
                                       r2*r2*B4(Z,a,rstart,rend)/64) ;
  Bror = gpt_mu0*I*(-B1(Z,a,rstart,rend)/2 + r2*B3(Z,a,rstart,rend)/16) ;

  BX = Bror*X ;
  BY = Bror*Y ;

  return( 1 ) ;
}



static double integ0(double zz, double rx)
{
  double nz = zz/rx ;
  return( nz*rx*log(rx*(1+sqrt(1+nz*nz))) ) ;
}

static double integ1(double zz, double rx)
{
  double nz = zz/rx ;
  return( 1 - 1/sqrt(1+nz*nz) + log(rx*(1+sqrt(1+nz*nz))) ) ;
}

static double integ2(double zz, double rx)
{
  double nz = zz/rx ;
  if( fabs(nz) < SMALLNZ )
    return( (3*nz*(4-5*nz*nz))/(8*rx) ) ;
  else
    return( (1 - 1/pow(1+nz*nz,1.5)) / (nz*rx) ) ;
}

static double integ3(double zz, double rx)
{
  double nz = zz/rx ;
  if( fabs(nz) < SMALLNZ )
    return( -(3*(-4+15*nz*nz))/(8*rx*rx) ) ;
  else
    return( (-1 + (1+4*nz*nz)/pow(1+nz*nz,2.5)) / (nz*nz*rx*rx) ) ;
}

static double integ4(double zz, double rx)
{
  double nz  = zz/rx ;
  double rx2 = rx*rx ;
  double nz2 = nz*nz ;
  if( fabs(nz) < SMALLNZ )
    return( (5*nz*(-9+35*nz2))/(4*rx2*rx) ) ;
  else
    return( (2 - (2+7*nz2+20*nz2*nz2)/pow(1+nz2,3.5)) / (nz2*nz*rx2*rx) ) ;
}

static double B0(double z, double a, double rstart, double rend )
{
  double tmp, amz, apz ;  

  amz = a-z ;
  apz = a+z ;
  tmp = integ0(amz,rend)+integ0(apz,rend)-integ0(amz,rstart)-integ0(apz,rstart) ;
  return( 1/(4*a*(rend-rstart)) * tmp ) ;    /* *mu0I */
}

static double B1(double z, double a, double rstart, double rend )
{
  double tmp, amz, apz ;  

  amz = a-z ;
  apz = a+z ;
  tmp = -integ1(amz,rend)+integ1(apz,rend)+integ1(amz,rstart)-integ1(apz,rstart) ;
  return( 1/(4*a*(rend-rstart)) * tmp ) ;    /* *mu0I */
}

static double B2(double z, double a, double rstart, double rend )
{
  double tmp, amz, apz ;  

  amz = a-z ;
  apz = a+z ;
  tmp = integ2(amz,rend)+integ2(apz,rend)-integ2(amz,rstart)-integ2(apz,rstart) ;
  return( 1/(4*a*(rend-rstart)) * tmp ) ;    /* *mu0I */
}

static double B3(double z, double a, double rstart, double rend )
{
  double tmp, amz, apz ;  

  amz = a-z ;
  apz = a+z ;
  tmp = -integ3(amz,rend)+integ3(apz,rend)+integ3(amz,rstart)-integ3(apz,rstart) ;
  return( 1/(4*a*(rend-rstart)) * tmp ) ;    /* *mu0I */
}

static double B4(double z, double a, double rstart, double rend )
{
  double tmp, amz, apz ;  

  amz = a-z ;
  apz = a+z ;
  tmp = integ4(amz,rend)+integ4(apz,rend)-integ4(amz,rstart)-integ4(apz,rstart) ;
  return( 1/(4*a*(rend-rstart)) * tmp ) ;    /* *mu0I */
}
