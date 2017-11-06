/* setcharg.c - Set the total charge of a set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void settotalcharge_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double Q, sum ;
  int neg, zero, pos ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,Q)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  Q    = gptgetargdouble(init,2) ;

  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set %s does not exist\n", name ) ;

  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  sum = 0.0 ;
  neg = zero = pos = 0 ;
  for( i=0 ; i<len ; i++ )
  {
    if( par[i].q < 0 ) neg++ ;
    if( par[i].q == 0 ) zero++ ;
    if( par[i].q > 0 ) pos++ ;
    sum += par[i].q ;
  }
  if( sum==0.0 )
    gpterror( "Total charge of elementary particles sums to zero\n" ) ;
  if( (Q<0 && sum>0) || (Q>0 && sum<0) )
    gpterror( "Sign of total charge of elementary particles is opposite to Q\n" ) ; 

  if( (Q<0 && pos) )
    gptwarning( "%d Elementary particle%s with opposite sign to Q encountered\n", pos, pos>1 ? "s" : "" ) ;  
  if( (Q>0 && neg) )
    gptwarning( "%d Elementary particle%s with opposite sign to Q encountered\n", neg, neg>1 ? "s" : "" ) ;  
  if( zero )
    gptwarning( "%d Elementary particle%s with zero charge encountered\n", zero, zero>1 ? "s" : "" ) ;

  for( i=0 ; i<len ; i++ ) par[i].n = Q/sum ;
}
