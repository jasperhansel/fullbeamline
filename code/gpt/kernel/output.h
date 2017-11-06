/* output.h - Output header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

struct gptoutvector
{
  double vec[3] ;
  int alive ;
} ;

struct gptoutscalar
{
  double scalar ;
  int alive ;
} ;

/* Should be used by user functions to avoid extensive memory mismanagement */
extern struct gptoutvector *gptoutvectors ;
extern struct gptoutscalar *gptoutscalars ;

void outcoords( gpttransform *tf, double t, double *dt, double *x, bool fields ) ;

/* Support for tout hooks */
typedef void (*toutfnc)(double t, double *dt, double *x, void *info) ;
void toutaddfncU( toutfnc fnc, void *info ) ;
#define toutaddfnc(a,b) toutaddfncU((toutfnc)(a),b)

/* Support for screen hooks */
struct screenparinfo
{
  double time ;
  double Wr[3] ;
  double GBr[3] ;

  double m, q ;
  double rmacro, n ;
  double ID ;

  bool alive ;
} ;

typedef void (*screenfnc)(struct screenparinfo *pars, int len, void *info) ;
void screenaddfnc(screenfnc fnc, void *info) ;

/* Public output routines */
void gpswparr( char *name, double *x ) ;
void gptoutvector(char *xname, char *yname, char *zname) ;
void gptoutscalar(char *name) ;

/* Interface between gpt and gdf */
void gptoutputdoublearray( char *name, double *x, int len ) ;
void gptoutputdouble( char *name, double x ) ;
void gptoutputdoublegroup( char *name, double x ) ;
void gptoutputendgroup( void ) ;
