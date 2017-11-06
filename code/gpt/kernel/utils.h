/* B&M: Header file for utilites for:
 *
 * Bessel functions
 * Erf functions
 * Double utilities
 * Random functions
 * Error handling
 * Spline support
 *
 * Remark: Break up when too large or diverse
 */


/* bessel.c */
double bessi0( double x ) ;
double bessi1( double x ) ;
double bessj0(double x) ;
double bessj1(double x) ;

/* erf.c */
double erf( double x) ;
double erfc( double x ) ;
double erfshape( double x ) ;

/* gamma.c */
double gamma(double x) ;
double lgamma(double x) ;

/* dblutil.c */
int dblequ( double a, double b ) ;
int fltequ( float a, float b ) ;
int dblroundtoint(double x) ;
double dblroundtointdouble(double x) ;
int dblisint(double x) ;
void dblswap(double *a,double *b) ;

int dblsolvequadratic(double a, double b, double c, double *result ) ;
int dblsolvecubic(  double a, double b, double c, double d, double *result ) ;

/* random.c */
int intpprand(void) ;
void intpprandsetseed(int i) ;
int intpprandgetseed(void) ;
double dblpprand(void) ;
double dblhammersley(int nbase, int i, int len) ;

/* errorh.c */
#ifndef _ERRMSG_DEFINED
#define _ERRMSG_DEFINED
struct errmsg
{
  void (**handler)(char *fmt, ... ) ;
  char *fmt ;
} ;
void errdisp( struct errmsg *errmsgs, unsigned int num, ... ) ;
void errdefwar(char *fmt, ...) ;
void errdeferr(char *fmt, ...) ;
#endif


/* spline.c */
#define NATSPLINE 1e30
int spline( double *x, double *y, double *y2, unsigned int n,
  double yp1, double ypn ) ;
int splint( double *xa, double *ya, double *y2a, unsigned int n,
  double x, double *y ) ;
int splintd( double *xa, double *ya, double *y2a, unsigned int n,
  double x, double *y, double *dy ) ;
int splintdd( double *xa, double *ya, double *y2a, unsigned int n,
  double x, double *y, double *dy, double *d2y ) ;
