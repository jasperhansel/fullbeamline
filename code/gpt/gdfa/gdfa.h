/* gdfa.h */

#define _CRT_SECURE_NO_DEPRECATE

typedef int (*gdfaprog)( double *result ) ;

struct prognames
{
  char *name ;
  gdfaprog prog ;
} ;

int    gdfmgetarr( const char *name, double **result, int *len ) ;
int    gdfmgetval( const char *name, double *result ) ;

double gdfamean(const double *n, const double *a, int num) ;
double gdfamean2(const double *n, const double *a, const double *b, int num) ;
double gdfastd(const double *n, const double *a, int num) ;
void   gdfasubavg(const double *n, double *dest, const double *src, int num) ;
int    gdfacmpdouble(const void *val1, const void *val2) ;

double stdsqrt( double x ) ;

/* Constants from gps.h */
#define gpt_c  299792458.0
#define gpt_e  2.71828182845904523536
#define gpt_pi 3.14159265358979323846
#define gpt_deg (180.0/gpt_pi)

#define gpt_ma  1.66057e-27
#define gpt_me  9.10953e-31
#define gpt_mp  1.67265e-27
#define gpt_qe -1.6021892e-19

#define gpt_mu0  (1e-7*4.0*gpt_pi)
#define gpt_eps0 (1/(gpt_c*gpt_c*gpt_mu0))
