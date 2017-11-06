/* odeman.h - ODE Management header file*/

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

/* ODEMAN List numbers */
#define ODE_TOTALNUM 6
#define ODE_INI  0  /* Initialization before a timestep */
#define ODE_FPR  1  /* Calculate derivates */
#define ODE_OUT  2  /* Possibility of output after the first derivates call */
#define ODE_ERR  3  /* Calculate error */
#define ODE_END  4  /* Called after a successfull RK step */
#define ODE_FAIL 5  /* Called after an un successfull RK step */

/* ODEMAN corresponding function types */
typedef void   (* odeinifnc)(double t, double **px, void *info) ;
typedef void   (* odefprfnc)(double t, double *x, double *p, void *info) ;
typedef void (* odemtfprfnc)(double t, double *x, double *p, void *info, int nOff, int nCPU) ;
typedef int    (* odeoutfnc)(double t, double *dt, double *x, void *info) ;
typedef double (* odeerrfnc)(double t, double dt, double *xstart, double *xend, double *xerr, void *info) ;
typedef void   (* odeendfnc)(double tstart, double tend, double *dt, double *xstart, double *xend, void *info, void *stepinfo) ;
typedef void  (* odefailfnc)(double t, double dt, double *xstart, double *xend, void *info) ;

/* Basic ODEMAN structure */
struct odeinfo
{
/* Initialized by odeinit */
  struct gptsyncset *syncset ;
  struct gptsyncarray *sax ;

/* Old variables, maintained for compatibility */
  double *x ;
  int num ;

/* Initialized by user */
  double tstart ;
  double tend ;
  double dtstart ;
  double dtmin ;
  double dtmax ;

/* Initialized by user or odeman */
  odeinifnc stepini ;
  void  (*stepfpr )( int nOff, int nCPU, double tstart, double tend, double t, double *x, double *p, void *info ) ;
  odeoutfnc  stepout ;
  odeerrfnc  steperr ;
  odeendfnc  stepend ;
  odefailfnc stepfail ;
} ;

/* ODEMAN functions */
void odeinit( struct odeinfo *odeinfo ) ;
int odeaddvar( struct odeinfo *odeinfo, int num, int *offset ) ;
void odeaddinifunction( int position, odeinifnc func, void *info ) ;
void odeaddfprfunction( int position, odefprfnc func, void *info ) ;
void odemtaddfprfunction( int position, odemtfprfnc func, void *info ) ;
void odeaddoutfunction( int position, odeoutfnc func, void *info ) ;
void odeadderrfunction( int position, odeerrfnc func, void *info ) ;
void odeaddendfunction( int position, odeendfnc func, void *info ) ;
void odeaddfailfunction( int position, odefailfnc func, void *info ) ;
int odedrv( struct odeinfo *odeinfo ) ;

/* List of stepper functions */
int rkdrv( struct odeinfo *rkinfo, void *locinfo ) ;
double *rkdrvinterpolate(double tstart, double tend, double t, void *stepinfo) ;
void rkdrvinterpolate(double tstart, double tend, double t, void *stepinfo, int num, int offset, double *out, double *doutdt ) ;


#define RKERR_NEM 1	   /* Not enough memory   */
#define RKERR_ISS 2        /* Step size too small */

/* For internal use */
extern double rk_minpow ;
extern double rk_maxpow ;
extern double rk_safety ;
extern double rk_minmul ;
extern double rk_maxmul ;

extern int mt_minomp ;
extern int mt_maxomp ;
extern unsigned int gptCPUs ;

/* Duplicates from rk.h */
#define ODEERR_NEM  1  /* Not enough memory */
#define ODEERR_ISS  2  /* Invalid (to small) step size */

/* ODEMAN error values */
#define ODEERR_IFN 11  /* Invalid function number (odeaddfnc) */
#define ODEERR_INI 12  /* Invalid initialization (no fpr or err specified) */
