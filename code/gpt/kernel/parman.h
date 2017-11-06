/* GPT: Particle manager */

typedef struct gptinitpar /* ap */
{
  /* Basic particle coordinates */
  double Wr[3] ;
  double GBr[3] ;
  double m, q, n, r ;
  double tstart ; /* Start time */

  /* Additional stuff */
  struct axis *paxis ;
  int ID ;
} gptinitpar ;

/* Particle set */
typedef struct gptparset /* pp */
{
  char *name ;
  gptinitpar *pars ;
  int len ;

  int ndisthammersley ;

  struct gptparset *next ;
} gptparset ;

/* Particle set functions */
gptparset *gptcreateparset(char *name) ;
void gptremoveparset(gptparset *ppp) ;
gptparset *gpttestparset(char *name) ;
gptparset *gptgetparset(char *name) ;
gptinitpar *gptgetparsetpars(gptparset *ppp, int *len) ;
gptinitpar *gptaddparmqnartid( gptparset *ppp, double *Wr, double *WGBr, double m, double q, double n, struct axis *paxis, double r, double tstart, int ID ) ;
gptinitpar *gptaddparmqnar( gptparset *ppp, double *Wr, double *WGBr, double m, double q, double n, struct axis *paxis, double r ) ;
gptinitpar *gptaddparmqn( gptparset *ppp, double *Wr, double *WGBr, double m, double q, double n ) ;
gptinitpar *gptaddpar( gptparset *ppp, double *Wr, double *WGBr ) ;
gptinitpar *gptaddparmq( gptparset *ppp, double *Wr, double *WGBr, double m, double q ) ;
gptinitpar *gptremovepar(gptparset *ppp, int n, int *len) ;

/* Simulation set variables and functions, mostly for internal use */
extern struct gptsyncset *gptparsyncset ;
void gptsimaddparmqnartid( gptparset *ppp, double *Wr, double *WGBr, double m, double q, double n, struct axis *paxis, double r, double tstart, int ID) ;
gptpar *gptgetsimpars(int *len) ;
void gptmoveparsets(double t) ;
void gptclearparsets(void) ;
void gptmain_lst( void *info ) ;
