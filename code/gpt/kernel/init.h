/* init.h - Initialisation header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#define LINELENGTH 256

extern char *globalfn ;
extern unsigned int globalln ;
extern FILE *fpin ;

#define GPTMAINFNC_STP 100  /* Kernel: Commandline, inputfile etc. */
#define GPTMAINFNC_INI 200  /* User  : Particles can be added to the sets */
#define GPTMAINFNC_LST 250  /* Kernel: Transform sets to master particle list */
#define GPTMAINFNC_DYN 275  /* Kernel: Start dynamic memorymanagement for all particle related components */
#define GPTMAINFNC_PAR 300  /* User  : All particles must be initialized and are stored in the master list */
#define GPTMAINFNC_SIM 400  /* Kernel: Simulation engine */
#define GPTMAINFNC_TER 500  /* User  : Termination routines */
#define GPTMAINFNC_EXT 600  /* Kernel: Close outputfile and general cleanup */
void gptaddmainfunction( int position, void (*func)(void *info), void *info ) ;

/* Element parameters */
#define GPTTYPE_INT    1
#define GPTTYPE_DOUBLE 2
#define GPTTYPE_STRING 3

char  *gptgetname(gptinit *init) ;
int    gptgetargnum(gptinit *init) ;
int    gptgetargint(gptinit *init, int argnum) ;
double gptgetargdouble(gptinit *init, int argnum) ;
char  *gptgetargstring(gptinit *init, int argnum) ;
int    gptgetargtype(gptinit *init, int argnum) ;

/* Error handlers */
void gpsinperr( char *fmt, ... ) ;
void gpsinpwar( char *fmt, ... ) ;
void gpsclnerr( char *fmt, ... ) ;
void gpsclnwar( char *fmt, ... ) ;
