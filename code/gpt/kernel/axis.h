/* axis.h - Axis header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

/* This header file contains:
 * - All global data in axis.c
 * - Axis data structure and function prototypes
 */

void gptbuildECS(gptinit *init) ;

struct axis *getaxis( char *axisname ) ;

#define GPTELEM_LOCAL  1
#define GPTELEM_GLOBAL 2
void gptaddEBelementU(gptinit *init, simfn simfn, exitfn exitfn, int type, void *info ) ;
#define gptaddEBelement(a,b,c,d,e) gptaddEBelementU(a,(simfn)(b),(exitfn)(c),d,e)
void gptaddboundaryelementU( gptinit *init, gptboundaryfnc boundfnc, exitfn exitfnc, int type, void *info ) ;
#define gptaddboundaryelement(a,b,c,d,e) gptaddboundaryelementU(a,(gptboundaryfnc)b,(exitfn)c,d,e) ;

void gptmain_axis(void *info) ;

void createidentitytransform( gpttransform *t ) ;
void gptconcattransform( gpttransform *t, gpttransform *t1, gpttransform *t2 ) ;
void createdirectionaltransform(struct gpttransform *t,double *p) ;

void gpttoWCS( gpttransform *t, double *rin, double *rout ) ;
void gptdirectiontoWCS( gpttransform *t, double *Bin, double *Bout ) ;
void gptadddirectiontoWCS( gpttransform *t, double *Bin, double *Bout ) ;
void gpttoUCS( gpttransform *t, double *rin, double *rout ) ;
void gptdirectiontoUCS( gpttransform *t, double *Bin, double *Bout ) ;

void gptr2carth(double Fr,double x,double y,double *Fx,double *Fy) ;
void gptrphi2carth(double Fr,double Fphi,double x,double y,double *Fx,double *Fy) ;
