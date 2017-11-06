/* scatter.h - Scatter model header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

double gpsparserr_trajectory( double t, double dt, double *xstart, double *xend, double *xerr, void *info) ;

struct scatter_parinfo
{
  double Qin, Qout ;
  double Ein, Eout ;
  double r[3] ;
  double inp ;
} ;

struct scatter_info
{
  /* Coordinate transform */
  gpttransform tf ;

  /* Memory management for scatter statistics, grand=exit, total=tout */
  struct scatter_parinfo *grandpars, *totalpars, *pars ;
  int numgrandpars, numtotalpars, numpars ;
} ;

/* Forward declarations */
void gptscatterinit(gptinit *init, struct scatter_info *scatinfo, char *name) ;
void gptscatterinitpar(struct scatter_info *scatinfo) ;
void gptscatteraddparmqn(struct scatter_info *scatinfo,gptparset *ppp, double *Wr, double *WGBr, double m, double q, double n) ;
void gptscatterremoveparticle(struct scatter_info *scatinfo,gpttrajectory *traj, gptpar *par) ;
