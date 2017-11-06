/* dists.h - Distributions header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

enum gptdistdim
{
  gptdist_1D=100,
  gptdist_2DR,
  gptdist_2DSin,
} ;

double *gptgetdistribution(gptinit *init, int len, int startparam, enum gptdistdim distdim,gptparset *set) ;

/* inverse integrated distributions */
double iiUniform(double y, double xc, double sigma) ;
double iiLinear(double y, double xc, double sigma, double hstart, double hend ) ;
double iiQuadratic(double y, double xc, double sigma ) ;
double iiGaussian(double y, double xc, double sigma, double sleft, double sright ) ;
double iiCosine(double y, double xc, double sigma, double alpha ) ;
double iiRUniform(double y, double rc, double sigma ) ;
double iiRLinear(double y, double rc, double sigma, double hstart, double hend ) ;
double iiRGaussian(double y, double rc, double sigma, double sleft, double sright ) ;
double iiRSphere(double y, double rc, double sigma ) ;
double iiSUniform(double y, double thetac, double sigma ) ;
