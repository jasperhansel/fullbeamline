/* compat.h - Compatability header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

/* gps.h */
#define SQR gptSQR

/* Output */
#define gpswarr gptoutputdoublearray
#define gpswval gptoutputdouble
#define gpswsgroup gptoutputdoublegroup
#define gpswegroup gptoutputendgroup

/* Init */
typedef void (*gpsfnc)(void *info) ;
#define gpsaddfnc gptaddmainfunction
#define GPSFNC_STP GPTMAINFNC_STP
#define GPSFNC_INI GPTMAINFNC_INI
#define GPSFNC_PAR GPTMAINFNC_PAR
#define GPSFNC_SIM GPTMAINFNC_SIM
#define GPSFNC_TER GPTMAINFNC_TER
#define GPSFNC_EXT GPTMAINFNC_EXT

/* Sim */
void gptaddparticle(gptinit *init, double *r, double *Br) ;
void gptaddparticleGB(gptinit *init, double *r, double *GBr) ;

/* Axis */
#define VECSQR gptVECSQR
#define LEN gptVECLEN
#define INP gptVECINP
#define BtoUCS gptdirectiontoUCS

/* Elem */
void readtabfile(gptinit *init, char *fn, double **px, double **py, unsigned int *pn ) ;
