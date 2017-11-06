/* sim.h - Simulation header file */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

/* This header file contains:
 * - All global data in sim.c
 * - All Function prototypes
 * - Position conventions for ODEMAN functions
 *
 * Needs:
 * odeman.h
 */

#define ODEFNC_INI 100 /* Kernel: Initilaization, No valid particles */
#define ODEFNC_USR 200 /* User  : User functions */
#define ODEFNC_INT 300 /* Kernel: Intermediate: Calculate forces */
#define ODEFNC_FOR 400 /* User  : Forces are prsent here */
#define ODEFNC_TER 500 /* Kernel: Termination, Ignored particle fields */

#define gptremoveparticle(X) (X)->tokill=1

void gptmain_sim( void *info ) ;
void gptsim_lst( void ) ;

extern struct odeinfo odeinfo ;

extern double tstart ;
extern double dtstart ;

extern double tottime ;
extern unsigned long numderivs ;
extern double simcputime ;

void gpssimerr( char *fmt, ... ) ;
void gpssimwar( char *fmt, ... ) ;
