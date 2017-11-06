/* parse.h */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

struct symbol
{
  char *name ;
  int type ; 
  union
  {
    double val ;
    int num ;
    double (*ptr)(double val) ;
    gptscatterfnc scatter ;
  } u ;
  void *info ;
  struct symbol *next ;
} ;

struct symbol *gptlookup( char *name ) ;
struct symbol *gptinstall( char *name, int t ) ;
struct symbol *gptinstalldouble( char *name, double d ) ;
struct symbol *gptinstalldoublefnc( char *name, double (*fnc)(double) ) ;
struct symbol *gptinstallscatterfnc( char *name, gptscatterfnc scatter, void *info ) ;

int gptgetvardouble(char *name, double *val) ;
int gptgetvarscatterfnc(char *name, gptscatterfnc *scatter, void **info ) ;
