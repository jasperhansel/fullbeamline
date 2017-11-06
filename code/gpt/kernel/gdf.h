/* gdf.h - General DataFile Format: Header file */

/* This file is part of GPS and belongs to:
 * S.B. van der Geer, M.J. de Loos and Roots Software
 */

#include <stdio.h>
#include <stdlib.h>
#include "gdffile.h"

#ifdef _MSC_VER
typedef __int64 off64_t ;
#endif

#ifdef __APPLE__
typedef off_t off64_t ;
#endif

/* File descriptor for GDF[rw] functions */
struct gdff
{
  FILE *fp ;
  const char *fn ;
  int append ;		     /* Only for gdfwmainhead */
} ;


/* Memory representation of DataBlock */
struct gdfdisk
{
  U32 type ;
  U32 size ;

  struct gdff *gdff ;
  off64_t foff ;
} ;

struct gdfdata
{
  unsigned int numchild ;
  struct gdfdata *next ;
  struct gdfdata *childs ;
  struct gdfdata *parent ;   /* NULL for virtual root directory */
  struct gdfhead oh ; /* Current information in memory */
  struct gdfdisk di ; /* Stored information on disk */
  void *buf ;
} ;


/* Memory representation of DataFile */
struct gdfmem
{
  struct gdfmainhead mh ;
  struct gdfdata ds ;	     /* Virtual root directory */
} ;


/* Write functions */
void gdfwmainhead( struct gdff *gdff, const char *creator, U8 cremaj, U8 cremin, const char *destin , U8 desmaj, U8 desmin ) ;
void gdfw( struct gdff *gdff, struct gdfhead *h, const void *buf ) ;
void gdfwascii( struct gdff *gdff, const char *name, U32 type, const char *string ) ;
void gdfwsval( struct gdff *gdff, const char *name, U32 type, const void *data ) ;
void gdfwarr( struct gdff *gdff, const char *name, unsigned int LEN, U32 type, const void *data ) ;
void gdfwclose( struct gdff *gdff ) ;
void gdfwmem( struct gdff *gdff, struct gdfmem *gm ) ;


/* Read functions */
#define GDFR_ABORTONERROR  0x1
#define GDFR_READARRAYS    0x2
#define GDFR_NOREADARRAYS  0x4
int gdfrgetfileMBs(struct gdff *gdff) ;
int gdfrmem( struct gdff *gdff, struct gdfmem *gm, int flags ) ;
void gdfrdata(struct gdfdata *ds) ;
void gdfrclose(struct gdff *gdff) ;


/* Memory functions */
struct gdfdata *gdfm( struct gdfdata *ds, struct gdfhead *h, const void *buf ) ;
struct gdfdata *gdfmascii( struct gdfdata *ds, const char *name, U32 type, const char *string ) ;
struct gdfdata *gdfmsval( struct gdfdata *ds, const char *name, U32 type, const void *data ) ;
struct gdfdata *gdfmarr( struct gdfdata *ds, const char *name, unsigned int LEN, U32 type, const void *data ) ;
void gdfmfillds( struct gdfdata *ds ) ;
void gdfmfillhead( struct gdfmainhead *mh,
		   const char *creator, U8 cremaj, U8 cremin,
		   const char *destin , U8 desmaj, U8 desmin ) ;
void gdfmfreechilds(struct gdfdata *ds) ;
void gdfmforall( struct gdfdata *ds,
		 void (*fnc)(struct gdfdata *ds, unsigned int level ),
		 unsigned int level ) ;
void gdfmforalldir( struct gdfdata *ds,
		    void (*fnc)(struct gdfdata *ds, unsigned int level ),
		    unsigned int level ) ;
void gdfmforalldirclear( struct gdfdata *ds,
		    void (*fnc)(struct gdfdata *ds, unsigned int level ),
		    unsigned int level ) ;
void gdfmforalldirinfo( struct gdfdata *ds,
		    void (*fnc)(struct gdfdata *ds, unsigned int level, void *info),
		    unsigned int level, void *info ) ;
void gdfmforalldirinfoclear( struct gdfdata *ds,
		    void (*fnc)(struct gdfdata *ds, unsigned int level, void *info),
		    unsigned int level, void *info ) ;
void gdfmcleardata( struct gdfdata *ds ) ;
void gdfmcleardatas( struct gdfdata *ds ) ;
struct gdfdata *gdfmmkdir( struct gdfdata *ds, struct gdfdata *dirds ) ;
struct gdfdata *gdfmschild( struct gdfdata *ds, const char *child ) ;
struct gdfdata *gdfmschildval( struct gdfdata *ds, struct gdfdata *dscmp ) ;
struct gdfdata *gdfmgetdata( const char *name, struct gdfdata *ds ) ;
double *gdfmgetdbl(struct gdfdata *ds, const char *name) ;
double *gdfmgetarrdbl(struct gdfdata *ds, const char *name, int *len) ;
double *gdfmgetarrdbllen(struct gdfdata *ds, const char *name, int points) ;

int gdfmhasarr(struct gdfdata *ds ) ;
struct gdfdata *gdfmnext(struct gdfdata *ds) ;
struct gdfdata *gdfmnextarr(struct gdfdata *ds) ;
struct gdfdata *gdfmprev(struct gdfdata *ds) ;
struct gdfdata *gdfmprevarr(struct gdfdata *ds) ;
char *gdfmmakeascii( U32 type, void *buf, int width, int verbose ) ;
char *gdfmgettitle(char *buf,struct gdfdata *ds) ;
void gdfmstorelocation(struct gdfdata *pointerds, struct gdfdata *outds) ;
void gdfmrestorelocation(struct gdfdata **pointerds, struct gdfdata *inds, struct gdfdata *rootds) ;

int gdfmtodbl( struct gdfdata *ds ) ;
void gdfmclearnotread( struct gdfdata *ds) ;


/* Support/statistics functions */
#ifdef _MSC_VER
off64_t ftello( FILE *stream ) ;
int fseeko( FILE *stream, off64_t offset, int origin );
#endif
S32 gdfsize( U32 type ) ;
void gdfswinit( const char *gdfoutfo, const char *gdfoutfa, struct gdff *gdff ) ;
void gdfsrinit( const char *gdfinf, struct gdff *gdff ) ;

struct gdfsimplestats
{
	// Array name
	char *array ;

	// Direct statistics
	int num ;
	double min ;
	double max ;
	double sum ;
	double sum2 ;

	// Derived statistics
	double avg ;
	double std ;
} ;

int gdfcalcsimplestats(struct gdfsimplestats *stats,struct gdfdata *ds,const char *array,int recurse) ;
int gdfstatsmin(struct gdfdata *ds,const char *array,double *result,int recurse) ;
int gdfstatsmax(struct gdfdata *ds,const char *array,double *result,int recurse) ;
int gdfstatsmaxmin(struct gdfdata *ds,const char *array,double *result,int recurse) ;
int gdfstatsavg(struct gdfdata *ds,const char *array,double *result,int recurse) ;
int gdfstatsstd(struct gdfdata *ds,const char *array,double *result,int recurse) ;
int gdfstatssum(struct gdfdata *ds,const char *array,double *result,int recurse) ;
int gdfstatsnum(struct gdfdata *ds,const char *array,double *result,int recurse) ;

extern struct gdfstatsfuncs
{
	const char *szName ;
	int (*func)(struct gdfdata *ds,const char *array,double *result,int recurse) ;
} statsfunctable[] ;

int gdfcalcstats(struct gdfdata *ds,const char *func,const char *array,double *result,int recurse) ;


/* Error handling */
#ifndef _ERRMSG_DEFINED
#include "../../utils/utils.h"
#endif
extern struct errmsg gdfem[] ;
extern void (*agdferr)( char *fmt, ... ) ;
extern void (*agdfwar)( char *fmt, ... ) ;

#include "gdferr.h"
