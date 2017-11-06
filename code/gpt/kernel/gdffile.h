/* gdffile.h - declarations for GDF file structure */

#define GDFMAJ 1
#define GDFMIN 1

#include "systyp.h"

#define GDFNAMELEN 16	     /* Length of the ascii-names */
#define GDFID 94325877L

/* Data types */
#define t_ascii  0x0001      /* Ascii string	      */
#define t_s32	 0x0002      /* Signed long	      */
#define t_dbl	 0x0003      /* Double		      */

#define t_undef  0x0000      /* Data type not defined */
#define t_nul	 0x0010      /* No data 	      */
#define t_u8	 0x0020      /* Unsigned char	      */
#define t_s8	 0x0030      /* Signed char	      */
#define t_u16	 0x0040      /* Unsigned short	      */
#define t_s16	 0x0050      /* Signed short	      */
#define t_u32	 0x0060      /* Unsigned long	      */
#define t_u64	 0x0070      /* Unsigned 64bit int    */
#define t_s64	 0x0080      /* Signed 64bit int      */
#define t_flt	 0x0090      /* Float		      */

/* Block types */
#define t_dir	 0x0100      /* Directory entry start */
#define t_edir	 0x0200      /* Directory entry end   */
#define t_sval	 0x0400      /* Single valued	      */
#define t_arr	 0x0800      /* Array		      */


/* Main header in DataFile */
struct gdfmainhead
{
  U32 ID ;		     /* ID so you can see that its a GDF  */
  U32 cretime ; 	     /* Creation time			  */
  U8  creator[GDFNAMELEN] ;  /* Creator of the DataFile 	  */
  U8  destin[GDFNAMELEN] ;   /* Destination, "" means "General"   */

  U8  gdfmaj ;		     /* Major version of GDF-software	  */
  U8  gdfmin ;		     /* Minor version of GDF-software	  */

  U8  cremaj ;		     /* Major version of creator	  */
  U8  cremin ;		     /* Minor version of creator	  */
  U8  desmaj ;		     /* Major version of destination or 0 */
  U8  desmin ;		     /* Minor version of destination or 0 */
  U8  dummy1 ;		     /* Alignment (on 32-bit boundary)    */
  U8  dummy2 ;
} ;


/* Header in DataFile */
struct gdfhead
{
  U8 name[GDFNAMELEN] ;
  U32 type ;
  U32 size ;
} ;
