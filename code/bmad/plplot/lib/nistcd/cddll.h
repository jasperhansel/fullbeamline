#ifndef __CD_DLL_H
#define __CD_DLL_H

#ifdef USINGDLL
  #if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
    #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
      #define CDDLLEXPORT    __declspec( dllexport )
      #define CDDLLIMPORT    __declspec( dllimport )
    #else
      #define CDDLLEXPORT
      #define CDDLLIMPORT
    #endif
  #elif defined ( __CYGWIN__ )
    #define CDDLLEXPORT    __declspec( dllexport )
    #define CDDLLIMPORT    __declspec( dllimport )
  #elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// CDDLLEXPORT to be visible.
    #define CDDLLEXPORT    __attribute__ ( ( visibility( "default" ) ) )
    #define CDDLLIMPORT
  #endif
#endif

// For an unknown compiler or static built we clear the macros
#ifndef CDDLLEXPORT
  #define CDDLLEXPORT
  #define CDDLLIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case). If
// cmake sets the corresponding macro xxxx_EXPORTS if the
// corresponding library is built DLLIMPEXP is set to DLLEXPORT
#if defined ( nistcd_EXPORTS )
  #define CDDLLIMPEXP    CDDLLEXPORT
  #define CDDLLIMPEXP_DATA( type )    CDDLLEXPORT type
#else
  #define CDDLLIMPEXP    CDDLLIMPORT
  #define CDDLLIMPEXP_DATA( type )    CDDLLIMPORT type
#endif

#endif // __CD_DLL_H
