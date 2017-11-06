#ifndef __CSA_DLL_H
#define __CSA_DLL_H

#ifdef USINGDLL
  #if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
    #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
      #define CSADLLEXPORT    __declspec( dllexport )
      #define CSADLLIMPORT    __declspec( dllimport )
    #else
      #define CSADLLEXPORT
      #define CSADLLIMPORT
    #endif
  #elif defined ( __CYGWIN__ )
    #define CSADLLEXPORT    __declspec( dllexport )
    #define CSADLLIMPORT    __declspec( dllimport )
  #elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// CSADLLEXPORT to be visible.
    #define CSADLLEXPORT    __attribute__ ( ( visibility( "default" ) ) )
    #define CSADLLIMPORT
  #endif
#endif

// For an unknown compiler or static built we clear the macros
#ifndef CSADLLEXPORT
#    define CSADLLEXPORT
#    define CSADLLIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case). If
// cmake sets the corresponding macro xxxx_EXPORTS if the
// corresponding library is built DLLIMPEXP is set to DLLEXPORT
#if defined ( csirocsa_EXPORTS )
  #define CSADLLIMPEXP    CSADLLEXPORT
  #define CSADLLIMPEXP_DATA( type )    CSADLLEXPORT type
#else
  #define CSADLLIMPEXP    CSADLLIMPORT
  #define CSADLLIMPEXP_DATA( type )    CSADLLIMPORT type
#endif

#endif // __CSA_DLL_H
