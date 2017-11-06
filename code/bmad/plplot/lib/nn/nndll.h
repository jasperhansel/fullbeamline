#ifndef __NN_DLL_H
#define __NN_DLL_H

#ifdef USINGDLL
  #if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
    #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
      #define NNDLLEXPORT    __declspec( dllexport )
      #define NNDLLIMPORT    __declspec( dllimport )
    #else
      #define NNDLLEXPORT
      #define NNDLLIMPORT
    #endif
  #elif defined ( __CYGWIN__ )
    #define NNDLLEXPORT    __declspec( dllexport )
    #define NNDLLIMPORT    __declspec( dllimport )
  #elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// NNDLLEXPORT to be visible.
    #define NNDLLEXPORT    __attribute__ ( ( visibility( "default" ) ) )
    #define NNDLLIMPORT
  #endif
#endif

// For an unknown compiler or static built we clear the macros
#ifndef NNDLLEXPORT
  #define NNDLLEXPORT
  #define NNDLLIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case). If
// cmake sets the corresponding macro xxxx_EXPORTS if the
// corresponding library is built DLLIMPEXP is set to DLLEXPORT
#if defined ( csironn_EXPORTS )
  #define NNDLLIMPEXP    NNDLLEXPORT
  #define NNDLLIMPEXP_DATA( type )    NNDLLEXPORT type
#else
  #define NNDLLIMPEXP    NNDLLIMPORT
  #define NNDLLIMPEXP_DATA( type )    NNDLLIMPORT type
#endif

#endif // __NN_DLL_H
