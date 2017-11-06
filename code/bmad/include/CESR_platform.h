#ifndef CESR_platform_H
#define CESR_platform_H

/************************************************************************/
/* File         :  CESR_platform.h                                      */
/*                                                                      */
/* Description  :  This header file provides system/platform            */
/*                 definitions.                                         */
/*                                                                      */
/* Author       :  M. Palmer  04/25/01                                  */
/*                                                                      */
/* Modifications:  Set up standard platform definitions.  MAP  4/30/01  */
/*                 Add generic UNIX specification.        MAP  9/20/01  */
/*                 Standardize for likely CESR platforms. MAP  9/29/01  */
/*                                                                      */
/************************************************************************/

/* UNIX */
#if defined(unix) || defined(__unix__)
#ifndef CESR_UNIX
#define CESR_UNIX
#endif
#endif

/* LINUX */
#if defined(__linux__) || defined(__linux) || defined(linux)
#ifndef CESR_LINUX
#define CESR_LINUX
#endif
#ifndef CESR_UNIX
#define CESR_UNIX
#endif
#endif

/* WINDOWS */
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(WINNT)
#ifndef CESR_WINCVF
#define CESR_WINCVF
#endif
#endif

/* ANALOG DEVICES SHARC DSP */
#if defined(__2106x__)
#ifndef SHARC_DSP
#define SHARC_DSP
#endif
#endif

#endif






