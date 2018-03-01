/*
 * $Source: /psa/share/repository/slide/src/utils/inc/slide_itimer.h,v $
 * $Revision: 1.1 $ 
 * $Author: vanvoor4 $
 * $Date: 2009/05/04 15:03:01 $

 * $Log: slide_itimer.h,v $
 * Revision 1.1  2009/05/04 15:03:01  vanvoor4
 * adding utils
 *
 * Revision 1.3  2007/09/28 18:33:49  toneroma
 * *** empty log message ***
 *
 * Revision 1.2  2006/09/19 18:53:22  vanvoor4
 * Added support for other timers besides profile and fixed handling of wrong
 * signal.
 *
 * Revision 1.1  2006/09/15 16:00:52  vanvoor4
 * Added support for profile timer
 *
 *
 */

#ifndef _TIMER_H_INCLUDED
#define _TIMER_H_INCLUDED
#include <time.h>

int slide_itimer_start(int interval_sec_in);

int slide_itimer_get(double *real, double *virt, double *prof);

char* slide_get_local_time(char *buf, size_t buflen);

/*! Wrapper to eliminate buggy gcc warning about using "%c" in strftime.*/
/*! From the man page for strftime:
 * Some buggy versions of gcc complain about the use of %c: warning: `%c'
 * yields only last 2 digits of year in some locales.  Of course
 * programmers are encouraged to use %c, it gives the preferred date
 * and time representation. One meets all kinds of strange obfuscations to
 * circumvent this gcc problem. A relatively clean one is to add this
 * intermediate function.
 *
 * @param s Character array of size max
 * @param max Size of the array s
 * @param fmt Format string
 * @param tm The tm structure to convert to a string.
 * @return The number of chars written to s if the time fits, else 0.
 */
size_t my_strftime(char *s, size_t max, const char  *fmt, const struct tm *tm);

#endif
