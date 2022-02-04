/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004-2012 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#include <mex.h>


#define BUF_LENGTH 4096

int fprintf(FILE *stream, const char *format, ...)
{
  va_list arg;
  va_start(arg,format);
  int return_size = 0;
  if (stream == stdout) {
    char tmp_buf[BUF_LENGTH];
    vsprintf(tmp_buf,format,arg);
    #if 0
    printf("tmp_buf = %s",tmp_buf);
    #endif
    mexPrintf(tmp_buf);
    mexEvalString("drawnow;"); /* to dump string.*/
    return_size = strlen(tmp_buf);
    if ( return_size >= BUF_LENGTH ) {
      mexPrintf("Too Long Message To PrintOut "
		"(some part might be truncated)");
    }
  }
  else {
    return_size = vfprintf(stream,format,arg);
  }
  va_end(arg);
  return return_size;
}

static int internal_fprintf(FILE *stream, const char *format, ...)
{
  va_list arg;
  va_start(arg,format);
  int return_size = 0;
  return_size = vfprintf(stream,format,arg);
  va_end(arg);
  return return_size;
}

size_t fwrite(const void *ptr, size_t size, size_t nmemb,
              FILE *stream)
{
  #if 1
  fprintf(stream,"%s",(const char*)ptr);
  #else
  // internal_fprintf does not work well here
  internal_fprintf(stream, (const char*) internal_fprintf);
  #endif
  return 1;
}

int fputc(int c, FILE* fp)
{
  if (fp == stdout) {
    mexPrintf("%c",(unsigned char)c);
  }
  else {
    internal_fprintf(fp,"%c",c);
    /*    fprintf(fp,"%c ",c);*/
  }
  return c;
}

#ifdef __GNUC__
/*Only GNU, to avoid
      warning: 'noreturn' function does return
*/
static void internal_exit()           __attribute__ ((noreturn));
extern void mexErrMsgTxt(const char*) __attribute__ ((noreturn));
       void exit()                    __attribute__ ((noreturn));
       void abort()                   __attribute__ ((noreturn));
#endif

static void internal_exit()
{
  mexWarnMsgTxt("SDPA exits with some error.");
  mexWarnMsgTxt("Matlab should be reboot to clear up memory space.");
  mexErrMsgTxt("SDPA exits with some error.");
}

void exit(int status)
{
  internal_exit();
}
void abort(void)
{
  internal_exit();
}
