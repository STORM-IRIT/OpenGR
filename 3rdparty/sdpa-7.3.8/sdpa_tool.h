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
/*--------------------------------------------------
  sdpa_tool.h
--------------------------------------------------*/

#ifndef __sdpa_tool_h__
#define __sdpa_tool_h__

#include "sdpa_right.h"  

#include <iostream>
#include <sys/time.h>
#include <string>
#include <cstring>


namespace sdpa {

// Note : only GNU has the macro __PRETTY_FUNCTION__
#define USE_PRETTY_FUNCTION 0
#if USE_PRETTY_FUNCTION
#define rMessage(message)  \
{ cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ \
  << " [ " <<  __PRETTY_FUNCTION__ << " ] "  << endl; }
#elif 1
#define rMessage(message) \
{cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ << endl; }
#else
#define rMessage(message) ;
#endif

#if USE_PRETTY_FUNCTION
#define rError(message) \
{ cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ \
  << " [ " << __PRETTY_FUNCTION__ << " ] "  << endl; \
exit(false);}
#else
#define rError(message) \
{cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ << endl; \
exit(false);}
#endif

#if 1
#define NewArray(val,type,number) \
  {val = NULL; \
    try{ val = new type[number]; } \
    catch(bad_alloc){ \
        rMessage("Memory Exhausted (bad_alloc)"); abort(); } \
    catch(...){ \
        rMessage("Fatal Error (related memory allocation"); abort(); } \
  }
#else
#define NewArray(val,type,number) \
  {rMessage("New Invoked"); \
   val = NULL; val = new type[number]; \
   if  (val==NULL) {rError("Over Memory");} \
  }
#endif

#define DeleteArray(val) \
  { if  (val!=NULL) { \
      delete[] val; \
      val = NULL; \
    } \
  }

#define REVERSE_PRIMAL_DUAL 1


// These are constant. Do NOT change
extern int IZERO   ; // =  0;
extern int IONE    ; // =  1;
extern int IMONE   ; // = -1;
extern double DZERO; // =  0.0;
extern double DONE ; // =  1.0;
extern double DMONE; // = -1.0;

class Time
{
public:
  static double rGetUseTime();
  static void   rSetTimeVal(struct timeval & targetVal);
  static double rGetRealTime(const struct timeval & start,
			     const struct timeval & end);
};

#define PROCESS_TIME 0
#if PROCESS_TIME // count time with process time
#include <sys/times.h>
#define TimeStart(START__) \
   static double START__; START__ = Time::rGetUseTime()
#define TimeEnd(END__) \
   static double END__;   END__ = Time::rGetUseTime()
#define TimeCal(START__,END__) (END__ - START__)
#else // count time with real time
#define TimeStart(START__) \
   static struct timeval START__; Time::rSetTimeVal(START__)
#define TimeEnd(END__) \
   static struct timeval END__; Time::rSetTimeVal(END__)
#define TimeCal(START__,END__) Time::rGetRealTime(START__,END__)
#endif

#define sdpa_dset(dset_length,dset_value,dset_pointer,dset_step) \
for (int dset_i=0,dset_index = 0; dset_i<dset_length; ++dset_i) { \
  dset_pointer[dset_index] = dset_value; \
  dset_index += dset_step; \
}


} // end of namespace 'sdpa'

#endif // __sdpa_tool_h__
