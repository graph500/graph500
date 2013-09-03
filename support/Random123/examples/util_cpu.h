/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef UTIL_CPU_H__
#define UTIL_CPU_H__ 1

#include "util.h"

/*
 * burn a little CPU by computing a logistic map function
 * before reading /proc/cpuinfo, might fool
 * energy-saving CPU into showing its true speed
 */
static double warmupCPU(long n){
    double d = 0.3;
    int i;
    for (i = 0; i < n; i++) {
	d = 3.6 * d * (1. - d);
    }
    dprintf(("logistic map produced %f\n", d));
    return d;
}

#if defined(_MSC_VER)
#define NOMINMAX    /* tells Windows.h to NOT define min() & max() */
#include <Windows.h>
static double clockspeedHz(int *ncores, char **modelnamep) {
    /*
     * To-Do: one could read clock speed from registry.  Or
     * maybe WMI?
     * http://www.codeproject.com/KB/system/Processor_Speed.aspx
     * http://msdn.microsoft.com/en-us/library/aa394373(v=vs.85).aspx
     */
    if (ncores) {
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	*ncores = sysinfo.dwNumberOfProcessors;
    }
    if (modelnamep) {
	/*
	 * To-Do: use __cpuid to get brand string see
	 * http://msdn.microsoft.com/en-us/library/hskdteyh%28v=vs.80%29.aspx
	 */
	*modelnamep = ntcsdup("Windows");
   }
   return 0.;
}
#elif defined(__APPLE__)
static double clockspeedHz(int *ncores, char **modelnamep){
    FILE *fp = popen("sysctl hw.cpufrequency", "r");
    double hz = 0.;
    double d = warmupCPU(100L*1000L*1000L);
    if( fscanf(fp, "%*s %lf", &hz)  != 1 )
 	return 0.;
    pclose(fp);
    if(ncores) *ncores = 1;
    if(modelnamep) *modelnamep = ntcsdup("Apple");
    return hz;
}
#elif defined(__SUNPRO_CC) || defined(__SUNPRO_C) || (defined(__GNUC__)&&defined(__sun__))
static double clockspeedHz(int *ncores, char  **modelnamep){
    FILE *fp = popen("kstat -p -s current_clock_Hz", "r");
    double hz = 0.;
    double d = warmupCPU(100L*1000L*1000L);
    /* To-do: get a model name from kstat too */
    if(modelnamep) *modelnamep = ntcsdup("Solaris");
    int nc = 0;
    while( fscanf(fp, "%*s %lf", &hz) == 1 ){
        nc++;
    }
    if(ncores) *ncores = nc;
    return hz;
}
#elif defined(__linux__)
/* Read the clock speed from /proc/cpuinfo - Linux-specific! */
static double clockspeedHz(int *ncores, char **modelnamep){
    char *s, buf[1024]; /* long enough for any /proc/cpuinfo line */
    double Mhz = 0.;
    double xMhz;
    int i;
    double d = warmupCPU(100L*1000L*1000L);
    FILE *fp;
    if ((fp = fopen("/proc/cpuinfo", "r")) == NULL) {
	if (ncores) *ncores = 1;
	if (modelnamep) *modelnamep = ntcsdup("unknown");
	return 0.;
    }
    if (ncores) *ncores = 0;
    while (fgets(buf, sizeof buf, fp) != NULL) {
	if (modelnamep && (s = strstr(buf, "model name")) != NULL) {
	    CHECKNOTZERO(s = strchr(s, ':'));
	    while (*++s == ' ')
		;
	    i = strchr(s, '\n') - s;
	    *modelnamep = (char *)malloc(i + 1);
	    memcpy(*modelnamep, s, i);
	    (*modelnamep)[i] = '\0';
	    dprintf(("raw modelname is %d bytes: %s\n", i, *modelnamep));
	    nameclean(*modelnamep);
	    dprintf(("cleaned modelname is %s\n", *modelnamep));
	}
	if ((s = strstr(buf, "cpu MHz")) || (s = strstr(buf, "clock"))) {
	    if (s[1] == 'p') // cpu MHz
		CHECKNOTZERO(sscanf(s, "cpu MHz : %lf %n", &xMhz, &i));
	    else // clock
		CHECKNOTZERO(sscanf(s, "clock : %lfMHz %n", &xMhz, &i));
	    dprintf(("parsed %f %d\n", xMhz, i));
	    if (xMhz > Mhz) Mhz = xMhz;
	    s += i;
	    if (ncores) *ncores += 1;
	}
    }
    d = Mhz*1e6;
    dprintf(("clockspeed is %f\n", d));
    return d;
}
#elif defined(__FreeBSD__)
static double clockspeedHz(int *nnodes, char **modelnamep){
    /* Seems to work with FreeBSD 8.2. */
    FILE *fp = popen("sysctl hw.ncpu hw.clockrate hw.model", "r");
    int ncpu, clockrate;
    if( fscanf(fp, "%*s%d%*s%d%*s", &ncpu, &clockrate) != 2 )
	return 0;
    if(nnodes) *nnodes = ncpu;
    if(modelnamep){
	char buf[256];
	if(fgets(buf, sizeof(buf), fp) == NULL) 
	    *modelnamep = ntcsdup("error reading sysctl");
	else
	    *modelnamep = ntcsdup(buf);
    }
    return 1.e6*clockrate;
}
#else
static double clockspeedHz(int *nnodes, char **modelnamep){
    if(nnodes) *nnodes = 1;
    if(modelnamep) *modelnamep = ntcsdup("unknown");
    return 0.;
}
#endif /* ! _MSC_VER */

#define uint unsigned int
typedef struct test_info {
    double hz;
    int ncores;
    char *cpuname;
} CPUInfo;

CPUInfo *cpu_init(const char *arg)
{
    CPUInfo *tp;
    tp = (CPUInfo*)malloc(sizeof(CPUInfo));
    tp->hz = clockspeedHz(&tp->ncores, &tp->cpuname);
    printf("%d cores, %.3f Ghz, cpu %s\n", tp->ncores, tp->hz*1e-9, tp->cpuname);
    if (arg) {
	int n = atoi(arg);
	if (n) {
	    printf("setting cores to %d\n", n);
	    tp->ncores = n;
	}
    }
    return tp;
}

void cpu_done(CPUInfo *tp)
{
    free(tp->cpuname);
    free(tp);
}
#endif /* UTIL_CPU_H__ */
